import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import plotly.express as px
import plotly.graph_objects as go
import snf
import os
from scipy.sparse import csgraph
from pyvis.network import Network
from colormap import rgb2hex
from sklearn.preprocessing import StandardScaler, MinMaxScaler
from sklearn.decomposition import PCA
from sklearn.cluster import SpectralClustering

from scipy.stats import mannwhitneyu, kruskal
from sklearn.metrics import adjusted_mutual_info_score
from statsmodels.sandbox.stats.multicomp import multipletests

from src.modules.py_deseq import py_DESeq2

import rpy2.robjects as ro
from rpy2.robjects.packages import importr
from rpy2.robjects import r, pandas2ri

from statannotations.Annotator import Annotator
from skbio.diversity import alpha_diversity
from itertools import combinations
import scikit_posthocs as sp
from scipy.stats import fisher_exact


random_state = 42

def eigengap_visual(affinity, res_dir):

    """
    Calculate the Laplacian matrix from the affinity matrix provided. The eigenvalues of the
    Laplacian are then calculated and plotted in ascending fashion.

    Parameters
    ----------
    affinity : np.array
        The similarity matrix for which to calculate the eigenvalues.
    res_dir : str
        The directory path where to save the generated plots.
    """

    laplacian = csgraph.laplacian(affinity, normed=True)
    
    # perform eigendecomposition and find eigengap
    eigenvalues = np.sort(np.linalg.eigvals(laplacian))
    plt.figure(figsize=(6, 4))
    plt.title("")

    plt.scatter(np.arange(1,16), eigenvalues[:15], s=14, color='black')
    # Adding annotations for the 4th and 6th eigenvalues with dashed lines of different colors
    #plt.axhline(y=eigenvalues[3]+0.02, linestyle='--', color='blue', linewidth=1, label='')
    #plt.axhline(y=eigenvalues[5]+0.02, linestyle='--', color='red', linewidth=1, label='')

    # Adding text annotations at the top of the lines
    plt.text(4, eigenvalues[3]+0.03, '4th', horizontalalignment='right', color='blue')
    plt.text(6, eigenvalues[5]+0.03, '6th', horizontalalignment='right', color='red')
    
        # Highlight the 4th and 6th points with different colors
    plt.scatter(4, eigenvalues[3], s=18, color='blue')
    plt.scatter(6, eigenvalues[5], s=18, color='red')

    # Customizing x and y ticks
    plt.xticks(np.arange(1, 16), fontsize=12, rotation=45)
    plt.yticks(fontsize=12, rotation=45)

    # Labels
    plt.xlabel('Order')
    plt.ylabel('Eigenvalue')

    # Adding legend
    #plt.legend()

    plt.savefig(res_dir + 'eigengap_visual.png', bbox_inches='tight', dpi=600)


def draw_clustermap(clus_labels, affinity, res_dir, scaling = 'standard', shape = (5,5)):
    
    """
    Function to draw the squared affinity matrix, sorted by cluster labels, thus highlighting
    the checkerboard structure along the main diagonal. 

    Parameters
    ----------
    clus_labels : pd.Series
        The cluster membership tags, one for each patient.
    affinity : np.array
        The squared similarity matrix.
    res_dir : str
        The directory path where to save the generated plot.
    scaling : str
        Variable that indicates which normalization to apply to the similarity matrix.
        This is purely visual and has no ripercussion on the analysis. 'no' does not apply
        a transformation, while 'standard' and 'minmax' respectively apply a sklearn
        StandardScaler and MinMaxScaler.
    shape : tuple, optional
        x and y dimensions of the plot.
    """
        
    # Sort the affinity matrix by cluster label 
    printable_labels = pd.Series([cl + f' N={list(clus_labels).count(cl)}' for cl in list(clus_labels)])
    indx = np.argsort(list(printable_labels))
    printable_labels = printable_labels[indx]
    data = affinity[:, indx]
    data = data[indx,:]

    # Fill diagonal with zeroes for a better plot
    np.fill_diagonal(data, 0)
    
    lut = dict(zip(sorted(printable_labels.unique()), sns.color_palette("colorblind")))
    row_colors = printable_labels.map(lut).to_numpy()

    # choose if you want to scale in someway the affinities
    if scaling == 'standard':
        data = StandardScaler().fit_transform(data)
    elif scaling == 'minmax':
        data = MinMaxScaler().fit_transform(data)
    elif scaling == 'no':
        pass

    # When normalized it is necessary to mirror the matrix over the diagonal again
    i_lower = np.tril_indices(data.shape[0], -1)
    data[i_lower] = data.T[i_lower] 

    img = sns.clustermap(data,
                        figsize = shape,
                        col_cluster = False,
                        row_cluster = False,
                        xticklabels = False,
                        yticklabels = False,
                        row_colors = row_colors,
                        cmap="Spectral_r",
                        dendrogram_ratio = (.1, .1),
                        cbar_pos=(1, .3, .02, .4),
                        )  
    for label in sorted(printable_labels.unique()):
        img.ax_col_dendrogram.bar(0, 0, color=lut[label], label=label, linewidth=0)
        img.ax_col_dendrogram.legend(title='Clusters', bbox_to_anchor=(-.2, .5), ncol=1)
        
    plt.savefig(res_dir + 'clustermap.png', bbox_inches='tight', dpi=600)    
    plt.show()


def features_ami(inputs, clus_labels, K, n_clusters, howmany = 20, plot_flag = False):

    """
    Function that computes the adjusted mutual information scores between the clustering using all
    the features available and each separate feature. These values can be assumed to be a sort
    of representation of the feature importance for the overarching clustering. 

    Parameters
    ----------
    inputs : tuple
        The first element is the feature matrix and the second is the metric to use to construct the
        separate affinity matrices that will be used for clustering. The implementation works for more
        than one dataset at a time (in case of integration).
    clus_labels : pd.Series
        The cluster membership tags, one for each patient.
    K : int
        Number of neighbors used to construct the similarity matrices.
    n_clusters : int
        Number of clusters that spectral clustering will generate.
    howmany : int, optional
        Number of top relevant features to plot.
    """
        

    ami = [np.empty(shape=(d.shape[-1])) for d, m in inputs]
    

    for ndtype, (dtype, metric) in enumerate(inputs):

        for nfeature, feature in enumerate(np.asarray(dtype).T):

            aff = snf.make_affinity(np.vstack(feature), K=K, mu=.5, metric=metric)

            aff = np.nan_to_num(aff)

            aff_labels = SpectralClustering(n_clusters=n_clusters, n_init = 100,
                        affinity = 'precomputed', assign_labels='cluster_qr',
                        random_state = random_state).fit_predict(aff)

            ami[ndtype][nfeature] = adjusted_mutual_info_score(clus_labels, aff_labels)

    
    for i,imp in enumerate(ami):
        imp = pd.Series(imp)
        
        scores = imp.copy()

        imp.index = inputs[i][0].columns
        imp = imp.sort_values(ascending=False).iloc[:howmany].round(4)
        imp = imp[::-1]
        fig = px.bar(imp, text_auto='.3',orientation='h', width=600, height=500)

        fig.update_layout(
            xaxis_title="Adjusted Mutual Information (AMI) score",
            yaxis_title="Feature",
            xaxis_range=[imp[0]-imp[0]/10,imp[-1]],
            font=dict(
                family="Helvetica",
                size=12,
            ),
            showlegend=False
        )

        if plot_flag:
            fig.show()
    
    return scores


def convert_pvalue_to_asterisks(pvalue):
    if pvalue <= 0.0001:
        return "****"
    elif pvalue <= 0.001:
        return "***"
    elif pvalue <= 0.01:
        return "**"
    elif pvalue < 0.05:
        return "*"
    return ""


def significant_reg(clusters, parameters, pvalues, res_dir):

    pvalues = pd.DataFrame(pvalues).iloc[:,:]

    for col in pvalues.columns:
        # added multiple testing
        pvalues[col] = multipletests(pvals = pvalues[col], alpha = .05, method = 'bonferroni')[1]
        pvalues[col] = pvalues[col].apply(convert_pvalue_to_asterisks)
        
    
    param_coefs = pd.DataFrame(parameters).iloc[:,:]
    
    plt.subplots(figsize=(6,3))
    
    print(MinMaxScaler(feature_range=(-1, 1)).fit_transform(param_coefs).T)
    sns.heatmap(StandardScaler().fit_transform(param_coefs).T, linewidth=2,
                linecolor='black', 
                yticklabels=param_coefs.columns, 
                xticklabels=list(clusters),
                center = 0,
                cmap = "RdBu_r",
                annot=pvalues.to_numpy().T, fmt="s")
    
    plt.savefig(res_dir + 'log_reg.png', bbox_inches='tight', dpi=600)
    
    return clusters, param_coefs, pvalues


def draw_pca(data, meta_data, color_variable, symbol_variable, name):

    transf = StandardScaler().fit_transform(data.to_numpy())

    pca = PCA(n_components=2)
    pca.fit(transf)
    pca_result = pca.transform(transf)

    df = pd.DataFrame()
    df['pca-one'] = pca_result[:,0]
    df['pca-two'] = pca_result[:,1]
    df[color_variable] = meta_data[color_variable]
    df[symbol_variable] = meta_data[symbol_variable]
    
    fig = px.scatter(df, x='pca-one', y='pca-two', color=color_variable, symbol=symbol_variable,
                     labels={
                         "pca-one": "PC1 (var = %.2f)" %
                                 pca.explained_variance_ratio_[0],
                         "pca-two": "PC2 (var = %.2f)" %
                                 pca.explained_variance_ratio_[1],
                     },
                     )
    fig.update_layout(font=dict(
                    size=10),
                    legend_title="",
                    width=600,
                    height=400)
    fig.update_traces(marker=dict(size=6,))
    fig.update_yaxes(showticklabels=False)
    fig.update_xaxes(showticklabels=False)
    #fig.write_image(results_dir + "pca.png", height=600, width=900, scale = 4)
    #fig.write_html(results_dir + "pca_all_genes.html")
    fig.show()


def mann_sign(data_tuples, clus_labels, res_dir):  

    returnable = None
    
    for tup in data_tuples:
        
        mann_pvalues = pd.DataFrame()
        l2fcs = []

        if tup[1] == 'Cytof' or tup[1] == 'TCR':
            

            for cluster in sorted(clus_labels.unique()):

                y_mod = pd.Series([1 if clus == cluster else 0 for clus in clus_labels])

                mann = []
                feat = []
                l2fc = []


                for col in tup[0].columns:

                    pos_comp = tup[0][col].loc[y_mod == 1].reset_index(drop=True)
                    neg_comp = tup[0][col].loc[y_mod == 0].reset_index(drop=True)


                    try:
                        _, p = mannwhitneyu(pos_comp, neg_comp, method="auto")
                    except Exception as e:
                        p = 1

                    mann.append(p)
                    feat.append(col)
                    l2fc.append(np.log2(np.mean(pos_comp)/np.mean(neg_comp)))

                    df = pd.DataFrame(list(zip(mann, feat)),
                       columns =['mann_pvalues', 'features'])
                    df['cluster'] = cluster

                mann_pvalues = pd.concat([mann_pvalues,df])
                l2fcs.append(l2fc)

            print(mann_pvalues)
            values = mann_pvalues.copy()
            values = values.pivot(columns='cluster', values='mann_pvalues')

            values = values.T
            values.columns = tup[0].columns
            

            for col in values:
                values[col] = multipletests(pvals = values[col], alpha = .05, method = 'bonferroni')[1]
                values[col] = values[col].apply(convert_pvalue_to_asterisks)
      
            
                
            plt.figure(figsize = (6,3))
            sns.heatmap(np.array(l2fcs).T,
                            linewidth=2,
                            linecolor='black',
                            xticklabels=values.index, 
                            yticklabels=values.columns,
                            center = 0,
                            cmap = "RdBu_r",
                            annot=values.T, fmt="s"
                                    )
            
            plt.savefig(res_dir + tup[1] + '_clusters_testing.png', bbox_inches='tight', dpi=600)

            if returnable:
                returnable.append((values, pd.DataFrame(l2fcs,columns = values.columns), tup[1]))
            else:
                returnable = [(values, pd.DataFrame(l2fcs,columns = values.columns), tup[1])]
                
    if returnable:
        return returnable
    

def micro_analysis(data, clus_labels, micro_taxa, res_dir):

    pandas2ri.activate()
    importr('multidiffabundance')
    ro.r('''
        source('src/modules/tools.r')
    ''')

    da = ro.globalenv['diff_abundance']
            
    people = data['Donor'].copy()
    micro_data = data.drop(columns = ['Donor']).T.copy()
    micro_data.columns = people
    count_data = data.copy()
    
    
    meta_data = pd.DataFrame(people)
    meta_data['cluster'] = clus_labels.copy()
    meta_data['cluster'] = [int(cl.split(' ')[1]) for cl in clus_labels]
    meta_data.set_index('Donor',inplace = True)
    meta_data = pd.get_dummies(meta_data,columns=['cluster'], prefix='cluster', prefix_sep='_')

    form = ['~ cluster_'+ str(c.split(' ')[1]) for c in clus_labels.unique()]

    res = da(count_data.set_index('Donor'), meta_data, form)
    res = ro.conversion.rpy2py(res)

    sel = res[res['qvalue'] < 0.05]
    sel = sel.groupby(['taxa','variable']).size().reset_index(name='counts')
    sel = sel[sel['counts'] >= 2]
    taxa = sel['taxa'].unique()

    
    if len(taxa) > 0:
        
        sel['taxa'] = [list(micro_taxa[micro_taxa['taxon_id'] == taxon]['taxon_name'])[0] for taxon in sel['taxa']]
        res['taxa'] = [list(micro_taxa[micro_taxa['taxon_id'] == taxon]['taxon_name'])[0] for taxon in res['taxa']]
        sel.reset_index(drop = True, inplace = True)
        
        sel['effectsize'] = 0
        for index, row in sel.iterrows():
            effectsize = float(res.loc[(res['variable'] == row['variable']) & (res['taxa'] == row['taxa']) & (res['method'] == 'DESeq2') ]['effectsize'])
            sel.loc[index, 'effectsize'] = effectsize

        

        sig  = sel.copy()
        sig['family'] = [list(micro_taxa[micro_taxa['taxon_name'] == taxon]['family'])[0] for taxon in sig['taxa']]
        sig['variable']=['C ' + g.split('_')[1] for g in sig['variable']]
        

        order = sorted(clus_labels.unique()).copy()

        fig = px.strip(sig.loc[abs(sig['effectsize']) > 1], x='effectsize', y='family', color='variable', stripmode = "group",
                    color_discrete_sequence=px.colors.qualitative.Pastel,
                    labels={
                        "effectsize": "Log2FoldChange",
                        "family": "Family",
                    },
                    category_orders={'variable': order},
                    )
    

        fig2 = px.strip(sig.loc[abs(sig['effectsize']) <= 1], x='effectsize', y='family', color='variable', stripmode = "group", 
                                color_discrete_sequence= len(sig.loc[abs(sig['effectsize']) <= 1]['variable'].unique()) * ['gray'])

        fig = go.Figure(data = fig.data + fig2.data)
            
        fig.update_traces(marker=dict(line_width=1, 
                                    symbol='circle',
                                    size=10))
        
        for i in range(len(clus_labels.unique()), 2*len(clus_labels.unique())):
            
            try:
                fig.data[i]['showlegend'] = False
                fig.data[i]['marker_symbol'] = 'circle-x'
            except IndexError:
                pass

        fig.update_layout(font=dict(
                        size=12),
                        legend_title="Cluster",
                        xaxis=dict(
                            #showgrid=False,
                            showline=True,
                            linecolor='rgb(102, 102, 102)',
                            tickfont_color='rgb(102, 102, 102)',
                            showticklabels=True,
                            #dtick=10,
                            #ticks='outside',
                            tickcolor='rgb(102, 102, 102)',
                        ),
                        xaxis_title="Log2FC",
                        yaxis_title="Family",
                        boxmode='group',
                        boxgroupgap=1,
                        height=1000,
                        width=800,
                        paper_bgcolor='white',
                        plot_bgcolor='white',
                        )
        
        fig.add_vline(x=0, line_width=2)
        fig.add_vline(x=1, line_width=1, line_dash='dash')
        fig.add_vline(x=-1, line_width=1, line_dash='dash')
        
        sig.to_csv(res_dir + 'sig_microbiome.csv', index = False)
        fig.write_image(res_dir + "LFC_taxa_scatter.png", height=1000, width=800, scale = 4, engine="kaleido")
        fig.show()





        ############### select only top 10 by effect size
        sel = sel.sort_values(by = ['effectsize'], key=abs, ascending = False).groupby('variable').head(10)
        sel.reset_index(drop = True, inplace = True)


        returnable = []
        for _, row in sel.iterrows():
            effectsize = float(res.loc[(res['variable'] == row['variable']) & (res['taxa'] == row['taxa']) & (res['method'] == 'DESeq2') ]['effectsize'])
            if abs(effectsize) >= 1:
                returnable.append({'source': row['taxa'], 'target': 'C ' + row['variable'].split('_')[1], 
                                        'weight': (effectsize/2), 'origin': 'Micro'})
        return returnable
    else:
        print('No Differentially Abundant Bacteria, select different parameters.')


def draw_network(network, clus_labels, res_dir):
        
    network = pd.DataFrame(network)

    icons = {
        "Cluster": "icons/people.png",
        "LR_reg": "icons/demo.png",
        "Cytof_clusters": "icons/pbmc.png",
        "TCR": "icons/antibody.png",
        "Micro": "icons/bacteria.png",
        "GSEA": "icons/rna.png",
        "GSEA_bacteria": "icons/bacteria_pl3.png",
    }

    images = {k: os.path.abspath(fname) for k, fname in icons.items()}


    net = Network(height="1440px", width="100%",
                  font_color="black",
                  notebook=True,
                  cdn_resources='in_line',
                  )
    

    net.barnes_hut(central_gravity=1, spring_length=100, overlap=1, damping=1)

    lut = dict(zip(set(network['target']), sns.color_palette("colorblind")))


    for _,s in enumerate(network['source'].unique()):
        ori = network[network['source'] == s]['origin'].iloc[0]
        net.add_node(s, s, 
                    title=s, size = 30,
                    font_size=100,
                    shape='image',
                    image = images[ori]
                    )

    for _,d in enumerate(network['target'].unique()):
        s = pd.Series(clus_labels).value_counts()[str(d)]*3
        net.add_node(d, d, color = rgb2hex(int(lut[d][0]*255), int(lut[d][1]*255),int(lut[d][2]*255)),
                    size = int(s)*1.5 if int(s) < 100 else 100, 
                    title=d,
                    font_size=100,
                    shape = 'image',
                    image = images['Cluster']
                    )
        



    sources = network['source']
    targets = network['target']
    weights = network['weight']

    edge_data = zip(sources, targets, weights)

    for i,e in enumerate(edge_data):
                    src = e[0]
                    dst = e[1]
                    w = e[2]
                    edge_color = 'red' if w > 0 else 'blue'
                    edge_width = abs(w)
                    net.add_edge(src, dst, color = edge_color, value=edge_width, edge_width = edge_width)



    #net.toggle_physics(False)
    
    net.set_edge_smooth('cubicBezier')
    net.show_buttons(filter_= True
                     #["physics", "layout"]
                   )
    net.save_graph(res_dir + 'graph.html')


    #from IPython.display import display, HTML
    # Read the contents of the HTML file
    #with open(res_dir + 'graph.html', 'r') as file:
    #    html_content = file.read() 
    # Display the HTML content
    #display(HTML(html_content))
    
    return(net)


def ami_boxplots(amis, features, res_dir, top_x_features = 20):

    scores = pd.DataFrame([list(el) for el in amis], columns = features)
    scores = scores.reindex(scores.mean().sort_values(ascending=False).index, axis=1)

    scores = scores.iloc[:,:top_x_features]

    sns.set_style("whitegrid")
    sns.set_palette("deep")

    plt.figure(figsize=(8,8))

            
    sns.swarmplot(data=scores,linewidth=.5, orient="h", size=2)

    #sns.stripplot(data=scores, color="0.25", size=2, orient="h")

    sns.boxplot(data=scores, dodge= False, saturation = 0, width=0.3, showfliers=False, orient="h")

    plt.tight_layout()
    plt.legend([],[], frameon=False)
    plt.ylabel(f'Top {top_x_features} Features')
    #plt.grid(axis='x')
    plt.grid(axis='y')
    plt.xlabel("AMI scores (Feature importance)")
    plt.savefig(res_dir + "AMI_boxplots.png", dpi=600, bbox_inches='tight')
    plt.show()


def alpha_diversity_microbiome(abundances, patients, clus_labels, res_dir):

    shannon = alpha_diversity('shannon', abundances, patients)
    simpson = alpha_diversity('simpson', abundances, patients)
    chao1 = alpha_diversity('chao1', abundances, patients)

    stat = pd.DataFrame()
    stat['Shannon index'] = shannon
    stat['Simpson index'] = simpson
    stat['Chao1 richness'] = chao1
    stat.reset_index(drop = True, inplace=True)
    stat['Cluster'] = clus_labels

    print(stat)


    fig, axes = plt.subplots(1, 3, figsize=(14, 8))
    
    for i, metric in enumerate(['Shannon index', 'Simpson index', 'Chao1 richness']):

        order = sorted(stat['Cluster'].unique())
        pairs = list(combinations(order,2))
        dunn = sp.posthoc_dunn([stat[stat['Cluster'] == uni][metric] for uni in order], p_adjust = 'bonferroni')

        dunn = [str(convert_pvalue_to_asterisks(dunn.iloc[int(pair[0].split(' ')[1]), int(pair[1].split(' ')[1])])) for pair in pairs]
        dunn = ['ns' if el == '' else el for el in dunn]

        sns.swarmplot(data=stat, y = metric, x = 'Cluster',hue = 'Cluster', linewidth=.5, size=4, order=order, ax=axes[i])
        sns.boxplot(data=stat, y = metric, x = 'Cluster',hue = 'Cluster', dodge= False, saturation = 0, width=0.5, showfliers=False, order=order, ax=axes[i])
        plt.tight_layout()
        axes[i].text(0.02, 0.99, f"Kruskal-Wallis pvalue = {kruskal(*[stat[stat['Cluster'] == uni][metric] for uni in order]).pvalue:.3e}", ha="left", va="top", transform=axes[i].transAxes)

        annotator = Annotator(axes[i], pairs, data=stat,  y = metric, x = 'Cluster', order=order)
        annotator.configure(loc='outside')
        annotator.set_custom_annotations(dunn)
        #annotator.configure(test='Mann-Whitney', text_format='star', loc='outside', comparisons_correction="Bonferroni", verbose=2)
        annotator.annotate()
    
        axes[i].get_legend().remove()

    #plt.legend([],[], frameon=False)

    fig.savefig(res_dir + 'alpha_diversity_microbiome.png', dpi=600, bbox_inches='tight')
    plt.show()
