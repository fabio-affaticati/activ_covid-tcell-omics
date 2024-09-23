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

from sklearn.manifold import MDS
from sklearn.cluster import KMeans
from scipy.spatial.distance import cdist



random_state = 42



def eigengap_visual(affinity, to_annotate, res_dir, plot_name = 'eigengap_visual.png'):

    """
    Calculate the Laplacian matrix from the affinity matrix provided. The eigenvalues of the
    Laplacian are then calculated and plotted in ascending fashion.

    Parameters
    ----------
    affinity : np.array
        The similarity matrix for which to calculate the eigenvalues.
    to_annotate : list
        The eigenvalues to annotate with dashed lines and text.
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
    
    for i,annotation in enumerate(to_annotate):
    # Adding text annotations at the top of the lines
        if annotation == 1:
            text = '1st'
        elif annotation == 2:
            text = '2nd'
        elif annotation == 3:
            text = '3rd'
        else:
            text = f'{annotation}th'
        
        if (len(to_annotate) == 2) and (i == 0):
            color = 'blue'
        else:
            color = 'red'
        plt.text(annotation, eigenvalues[annotation-1]+0.03, text, horizontalalignment='right', color=color)
        plt.axhline(y=eigenvalues[annotation-1]+0.02, linestyle='--', color=color, linewidth=1, label='')
        plt.scatter(annotation, eigenvalues[annotation-1], s=18, color=color)

    # Customizing x and y ticks
    plt.xticks(np.arange(1, 16), fontsize=12, rotation=45)
    plt.yticks(fontsize=12, rotation=45)

    # Labels
    plt.xlabel('Order')
    plt.ylabel('Eigenvalue')

    plt.savefig(res_dir + plot_name, bbox_inches='tight', dpi=600)


def draw_clustermap(clus_labels, affinity, res_dir, scaling = 'standard', shape = (5,5), plot_name = 'clustermap.png'):
    
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
        
    plt.savefig(res_dir + plot_name, bbox_inches='tight', dpi=600)    
    plt.show()



def plot_coclustering_probability(stability_df, comparison_df, clus_labels, res_dir):

    comparison_df['Clusters'] = clus_labels
    comparison_df = comparison_df[['Clusters', 'Donor']]
    stability_df.rename(columns={'level_0': 'Subjectnr', 0:'Co-ClusteringProbability'}, inplace=True)
    stability_df = stability_df.query('Subjectnr != Donor')
    stability_df['SortedCombination'] = stability_df.apply(lambda row: ''.join(sorted([row['Subjectnr'], row['Donor']])), axis=1)

    stability_df = stability_df.drop_duplicates(subset='SortedCombination')
    stability_df.sort_values(by='SortedCombination', inplace=True)
    stability_df = pd.merge(stability_df, comparison_df, on = 'Donor')
    stability_df.sort_values(by='Clusters', inplace=True)
    ax = sns.displot(data=stability_df, x="Co-ClusteringProbability", kind="kde", hue = 'Clusters')
    plt.grid(axis = 'y')
    plt.savefig(res_dir + 'CoClusteringProbabilities.png', dpi=600)

    

def plot_mds_projection_cooccurrence(cooccurrence, clus_labels, clus_labels_RNA, res_dir, random_state = 42):

    # Convert similarity matrix to dissimilarity matrix (required for MDS)
    dissimilarity_matrix = 1 - cooccurrence

    # Apply Multidimensional Scaling (MDS)
    mds = MDS(n_components=2, dissimilarity='precomputed', random_state=42, max_iter = 1000)
    embedding = mds.fit_transform(dissimilarity_matrix)

    # Create a DataFrame for Seaborn
    df = pd.DataFrame({'X': embedding[:, 0], 'Y': embedding[:, 1], 'Cluster': clus_labels})
    df.sort_values(by='Cluster', inplace=True)

    # Plot the scatter plot with Seaborn
    plt.figure(figsize=(8, 6))
    sns.scatterplot(x='X', y='Y', hue='Cluster', data=df, palette='Set2', edgecolor='k', style = 'Cluster', s=100)


    # Show the plot
    plt.title('MDS stability - spectral clustering')
    plt.xlabel('Dimension 1')
    plt.ylabel('Dimension 2')
    plt.savefig(res_dir + 'MDS_stability.png', bbox_inches='tight', dpi=600)
    
    
    
    df = pd.DataFrame({'X': embedding[:, 0], 'Y': embedding[:, 1], 'Cluster': clus_labels_RNA})
    kmeans = KMeans(n_clusters=7, random_state=0, n_init="auto").fit(df[['X', 'Y']])
    df['kmeans_labels'] = kmeans.labels_


    # Compute Euclidean distances to assigned centroids
    distances = cdist(df[['X', 'Y']], kmeans.cluster_centers_, metric='euclidean')
    distances = pd.Series([d[l] for l, d in zip(kmeans.labels_, distances)])

    # Set a threshold for identifying distant points
    threshold = .2  # Adjust as needed

    # Identify distant points
    df.loc[df[distances > threshold].index, 'Cluster'] = '-1'
    df['Cluster'].replace({'-1': 'Unstable'}, inplace=True)
    df.sort_values(by='Cluster', inplace=True)

    plt.figure(figsize=(8, 6))
    sns.scatterplot(x='X', y='Y', hue='Cluster', data=df, palette='Set2', edgecolor='k', style = 'Cluster', s=100)


    # Show the plot
    plt.title('MDS stability - Kmeans')
    plt.xlabel('Dimension 1')
    plt.ylabel('Dimension 2')
    plt.savefig(res_dir + 'MDS_stability_kmeans.png', bbox_inches='tight', dpi=600)
    
    

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



def mann_sign(data_tuples, clus_labels, res_dir):  

    returnable = None
    
    for tup in data_tuples:
        
        mann_pvalues = pd.DataFrame()
        l2fcs = []

        if tup[1] == 'Cytof' or tup[1] == 'TCR':
            

            for cluster in sorted(set(clus_labels)):

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
    meta_data['cluster'] = clus_labels
    meta_data['cluster'] = [int(list(clus_labels.unique()).index(cl)) for cl in clus_labels]

    meta_data.set_index('Donor',inplace = True)
    meta_data = pd.get_dummies(meta_data,columns=['cluster'], prefix='cluster', prefix_sep='_')

    form = ['~ cluster_'+ str(list(clus_labels.unique()).index(c)) for c in clus_labels.unique()]

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
        print([list(clus_labels.unique())[int(g.split('_')[1])] for g in sig['variable']])
        sig['variable'] = [list(clus_labels.unique())[int(g.split('_')[1])] for g in sig['variable']]
        sig.rename(columns = {'counts': 'number_of_positive_tests '}, inplace = True)
        

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
                            showline=True,
                            linecolor='rgb(102, 102, 102)',
                            tickfont_color='rgb(102, 102, 102)',
                            showticklabels=True,
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
        fig.write_image(res_dir + "LFC_taxa_scatter.png", height=1000, width= 1800, scale = 2, engine="kaleido")
        #fig.show()


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


    fig, axes = plt.subplots(1, 3, figsize=(14, 10))
    
    for i, metric in enumerate(['Shannon index', 'Simpson index', 'Chao1 richness']):

        order = list(stat['Cluster'].unique())
        print(order)
        pairs = list(combinations(order,2))
        print(pairs)
        dunn = sp.posthoc_dunn([stat[stat['Cluster'] == uni][metric] for uni in order], p_adjust = 'bonferroni')
        print(dunn)
        dunn = [str(convert_pvalue_to_asterisks(dunn.iloc[int(order.index(pair[0])), int(order.index(pair[1]))])) for pair in pairs]
        dunn = ['ns' if el == '' else el for el in dunn]

        sns.swarmplot(data=stat, y = metric, x = 'Cluster',hue = 'Cluster', linewidth=.5, size=4, order=order, ax=axes[i])
        sns.boxplot(data=stat, y = metric, x = 'Cluster',hue = 'Cluster', dodge= False, saturation = 0, width=0.5, showfliers=False, order=order, ax=axes[i])
        plt.tight_layout()
        axes[i].text(0.02, 0.99, f"Kruskal-Wallis pvalue = {kruskal(*[stat[stat['Cluster'] == uni][metric] for uni in order]).pvalue:.3e}", ha="left", va="top", transform=axes[i].transAxes)

        annotator = Annotator(axes[i], pairs, data=stat,  y = metric, x = 'Cluster', order=order)
        annotator.configure(loc='outside')
        annotator.set_custom_annotations(dunn)
        annotator.annotate()
        # rotate axis labels by 45 degrees
        
        axes[i].get_legend().remove()
        axes[i].set_xticklabels(axes[i].get_xticklabels(), rotation=45)

    fig.savefig(res_dir + 'alpha_diversity_microbiome.png', dpi=600, bbox_inches='tight')
    plt.show()
