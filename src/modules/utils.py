import pandas as pd
import numpy as np
import re
import os

import matplotlib.pyplot as plt
import seaborn as sns

from pybiomart import Server

from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import KFold
from sklearn.cluster import SpectralClustering

import statsmodels.api as sm
from scipy.stats import mannwhitneyu
from scipy.sparse import csgraph
from statsmodels.sandbox.stats.multicomp import multipletests
from statsmodels.tools.tools import add_constant
import statistics

import snf

import warnings

from gseapy import dotplot, gsea

from src.modules.py_deseq import py_DESeq2
from src.modules.plotting import significant_reg, convert_pvalue_to_asterisks



random_state = 42


def eigentaxa_encoding(micro_fused_rel, eigentaxa):
    
    pca = PCA(n_components=1)
    scaler = StandardScaler()
    
    dim_red_bacteria = micro_fused_rel.drop(columns = ['Donor']).T.copy()
    dim_red_bacteria.columns = micro_fused_rel['Donor']
    dim_red_bacteria['taxaname'] = dim_red_bacteria.index
    dim_red_bacteria = dim_red_bacteria.replace({"taxaname": dict(zip(eigentaxa['microbe'], eigentaxa['eigentaxa']))})

    dim_red_bacteria.reset_index(drop=True,inplace=True)
    agg = dim_red_bacteria.groupby('taxaname', as_index = False)

    pose = pd.DataFrame()

    for group_name, df_group in agg:
        df_group = pca.fit_transform(scaler.fit_transform(df_group.drop(columns='taxaname').to_numpy().T)).T
        pose[group_name] = df_group.flatten().tolist()

    return pose


def testing_metrics(data, mets, dimensions, dataset_name):

    warnings.simplefilter(action='ignore', category=FutureWarning)

    def snf_mod(data, metric, K, mu = .5):

        transf = data.copy() 
        
        affinity = snf.make_affinity(transf, 
                                        metric=metric,
                                        K=K, 
                                        mu = mu,
                                        normalize = True
                                        )
        if len(metric) > 1:
            affinity = snf.snf(affinity, K=K)

        return affinity

    silhouettes = pd.DataFrame(columns=['metric', 'dimension', 'n_clusters','silhouette','dataset'])

    for n_clusters in [2,3,4,5,6,7,8]:

        for met in mets:
            #print(f'Metric {met}')

            for dimension in dimensions:

                affinity = snf_mod(data, [met], dimension)
                
                graph = affinity.copy()
                laplacian = csgraph.laplacian(graph, normed=True)
                
                eigenvalues = np.sort(np.linalg.eigvals(laplacian))
                
                clus_labels = utils.spectral_clustering_custom(affinity, n_clusters)
                silhouette = snf.metrics.silhouette_score(affinity, clus_labels)
                silhouettes = silhouettes.append({'metric': met,'dimension':dimension,'n_clusters':n_clusters,'silhouette':silhouette, 'dataset':dataset_name},ignore_index = True)

    g = sns.FacetGrid(silhouettes, col="n_clusters", hue="metric", sharey = True)
    g.map(sns.scatterplot, "dimension", "silhouette", alpha=.7)
    g.add_legend()    
    
    return silhouettes



def deseq2(count_matrix, design_matrix, gene_column, design_formula, dir_path):
    
    """
    Function to normalize the raw transcript counts for further analysis.
    The design formula for DESeq2 could potentially be a true comparison between phenotypes but
    in our case it will be only the 'intercept' as we do not need a DE analysis. 

    Parameters
    ----------
    micro_path : str
        Location path of the microbiome data view.
    percent : float
        Percentage of allowed missing values (zeros) for the bacteria that will be kept.
    flag : str, optional
        Parameter that allows to choose between relative and absolute abundances to import.
        Choose between 'rel_abundance' and 'abundance'.

    Returns
    ----------
    DataFrame of the normalized transcript counts using DESeq's median of ratios.

    """

    dds = py_DESeq2(count_matrix = count_matrix,
                    design_matrix = design_matrix,
                    design_formula = design_formula,
                    gene_column = gene_column)

    dds.run_deseq() 
    dds.get_deseq_result()
    dds_res = round(dds.normalized_count())
    dds_res.reset_index(drop = True, inplace = True)
    count_matrix.iloc[:,1:] = dds_res.iloc[:,:-1]
    count_matrix.to_csv(dir_path, sep='\t')
    return count_matrix



def gene_mapping(normalized_counts, dir_path):

    """
    Function that connects to the Ensembl website and replaces the gene ids with their human interpretable
    gene names. 

    Parameters
    ----------
    normalized_counts : pd.DataFrame
        DataFrame of normalized transcript counts. The first column has to be called 'genename'. It has to contain 
        the Ensembl ids of the genes. The subsequent columns are named after the samples.
    results_dir : str
        Path to directory where to save the processed data to avoid rerunning this function at every iteration.

    Returns
    ----------
    DataFrame of the normalized transcript counts in which the 'genename' column now contains the readable gene names.

    """    

    server = Server(host='http://www.ensembl.org')

    dataset = (server.marts['ENSEMBL_MART_ENSEMBL']
                      .datasets['hsapiens_gene_ensembl'])

    res = dataset.query(attributes=['ensembl_gene_id', 'external_gene_name'],)
    
    print('Finished querying!')
    
    normalized_counts.replace({"genename": dict(zip(res['Gene stable ID'], res['Gene name']))}, inplace=True)
    
    print('Finished mapping!')
    
    normalized_counts.dropna(how='any', inplace=True)
    normalized_counts.reset_index(drop = True, inplace = True)
    #normalized_counts.to_csv(dir_path + 'normalized_counts_renamed.csv', sep='\t')
    
    return normalized_counts


def aggregate_modules(typ, normalized_count, data_dir, dir_path):

    """
    Function that uses a user defined operation to aggregate the transcript counts into the 
    382 transcriptional modules of BloodGen3 from Altman et al, 2021  (https://doi.org/10.1038/s41467-021-24584-w). 

    Parameters
    ----------
    typ : str
        Identifier of the operation to use. 'pca' is the trasformation suggested by the authors. In this case
        the value of the first PC obtained from the genes comprising a module is assigned to the module itself.
    normalized_counts : pd.DataFrame
        DataFrame of normalized transcript counts. The first column has to be called 'genename'. It has to contain 
        the Ensembl ids of the genes. The subsequent columns are named after the samples.
    data_dir : str
        Path to the directory of the module_list file for the mapping of the genes into the respective modules.
    results_dir : str
        Path to directory where to save the processed data to avoid rerunning this function at every iteration.

    Returns
    ----------
    DataFrame of the aggregated gene counts.

    """    
    
    module_list = pd.read_csv(data_dir + 'module_list.csv')
    
    normalized_count['genename'] = [x.upper() for x in normalized_count['genename']]
    renamed = normalized_count.replace({"genename": dict(zip(module_list['Gene'], module_list['Module'] + '.' + module_list['Function']))})
    modules = [x for x in renamed['genename'].unique() if re.match('M[0-9]{1,3}\.[0-9]{1,3}', x)]
    
    # Keep only data present in the modules 
    renamed = renamed.loc[renamed['genename'].isin(modules)]
    
    if typ.startswith('pca') | typ.endswith('pca'):
        
        pca = PCA(n_components=1)
        scaler = StandardScaler()
        agg = renamed.groupby('genename', as_index = False)

        pose = pd.DataFrame()
        for group_name, df_group in agg:
            df_group = pca.fit_transform(scaler.fit_transform(df_group.drop(columns='genename').to_numpy().T)).T
            pose[group_name] = df_group.flatten().tolist()

    pose.to_csv(dir_path + 'pca_aggregated.csv', sep='\t')

    return pose


def spectral_clustering_custom(affinity, n_clusters):

    """
    Wrapper function for the spectral clustering algorithm. It uses the precomputed similarity matrix provided.
    
    Parameters
    ----------
    affinity : np.array
        The similarity matrix for which the clustering has to be performed.
    n_clusters : int
        Number of clusters that spectral clustering will generate.
    Returns
    ----------
    pd.Series of length equal to the number of the patients for which a similarity was calculated,
    each element is a cluster label in the form 'C x'.

    """   

    clus_labels = SpectralClustering(n_clusters=n_clusters, n_init = 100,
                    affinity = 'precomputed', assign_labels='cluster_qr',
                    random_state = random_state).fit_predict(affinity)
    
    clus_labels = ['C ' + str(el) for el in list(clus_labels)]

    print(f'Silhouette score: {snf.metrics.silhouette_score(affinity, clus_labels)}')

    return pd.Series(clus_labels)


def lr_custom(data, meta, clus_labels, res_dir): 
    
    X = data.copy()
    X['labels'] = clus_labels

    X_meta = pd.merge(X, meta, how = 'inner', left_on = 'Donor', right_on = 'Subjectnr')
    print(f'Size of metadata: {X_meta.shape}')
    
    X_meta.drop(columns = ['Donor', 'Subjectnr'], inplace=True)
    X_meta = X_meta.dropna()
    X_meta.reset_index(drop = True, inplace=True)
    print(f'Size of metadata after NaNs removal: {X_meta.shape}')

    y = X_meta['labels']
    X_meta.drop(columns = ['labels'], inplace=True)
    common_cols = [col for col in set(X_meta.columns).intersection(meta.columns)]
    X_meta = X_meta[common_cols]
    
    X_meta['smoking'].replace({1.0: 'Yes', 0.0: 'No'}, inplace=True)
    
    # Get dummy variables
    for categorical_col in ['Gender', 'smoking', 'Group']:
        if categorical_col == 'Group':
            X_meta = pd.get_dummies(X_meta,columns=['Group'])
            X_meta.drop(columns = ['Group_control'], inplace=True)
        else:
            X_meta = pd.get_dummies(X_meta,columns=[categorical_col] ,drop_first=True)
  

    X = X.iloc[X_meta.index]
    X.drop(columns = ['Donor', 'labels'], inplace=True)
    X.reset_index(drop = True,inplace=True)
    X_meta.reset_index(drop = True,inplace=True)
    

    print(f'Shape of the regression input: {X_meta.shape}')
    
    # Normalize inputs aside for the binary features
    transf = X_meta.copy()
    transf.loc[:, ~X_meta.columns.isin(['Gender_m', 'Group_household member', 'Group_patient', 'smoking_Yes'])] = \
        StandardScaler().fit_transform(X_meta.loc[:, ~X_meta.columns.isin(['Gender_m',
                                                    'Group_household member', 'Group_patient', 'smoking_Yes'])].to_numpy())
    
    transf.columns = X_meta.columns
    
    
    clusters = []
    parameters = []
    pvalues = []
    

    y, ordered = zip(*sorted(zip(list(y), transf.index)))
    
    y = pd.Series(y)
    ordered = list(ordered)

    transf = transf.reindex(ordered, axis=0)
    transf.reset_index(drop=True,inplace=True)
    
    for cluster in y.unique():
        
        print(f"{cluster}")
        y_mod = pd.Series([1 if clus == cluster else 0 for clus in y])
        lr = sm.Logit(y_mod, add_constant(transf, has_constant = 'add')).fit(method = 'bfgs', maxiter = 1000, full_output = True)
        print(lr.summary())
            
        clusters.append(cluster)
        parameters.append(lr.params)
        pvalues.append(lr.pvalues)

        
    return significant_reg(clusters, parameters, pvalues, res_dir)



def add_to_network(clusters, parameters, pvalues, origin, network):

    parameters = parameters[pvalues.replace('',np.nan).notnull()]
    parameters.index = clusters
    parameters.dropna(how = 'all')

    for index, row in parameters.iterrows():
        for i, cell in enumerate(row):
            if (not np.isnan(cell)) and (row.index[i] != 'const'):
                network.append({'source': row.index[i], 'target': index, 'weight': cell, 'origin': origin})

    return network     



# not used after development
def gsea_analysis_multi(counts, clus_labels, res_dir, path_sets):

        for cluster in sorted(clus_labels.unique()):
                
                results_pos = []
                results_neg = []

                #print(f'Cluster {cluster}')
                rna_data = counts.copy()

                ### done to make have 'genename' as first column, otherwise GSEA will not run properly
                rna_data.set_index('genename', inplace =True)
                rna_data.reset_index(inplace = True)
                
                clusters = pd.DataFrame()
                clusters['Cluster'] = clus_labels
                clusters.index = rna_data.drop(columns = 'genename').columns
                clusters['Cluster'] = [clus if clus == cluster else 'Cluster Else' for clus in clusters['Cluster']]
                clusters['Cluster'] = clusters['Cluster'].str.replace(' ', '_')

                with open('temp.cls', 'w') as cl:
                        line = f"{len(list(clusters['Cluster']))} 2 1\n# {cluster.replace(' ','_')} Cluster_Else\n" 
                        cl.write(line) 
                        cl.write(' '.join(list(clusters['Cluster'])) + '\n')


                for sets in path_sets:
                        gs_res = gsea(data=rna_data,
                                gene_sets=sets,
                                cls = 'temp.cls',
                                permutation_type='phenotype',permutation_num=100,
                                outdir=None,
                                method='signal_to_noise',
                                verbose = False,
                                threads = 4, seed = random_state)
                        
                        gsea_results_pos = gs_res.res2d[(gs_res.res2d['NES'] > 0) & (gs_res.res2d['FDR q-val'] < 0.25)].reset_index(drop=True).sort_values(by=['NES', 'FDR q-val'], ascending=[False,True]).iloc[:10,:]
                        
                        gsea_results_pos['Pathway'] = sets
                        gsea_results_pos['cluster'] = cluster
                        
                        gsea_results_pos['Term']=[g.split(' (GO')[0] for g in gsea_results_pos['Term']]
                        results_pos.append(gsea_results_pos)

                        gsea_results_neg = gs_res.res2d[(gs_res.res2d['NES'] < 0) & (gs_res.res2d['FDR q-val'] < 0.25)].reset_index(drop=True).sort_values(by=['NES', 'FDR q-val'], ascending=[True,True]).iloc[:10,:]
                        
                        gsea_results_neg['Term']=[g.split(' (GO')[0] for g in gsea_results_neg['Term']]
                        gsea_results_neg['Pathway'] = sets
                        gsea_results_neg['cluster'] = cluster
                        results_neg.append(gsea_results_neg)


                if os.path.isfile('temp.cls'):
                        os.remove('temp.cls')
          
                results_pos = pd.concat(results_pos)
                results_neg = pd.concat(results_neg)
                
                try:
                    dotplot(pd.concat([results_pos, results_neg], axis=0).reset_index(drop=True),
                        column="FDR q-val",
                        x = 'Pathway',
                        title=f'Cluster {cluster}',
                        cmap=plt.cm.viridis,
                        top_term=10,
                        size=3,
                        xticklabels_rot=45,
                        ofname = res_dir + f'Cluster {cluster} GSEA dotplot multisets',
                        show_ring = True,
                        figsize=(8,20), cutoff=0.25)
                
                except ValueError:
                        pass
        

def gsea_analysis(counts, clus_labels, res_dir, path_sets):

        results_pos = []
        results_neg = []

        for cluster in sorted(set(clus_labels)):

                #print(f'Cluster {cluster}')
                rna_data = counts.copy()

                ### done to make have 'genename' as first column, otherwise GSEA will not run properly
                rna_data.set_index('genename', inplace =True)
                rna_data.reset_index(inplace = True)
                
                clusters = pd.DataFrame()
                clusters['Cluster'] = clus_labels
                clusters.index = rna_data.drop(columns = 'genename').columns
                clusters['Cluster'] = [clus if clus == cluster else 'Cluster Else' for clus in clusters['Cluster']]
                clusters['Cluster'] = clusters['Cluster'].str.replace(' ', '_')

                with open('temp.cls', 'w') as cl:
                        line = f"{len(list(clusters['Cluster']))} 2 1\n# {cluster.replace(' ','_')} Cluster_Else\n" 
                        cl.write(line) 
                        cl.write(' '.join(list(clusters['Cluster'])) + '\n')


                for sets in path_sets:
                        gs_res = gsea(data=rna_data,
                                gene_sets=sets,
                                cls = 'temp.cls',
                                permutation_type='phenotype', 
                                permutation_num=100,
                                outdir=None,
                                method='signal_to_noise',
                                verbose = False,
                                threads =4 , seed = random_state)
                        
                        gsea_results_pos = gs_res.res2d[(gs_res.res2d['NES'] > 0) & (gs_res.res2d['FDR q-val'] < 0.25)].reset_index(drop=True).sort_values(by=['NES', 'FDR q-val'], ascending=[False,True]).iloc[:10,:]
                        
                        gsea_results_pos['Pathway'] = sets
                        gsea_results_pos['cluster'] = cluster

                        gsea_results_pos['Term']=[g.split(' (GO')[0] for g in gsea_results_pos['Term']]
                        results_pos.append(gsea_results_pos)

                        gsea_results_neg = gs_res.res2d[(gs_res.res2d['NES'] < 0) & (gs_res.res2d['FDR q-val'] < 0.25)].reset_index(drop=True).sort_values(by=['NES', 'FDR q-val'], ascending=[True,True]).iloc[:10,:]
                        
                        gsea_results_neg['Term']=[g.split(' (GO')[0] for g in gsea_results_neg['Term']]
                        gsea_results_neg['Pathway'] = sets
                        gsea_results_neg['cluster'] = cluster
                        results_neg.append(gsea_results_neg)
                        

                if os.path.isfile('temp.cls'):
                        os.remove('temp.cls')
          

                try:
                    dotplot(pd.concat([gsea_results_pos, gsea_results_neg], axis=0).reset_index(drop=True),
                        column="FDR q-val",
                        title=f'Cluster {cluster} {sets} only',
                        cmap=plt.cm.viridis,
                        top_term=40,
                        size=4,
                        ofname = res_dir + f'Cluster {cluster} GSEA dotplot',
                        show_ring = True,
                        figsize=(8,20), cutoff=0.25)
                
                except ValueError:
                        pass


        results_pos = pd.concat(results_pos)
        results_neg = pd.concat(results_neg)
        
        returnable = []
        for gs in results_pos.reset_index(drop=True).to_dict(orient='records'):
                returnable.append({'source': gs['Term'], 'target': gs['cluster'], 
                                        'weight': (gs['NES']/2), 'origin': 'GSEA',
                                        'genes': gs['Lead_genes']
                                        })
        for gs in results_neg.reset_index(drop=True).to_dict(orient='records'):
                returnable.append({'source': gs['Term'], 'target': gs['cluster'], 
                                        'weight':  (gs['NES']/2), 'origin': 'GSEA',
                                        'genes': gs['Lead_genes']
                                        })
        return returnable

        

def cytof_clusters_analysis(cytof, clus_labels, res_dir):
        
    mann_pvalues = pd.DataFrame()
    l2fcs = []
    #setup = cytof.iloc[:,1:].copy()
    print(f'Cytof number of columns {len(cytof.columns)}')
    setup = cytof.copy()
    try:
        setup.drop(columns = 'Donor', inplace = True)
    except:
        print('No need to drop Donor')

    for cluster in sorted(clus_labels.unique()):

        y_mod = pd.Series([1 if clus == cluster else 0 for clus in clus_labels])

        mann = []
        feat = []
        l2fc = []

        for col in setup.columns:

            pos_comp = setup[col].loc[y_mod == 1].reset_index(drop=True)
            neg_comp = setup[col].loc[y_mod == 0].reset_index(drop=True)

            try:
                _, p = mannwhitneyu(pos_comp, neg_comp, method="auto")
            except Exception:
                p = 1

            mann.append(p)
            feat.append(col)
            l2fc.append(np.log2(np.mean(pos_comp)/np.mean(neg_comp)))

            df = pd.DataFrame(list(zip(mann, feat)),
                columns =['mann_pvalues', 'features'])
            df['cluster'] = cluster

        mann_pvalues = pd.concat([mann_pvalues,df])
        l2fcs.append(l2fc)


    values = mann_pvalues.copy()
    values = values.pivot(columns='cluster', values='mann_pvalues')



    values = values.T
    values.columns = setup.columns
    l2fcs = pd.DataFrame(l2fcs,columns = values.columns)

    for col in values:
        values[col] = multipletests(pvals = values[col], alpha = .05, method = 'bonferroni')[1]
        for i, _ in enumerate(values[col]):
            values[col][i] = convert_pvalue_to_asterisks(values[col][i])

    plt.figure(figsize = (6,12))
    sns.heatmap(np.array(l2fcs.T), linewidth=2,
                linecolor='black',
                xticklabels=values.index, 
                yticklabels=values.columns,
                center = 0,
                cmap = "RdBu_r",
                cbar_kws={"shrink": 0.5},
                annot=values.T, fmt="s"
                        )

    plt.savefig(res_dir + 'Cytof_metaclusters_testing.png', bbox_inches='tight', dpi=600)

    return (values, l2fcs, 'Cytof_clusters')



def gsea_analysis_bacteria(counts, clus_labels, bubble_dir, res_dir, path_sets):
        

        results_pos = []
        results_neg = []
        
        
        for cluster in sorted(clus_labels.unique()):
                
                if os.path.isfile('temp.cls'):
                        os.remove('temp.cls')
                        
                print(f'Cluster {cluster}')
                micro_data = counts.copy()
                ### done to make have 'taxaname' as first column, otherwise GSEA will not run properly
                micro_data.set_index('taxaname', inplace =True)
                micro_data.reset_index(inplace = True)
                
                clusters = pd.DataFrame()
                clusters['Cluster'] = clus_labels
                clusters.index = micro_data.drop(columns = 'taxaname').columns
                clusters['Cluster'] = [clus if clus == cluster else 'Cluster_Else' for clus in clusters['Cluster']]
                clusters['Cluster'] = clusters['Cluster'].str.replace(' ', '_')
                
                with open('temp.cls', 'w') as cl:
                        line = f"{len(list(clusters['Cluster']))} 2 1\n# {cluster.replace(' ','_')} Cluster_Else\n" 
                        cl.write(line) 
                        cl.write(' '.join(list(clusters['Cluster'])) + '\n')
                        
                for sets in path_sets:
                        gs_res = gsea(data=micro_data,
                                gene_sets=sets,
                                cls = 'temp.cls',
                                permutation_type='phenotype', 
                                permutation_num=1000,
                                outdir=None,
                                method='signal_to_noise',
                                verbose = False,
                                threads =4 , seed = random_state)
                        
                        gsea_results_pos = gs_res.res2d[(gs_res.res2d['NES'] > 0) & (gs_res.res2d['FDR q-val'] < 0.25)].reset_index(drop=True).sort_values(by=['NES', 'FDR q-val'], ascending=[False,True]).iloc[:10,:]
                        
                        gsea_results_pos['cluster'] = cluster

                        results_pos.append(gsea_results_pos)

                        gsea_results_neg = gs_res.res2d[(gs_res.res2d['NES'] < 0) & (gs_res.res2d['FDR q-val'] < 0.25)].reset_index(drop=True).sort_values(by=['NES', 'FDR q-val'], ascending=[True,True]).iloc[:10,:]
                        
                        gsea_results_neg['cluster'] = cluster
                        
                        results_neg.append(gsea_results_neg)
                        pd.DataFrame(gs_res.res2d).to_csv(f'{bubble_dir}bubbleplot_{cluster}.csv', sep='\t')

                if os.path.isfile('temp.cls'):
                        os.remove('temp.cls')
          

                try:
                    dotplot(pd.concat([gsea_results_pos, gsea_results_neg], axis=0).reset_index(drop=True),
                        column="FDR q-val",
                        title=f'Cluster {cluster}',
                        cmap=plt.cm.viridis,
                        top_term=40,
                        size=4,
                        ofname = res_dir + f'Cluster {cluster} GSEA bacteria dotplot',
                        show_ring = True,
                        figsize=(8,8), cutoff=0.25)
                
                except ValueError:
                        pass


        results_pos = pd.concat(results_pos)
        results_neg = pd.concat(results_neg)
        
        returnable = []
        for gs in results_pos.reset_index(drop=True).to_dict(orient='records'):
                returnable.append({'source': gs['Term'], 'target': gs['cluster'], 
                                        'weight': (gs['NES']/2), 'origin': 'GSEA_bacteria',
                                        'genes': gs['Lead_genes']
                                        })
        for gs in results_neg.reset_index(drop=True).to_dict(orient='records'):
                returnable.append({'source': gs['Term'], 'target': gs['cluster'], 
                                        'weight':  (gs['NES']/2), 'origin': 'GSEA_bacteria',
                                        'genes': gs['Lead_genes']
                                        })
        return returnable



def plot_heatmaps_gsea_bacteria(res, micro_data, clus_labels, micro_taxa, res_dir):


        rel_counts = micro_data.copy()

        taxa_list = [list(micro_taxa[micro_taxa['taxon_id'] == taxon]['taxon_name'])[0] for taxon in rel_counts.drop(columns = ['Donor']).columns]
        rel_counts.drop(columns = ['Donor'], inplace=True)
        rel_counts.columns = taxa_list
        
        for family in res:
                print(family)
                rel_counts['labels'] = clus_labels.copy()
                rel_counts['labels'] = np.where(rel_counts['labels'] != family['target'], 'C Else', family['target'])

                rel_counts.sort_values(by='labels', inplace=True)

                lut = dict(zip(rel_counts['labels'].unique(), sns.color_palette("colorblind")))
                row_colors = rel_counts['labels'].map(lut).to_numpy()

                img = sns.clustermap(rel_counts.iloc[:,rel_counts.columns.isin(family['genes'].split(';'))].T.values,
                                        standard_scale=0,
                                        figsize= (8,len(rel_counts.iloc[:,rel_counts.columns.isin(family['genes'].split(';'))].columns)/4),
                                        col_cluster = False,
                                        row_cluster = False,
                                        col_colors = row_colors,
                                        yticklabels=rel_counts.iloc[:,rel_counts.columns.isin(family['genes'].split(';'))].columns,
                                        xticklabels=False,
                                        cmap="Spectral_r",
                                        #metric='correlation',
                                        dendrogram_ratio=(.05, .05),
                                        cbar_pos=(1, .3, .02, .4),
                                        )
                img.figure.suptitle(family['source'])

                for label in sorted(rel_counts['labels'].unique()):
                        img.ax_col_dendrogram.bar(0, 0, color=lut[label],
                                                label=label, linewidth=0)
                        img.ax_col_dendrogram.legend(title='Clusters', bbox_to_anchor=(-.1, .3), ncol=1)

                plt.savefig(res_dir + 'heatmap_'+str(family['source'])+'_'+str(family['target'])+'.png', bbox_inches='tight', dpi=600)    



def stability_analysis(data, individuals, metric, dim , n_clusters, n_iterations = 10, n_ks = 10):

    silhouettes = []
    cooccurrence = pd.DataFrame(0, index=np.arange(len(individuals)), columns=list(individuals))
    cooccurrence['Donor'] = individuals

    for j in range(n_iterations):

        print(f"Run {j}:")
        
        kf = KFold(n_splits=n_ks, random_state=random_state+j, shuffle=True)

        for _, (train_index, _) in enumerate(kf.split(data)):
            
            
            run = data.copy()
            run['Donor'] = individuals
            run = run.iloc[train_index,:]
            run.reset_index(drop=True,inplace = True)

            affinity = snf.make_affinity(run.drop(columns='Donor'), metric = metric, K = dim, mu = .5, normalize = True)
            clus_labels = spectral_clustering_custom(affinity, n_clusters)

            silhouettes.append(snf.metrics.silhouette_score(affinity, clus_labels))

            one_hot = pd.DataFrame()
            one_hot['labels'] = clus_labels
            one_hot = pd.get_dummies(one_hot['labels'])
            one_hot['Donor'] = run['Donor']
            one_hot.set_index(keys = ['Donor'], inplace=True)
            one_hot = one_hot.T.astype(int)
            one_hot = one_hot.T.dot(one_hot)


            ### Could concat them all at the end to speed up execution
            cooccurrence = pd.concat([cooccurrence, one_hot.reset_index()]).groupby('Donor',as_index=False).sum()
            #print(one_hot.reset_index())


    cooccurrence.set_index(keys = ['Donor'], inplace = True)
    cooccurrence = cooccurrence/(n_iterations*(n_ks-1))
    np.fill_diagonal(cooccurrence.values, 0)

    print(f'\n\nMean silhouette over the runs: {statistics.mean(silhouettes)}, \
          standard deviation: {statistics.stdev(silhouettes)}')

    return cooccurrence