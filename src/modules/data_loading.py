import numpy as np
import pandas as pd
import glob
import os


def metadata_import(meta_dir):

    """
    Import and process metadata view.

    Parameters
    ----------
    meta_dir : str
        Location path of the metadata dataset.

    Returns
    ----------
    DataFrame of processed metadata information.

    """
    
    replacements_subjectnr = {'Covid_':'Covid ', 'ERC_':'ERC ', 'HH_HA':'HHHA','HH HA':'HHHA', '_':'.'}
    replacements_CMV = {'<':'', '>':'', ',':'.', '-' : np.nan}
    replacements_groups = {'Covid HCW group 1':'household member', 'Covid HA group 5':'household member',
                                'Covid HCW group 3':'household member', 'Covid HCW group 5':'household member',
                               'Covid HCW group 4':'household member', 'Covid HCW group 2':'household member',
                                'Covid HA group 2':'household member', 'Covid HA group 1':'household member',
                               'Covid HA group 4':'household member', 'Covid HA group 3':'household member',
                                'mild patient' : 'patient', 'moderate patient' : 'patient',
                                'severe patient' : 'patient', 'possible patient' : 'patient'}
    
    # Load raw data
    meta = pd.read_excel(meta_dir)

    meta['Subjectnr'].replace(replacements_subjectnr, regex=True, inplace=True)

    # Replace marker with 1 - Long covid 
    meta['Long covid'].replace({'þ':1}, inplace=True)

    # Missing values for depression and anxiety scores
    meta[['BDI T-score (depression)', 'BAI T-score (anxiety)']] = meta[['BDI T-score (depression)', 'BAI T-score (anxiety)']].replace({'-': np.nan})

    meta['CMV IgG (U/mL)'] = pd.to_numeric(meta['CMV IgG (U/mL)'].replace(replacements_CMV, regex=True))
    meta['25-OH-Vit D'] = pd.to_numeric(meta['25-OH-Vit D'].replace({'-': np.nan}))
    meta['Covid 19 (1/0)'] = pd.to_numeric(meta['Covid 19 (1/0)'].replace({'-': np.nan}))

    # Drop missing BMI
    meta = meta[meta['BMI'] != 'not applicable']
    meta['BMI'] = pd.to_numeric(meta['BMI'])

    # Order by patient ID (Subjectnr)
    meta['Group'].replace(replacements_groups,inplace=True)
    meta.sort_values(by = ['Subjectnr'], inplace = True)
    meta.reset_index(drop=True, inplace=True)

    # Ignore columns
    return meta
    


def tcr_stats_import(tcr_stats_dir):

    """
    Import and process TCR summary statistics data used for TCR-seq dataset testing.

    Parameters
    ----------
    tcr_stats_dir : str
        Location path of the TCR statistics data.

    Returns
    ----------
    DataFrame of processed TCR stats information.

    """

    replacements_description = {'Covid_':'Covid ', 'ERC_':'ERC ', 'HH.HA':'HHHA', '_':'.', 'Covid HH0801':'Covid HHHA031'}
    
    TCR_stats = pd.read_csv(tcr_stats_dir + 'TCR_stats.txt', sep = ',', index_col=0)

    TCR_stats['totalcountCD4'] = np.where(TCR_stats['type']=='CD4', TCR_stats['totalcount'], 0)
    TCR_stats['entropyCD4'] = np.where(TCR_stats['type']=='CD4', TCR_stats['entropy'], 0)
    TCR_stats['tcrcovidmatchCD4'] = np.where(TCR_stats['type']=='CD4', TCR_stats['tcrcovidmatch'], 0)
    TCR_stats['tcrcovidfreqCD4'] = np.where(TCR_stats['type']=='CD4', TCR_stats['tcrcovidfreq'], 0)

    TCR_stats['totalcountCD8'] = np.where(TCR_stats['type']=='CD8', TCR_stats['totalcount'], 0)
    TCR_stats['entropyCD8'] = np.where(TCR_stats['type']=='CD8', TCR_stats['entropy'], 0)
    TCR_stats['tcrcovidmatchCD8'] = np.where(TCR_stats['type']=='CD8', TCR_stats['tcrcovidmatch'], 0)
    TCR_stats['tcrcovidfreqCD8'] = np.where(TCR_stats['type']=='CD8', TCR_stats['tcrcovidfreq'], 0)
    
    TCR_stats = TCR_stats.groupby(['description'], axis=0).sum()[['totalcountCD4','totalcountCD8',
                                                                  'entropyCD4','entropyCD8',
                                                                 'tcrcovidmatchCD4', 'tcrcovidmatchCD8',
                                                                 'tcrcovidfreqCD4', 'tcrcovidfreqCD8']]
    # make description a column and not an index
    TCR_stats.reset_index(inplace = True)

    for orig, rep in replacements_description.items():
        TCR_stats['description'] = TCR_stats['description'].str.replace(orig, rep)


    # Drop problematic individual
    TCR_stats = TCR_stats[TCR_stats['description'] != 'Covid HH053']

    TCR_stats.rename(columns = {'description': 'Donor'}, inplace=True)
    TCR_stats.sort_values(by = ['Donor'], inplace = True)
    TCR_stats.reset_index(drop=True,inplace=True)
    
    return TCR_stats


def import_all_cytof(cytof_dir):

    """
    Function that fuses the corresponding CD4 and CD8 cluster information for the patients, with their summary
    statistics.

    Parameters
    ----------
    cytof_dir : str
        Location path of the Flow Cytometry data view.
  

    Returns
    ----------
    Merged DataFrame of the cell counts and the relative metadata that will be used for testing.

    """
    
    replacements = {'Covid_':'Covid ', 'ERC_':'ERC ', 'HH_HA':'HHHA','HH HA':'HHHA', '_':'.'}
    
    cytof_meta = pd.read_csv(cytof_dir + 'CellType_mapping/celltype_mapping.csv', sep = ',')
    cytof_meta['Group'].replace({'control':'CO', 'preCovid':'PreCovid', 'household':'HH',
                                'patient mild':'PT_mild', 'patient moderate':'PT_moderate',
                                'patient severe':'PT_severe'},inplace=True)
    
    
    cytof_meta.rename(columns={'Unnamed: 0':'Donor'},inplace=True)
    cytof_meta['Donor'] = [('ERC ' + donor) if donor.startswith('CO') else ('Covid '+ donor) for donor in cytof_meta['Donor']]
    cytof_meta['Donor'].replace(replacements, regex=True, inplace=True)


    #cytof = pd.read_csv(cytof_dir + 'FlowSOM/freq_percent_allcells.csv', sep = ',', index_col = 0)
    cytof = pd.read_csv(cytof_dir + 'FlowSOM/freq_percent_allcells_allmarkers_agg10m_v2.csv', sep = ',', index_col = 0)
    
    cytof = cytof[cytof.Name != 'FlowSOM_c18_COVID-19 STUDY CyTOF FIXED 221121 HKM_CO161.fcs']
    cytof = cytof[cytof.Name != 'FlowSOM_c09_COVID-19 STUDY CyTOF FIXED 281021 HKM_CO203.fcs']
    cytof = cytof[cytof.Name != 'FlowSOM_c20_COVID-19 STUDY CyTOF FIXED 2 41121 HKM_H40d365.fcs']

    cytof.reset_index(drop = True, inplace=True)
    cytof['Donor'] = cytof['Name'].str.split('HKM').str[1]
    cytof['Donor'] = cytof['Donor'].str.split('_').str[-1]
    cytof['Donor'] = cytof['Donor'].str.split('.fcs').str[0]
    cytof.drop(columns = 'Name', inplace = True)

    cytof['Donor'] = [('ERC ' + donor) if donor.startswith('CO') else ('Covid '+ donor) for donor in cytof['Donor']]
    cytof['Donor'].replace(replacements, regex=True, inplace=True)

    # Translate metacluster tags into cell subpopulations
    #cytof.rename(columns={
    #                'Metacluster 1': 'CD3+CD4-CD8+CD45RA-CCR7+',
    #                'Metacluster 2': 'CD3+CD4-CD8+CD45RA-CCR7-',
    #                'Metacluster 3': 'CD3-CD19+CCR6+',
    #                'Metacluster 4': 'CD3-CD56+',
    #                'Metacluster 5': 'CD3-CD11c+CD14+',
    #                'Metacluster 6': 'CD3+CD4+CD45RA-CCR7+',
    #                'Metacluster 7': 'CD3+CD4+CD45RA-CCR7-'},
    #                  inplace = True)
    cytof.rename(columns={
                   'Metacluster 1': 'Sen_CD8+Temra',
                   'Metacluster 2': 'Sen_CD4+Tem',
                   'Metacluster 3': 'Exh_CD4+Tef',
                   'Metacluster 4': 'CD8+Tcm',
                   'Metacluster 5': 'Exh_CD4+Tcm',
                   'Metacluster 6': 'CD4+Tregs',
                   'Metacluster 7': 'NKT',
                   'Metacluster 8': 'Exh_CD4+Th2cm',
                   'Metacluster 9': 'NK',
                   'Metacluster 10': 'CD8+Tem',
                   'Metacluster 11': 'CD11cint',
                   'Metacluster 12': 'Exh_CD8+Eef|em',
                   'Metacluster 13': 'Act_NK',
                   'Metacluster 14': 'CD4+Tna',
                   'Metacluster 15': 'Fas+Tim3+CD11c+CD14+',
                   'Metacluster 16': 'Act_CCR4+CD11cintCD14int',
                   'Metacluster 17': 'CD8+Tna',
                   'Metacluster 18': 'DNT',
                   'Metacluster 19': 'CCR7+CD45RA+CCR6+CXCR5+LAG3+CD25+B',
                   'Metacluster 20': 'Act_CD4+Tef',
                   'Metacluster 21': 'CCR7-CD45RA+CCR6-CXCR5-Fas+B',
                   'Metacluster 22': 'CCR7+CD45RA+CCR6+CXCR5+LAG3-CD25-B',
                   },
                     inplace = True)
    

    cytof_meta.drop(columns={'Group'}, inplace=True)
    return pd.merge(cytof, cytof_meta, on ='Donor', how='inner').dropna(how='any')




def microbiome_import(micro_path, percent, flag = 'relative'):

    """
    Function to import and preprocess the gut microbiome abundances of the patients.

    In a secondary DataFrame, information for the bacteria is also reported.
    This include 'kingdom', 'phylum', 'class', 'order', 'family', 'genus, 'species' of the
    bacteria, as well as an identification id, 'taxon_id'.

    A choice beween relative and absolute abundances is given, as well as the percentage of individuals
    allowed to not be presenting the bacteria before it is discarded.

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
    A DataFrame of the processed abundances and a DataFrame for the classification of the bacteria.

    """
    
    
    micro_samples = pd.read_csv(micro_path + 'samples.csv', sep = ',', index_col=0)
    micro_samples.columns = [col.split('samples.')[1] for col in micro_samples.columns]

    micro_taxa = pd.read_csv(micro_path + 'taxa.csv', sep = ',', index_col=0)
    micro_taxa.columns = [col.split('taxa.')[1] for col in micro_taxa.columns]
    micro_taxa.drop(columns={'arrange_by_me'},inplace=True)
    
    micro_abundances = pd.read_csv(micro_path + 'abundances.csv', sep = ',', index_col=0)
    micro_abundances.columns = [col.split('abundances.')[1] for col in micro_abundances.columns]
  
    micro_taxa['genus'] = micro_taxa['genus'].fillna(micro_taxa['taxon_name'].str.split(' ').str[0])



    abundances = micro_abundances.pivot(index = 'sample_id', columns = 'taxon_id', values=flag)


    abundances.columns = [col for col in abundances.columns]
    abundances.reset_index(inplace=True)

    micro_fused = pd.merge(micro_samples, abundances, on = 'sample_id', how='inner')

    micro_fused['description'] = [row.split('-F-T')[0] for row in micro_fused['description']]
    micro_fused['description'] = micro_fused['description'].str.replace('COVHH-','Covid HH')
    micro_fused['description'] = micro_fused['description'].str.replace('ERCCO-','ERC CO')
    micro_fused['description'] = micro_fused['description'].str.replace('COVPT-','Covid PT')


    micro_fused.drop_duplicates(subset=['description'], keep='first', inplace=True)
    micro_fused.drop(columns=['sample', 'participant','location','timepoint','plate','lib_size', 'dna_pcr_reads', 'pool', 'dna_pcr_qubit', 'run','sample_id','control','project','condition', 'passes_qc'],inplace=True)

    micro_fused = micro_fused.loc[:, micro_fused.isnull().sum(axis=0) <= percent*micro_fused.shape[0]]
    micro_fused.fillna(0, inplace=True)

    micro_taxa[micro_taxa['taxon_id'].isin([col.split('_')[0] for col in micro_fused.columns[1:]])]
   
    micro_fused.rename(columns={'description': 'Donor'}, inplace = True)
    micro_fused.sort_values(by = ['Donor'], inplace = True)
    micro_fused.reset_index(drop=True,inplace = True)

    return micro_fused, micro_taxa


def rna_import(rna_path):

    """
    Function that merges the separate RNA-seq readcounts textual files .

    Parameters
    ----------
    rna_path : str
        Location path of the RNA-seq raw data view.
  
    Returns
    ----------
    Merged DataFrame of the RNA readcounts.

    """

    all_files = glob.glob(os.path.join(rna_path , "readcounts_*.txt"))
    start = True
    print(all_files)
    print(f'\nNumber of files: {len(all_files)}')
    for f in all_files:
        if start:
            start = False
            raw_rna = pd.read_csv(f, sep = '\t')
            print('Started!')
        else:
            raw_rna = pd.merge(raw_rna, pd.read_csv(f, sep = '\t'), how="inner", on='genename')
    print('Ended!')
    return raw_rna


def cleanup_rna(rna):

    """
    Final clean up of the RNA-seq dataset before gene count normalisation.
    Primarily the identifiers of the patients are here homogenized.

    Parameters
    ----------
    rna : pd.DataFrame
        DataFrame of the raw transcript counts that require a normalisation of their notation.
  
    Returns
    ----------
    DataFrame of the RNA readcounts with the cleaned id tags of the patients.
    A secondary auxiliary DataFrame is also returned that will be used to normalise the
    transcript counts with DESeq2 normalisation implementation.
    """
    
    # prepare name of columns
    rna.drop(rna.columns[rna.columns.str.startswith('Unnamed')], axis = 1, inplace=True)
    rna.columns = [col.split('/tmp_files/')[-1] for col in rna.columns]
    rna.columns = rna.columns.str.rstrip('.sam')


    rna.columns  = rna.columns.str.replace('_', ' ')
    rna.columns  = rna.columns.str.replace('-', '.')
    rna.columns = rna.columns.str.replace('Covid HH021  ', 'Covid HH021 ')
    rna.columns = rna.columns.str.replace('Covid HHHA021  ', 'Covid HHHA021 ')
    rna = rna[rna.genename != '_ambiguous']
    rna = rna[rna.genename != '_no_feature']
    rna = rna[rna.genename != '_unmapped']
    unique = [col.split(' L00')[0].rsplit(' S')[0] for col in rna.columns]
    unique = list(dict.fromkeys(unique))
    
    
    grouper = [next(p for p in unique if p in c) for c in rna.columns]
    results = rna.groupby(grouper, axis=1).sum()
    results['genename'] = [row.split('.')[0] for row in results['genename']]
    print(f"Number of unique transcripts: {len(results['genename'].unique())}")
    results = results.groupby(['genename'], axis=0).sum().reset_index()
    print(f"Number of patients: {len(results.columns)-1}")
    

    # placeholder column that will be used when normalising gene count with DESeq2.
    # DESeq2 asks in input the phenotypes for the samples to perform DE analysis but we are
    # only interested in the normalisation for the moment.
    design = ['place'] * (len(results.columns)-1)
    # save the list of the patients/samples for DESeq2
    colnames = results.columns[1:]

    return results, pd.DataFrame({'place' : design, 'colnames' : colnames })
    