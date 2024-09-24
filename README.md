# Bridging immunotypes and enterotypes using a systems immunology approach.

All data and scripts to reproduce the results from the 'Bridging immunotypes and enterotypes using a systems immunology approach' study.

## 1. Summary:

Recent advances in multi-omics analyses have revolved mainly around well established, directly
determinable molecular relationships such as the protein synthesis cascade and epigenetic
mechanisms. However, determining the broad effects of disease and health maintenance must
occur in a truly holistic way. In this study, including 394 individuals, we found direct linkage of
systemic branches spanning human biological functions often not studied in conjunction, using
clinical data, gut microbial abundances, blood immune cell repertoires, blood bulk transcriptomic
and blood bulk T cell receptor data. Aggregation of patient similarities via similarity network
fusion was used to promote a data agnostic averaging of the distance across modalities. The gut
microbiome and blood immune cell repertoire were novelly found to be orthogonal datasets, only
bridged via the blood transcriptome. Even though prior CMV and SARS-CoV-2 infection co-shaped
the immune adaptive compartment, inflammation was highlighted as the most diversifying
marker.




## 2. Dependencies


### 2.1 Python dependencies
**Python Version:** 3.11.9

Data manipulation:
- numpy==1.24.2
- pandas==1.5.3
- openpyxl==3.1.5

Statistics:
- scipy==1.14.1
- statsmodels==0.14.3
- gseapy==1.1.3
- scikit-bio==0.6.2
- scikit-learn==1.5.2
- scikit-posthocs==0.9.0
- statannotations==0.6.0

Data visualisation:
- matplotlib==3.7.1
- seaborn==0.11.2
- plotly==5.24.1
- colormap==1.1.0
- networkx==3.3
- pyvis==0.3.2

Misc:
- ipython==8.27.0
- ipykernel==6.29.5
- rpy2==3.5.16
- pybiomart==0.2.0
- requests==2.32.3
- snfpy==0.2.2


### 2.2 R dependencies
**R Version:** 4.3.3
Data manipulation:
- tidyverse==2.0.0
- tidyr==1.3.1
- dplyr==1.1.4
- reshape2==1.4.4
- Matrix==1.6-5

Modeling:
- ROCR==1.0-11
- vegan==2.6-8
- kernlab==0.9-33
- multidiffabundance==0.0.1

Data visualisation:
- circlize==0.4.16
- ggplot2==3.5.1
- ggvenn==0.1.10
- RColorBrewer==1.1-3
- gridExtra==2.3
