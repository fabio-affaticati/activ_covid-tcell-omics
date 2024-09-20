library(tidyverse)
library(RColorBrewer)
devtools::install_github("microbiome/mia")

workingDir = "/Users/fabioaffaticati/Desktop/Work/activ_covid-and-omics"
setwd(workingDir) 

### pass relative abundances
Tab = read.csv('results/relative_abundances_forplot.csv', sep = '\t', row.names = 1)

### read or pass the taxa information
Tax = read.csv('data/microbiome_analysis/taxa.csv', sep = ',', row.names = 1)

### pass the cluster information
Meta = read.csv('data/clus_labels.csv', sep = '\t', row.names = 1)


Tab$taxa.taxon_id <-rownames(Tab)


dat <- Tab %>%
  pivot_longer(-taxa.taxon_id, names_to = "sample_id", values_to = 'rel_ab')

dat <- dat %>%
  left_join(Tax, by='taxa.taxon_id')

dat$sample_id <-gsub(".", " ", dat$sample_id, fixed=TRUE)

dat <- dat %>%
  left_join(Meta, by='sample_id')



#dat <- dat[dat$taxa.family %in% c('Bacteroidaceae', 'Lachnospiraceae', 'Prevotellaceae', 'Ruminococcaceae', 'Christensenellaceae'), ]

family_order <- c('Other',
                  'Bacteroidaceae', 
                  'Lachnospiraceae', 
                  'Prevotellaceae', 
                  'Ruminococcaceae',
                  #'Rikenellaceae',
                  'Christensenellaceae')

shouldBecomeOther <- !(dat$taxa.family %in% family_order)
dat$taxa.family[shouldBecomeOther]<- "Other"

dat <- dat %>%
  mutate(taxa.family = factor(taxa.family, levels = family_order))


dat %>%
  ggplot(aes(x = sample_id, y = rel_ab)) + 
  facet_grid(cols = vars(Clusters), scales = 'free_x', space = 'free_y') +
  geom_bar(aes(fill = taxa.family), stat = 'identity', position = 'fill', width = 1) +
  scale_fill_brewer(palette = 'Paired') +
  scale_y_continuous(name = 'Relative abundance',
                     labels = scales::percent) + 
  theme(axis.text.x = element_text(angle = 90, size = 6),
        axis.title.y = element_text(color = 'black'),
        strip.text = element_text(face = 'bold'),
        strip.background = element_blank())
















Tab = read.csv('results/Cluster_RNA_Cytof_tab.csv', sep = '\t', row.names = 1)
Tab$feature <-rownames(Tab)

### read or pass the taxa information
Meta = read.csv('results/Cluster_RNA_Cytof_meta.csv', sep = '\t', row.names = 1)
Meta = Meta[,c("Donor", "Clusters")]

dat <- Tab %>%
  pivot_longer(-feature, names_to = "Donor", values_to = 'per_of_cells')
dat$Donor <-gsub(".", " ", dat$Donor, fixed=TRUE)
dat <- dat %>%
  left_join(Meta, by='Donor')


colourCount = length(unique(Tab$feature))
getPalette = colorRampPalette(brewer.pal(9, "Paired"))



dat %>%
  ggplot(aes(x = Donor, y = per_of_cells)) + 
  facet_grid(cols = vars(Clusters), scales = 'free_x', space = 'free_y') +
  geom_bar(aes(fill = feature), stat = 'identity', position = 'fill', width = 1) +
  scale_fill_manual(values = getPalette(colourCount)) +
  scale_y_continuous(name = '% of cells',
                     labels = scales::percent) + 
  theme(axis.text.x = element_text(angle = 90, size = 6),
        axis.title.y = element_text(color = 'black'),
        strip.text = element_text(face = 'bold'),
        strip.background = element_blank())

library(mia)
library(miaTime)
library(miaViz)

data(GlobalPatterns, package="mia")
tse <- GlobalPatterns
tse <- transformAssay(tse, MARGIN = "samples", method="relabundance")

top_taxa <- getTopFeatures(tse,
                           method = "mean",
                           top = 5,
                           assay.type = "counts", useNames = TRUE)
top_taxa
