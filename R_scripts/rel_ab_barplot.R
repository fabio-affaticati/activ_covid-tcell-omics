library(rstudioapi)

# Get the path of the script's directory
script_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)

# Set working directory one level above the script's directory
setwd(dirname(script_dir))


# Define the directory path
res_path <- "results/R_scripts_plots/"

# Check if the directory exists; if not, create it
if (!file.exists(res_path)) {
  dir.create(res_path, recursive = TRUE)
  cat("Directory created:", res_path, "\n")
} else {
  cat("Directory already exists:", res_path, "\n")
}


library(tidyverse)
library(RColorBrewer)




### pass relative abundances
Tab = read.csv('data/processed_data/relative_abundances_forplot.csv', sep = '\t', row.names = 1)

### read or pass the taxa information
Tax = read.csv('data/microbiome_analysis/taxa.csv', sep = ',', row.names = 1)

### pass the cluster information
Meta = read.csv('data/processed_data/chord_plot_data/Micro_clusters.csv', sep = '\t', row.names = 1)
Meta <- Meta %>% rename(sample_id = Donor)

Tab$taxa.taxon_id <-rownames(Tab)


dat <- Tab %>%
  pivot_longer(-taxa.taxon_id, names_to = "sample_id", values_to = 'rel_ab')

dat <- dat %>%
  left_join(Tax, by='taxa.taxon_id')

dat$sample_id <-gsub(".", " ", dat$sample_id, fixed=TRUE)

dat <- dat %>%
  left_join(Meta, by='sample_id')



family_order <- c('Other',
                  'Bacteroidaceae', 
                  'Lachnospiraceae', 
                  'Prevotellaceae', 
                  'Ruminococcaceae',
                  'Christensenellaceae')

shouldBecomeOther <- !(dat$taxa.family %in% family_order)
dat$taxa.family[shouldBecomeOther]<- "Other"

dat <- dat %>%
  mutate(taxa.family = factor(taxa.family, levels = family_order))

p <-  ggplot(dat, aes(x = sample_id, y = rel_ab)) + 
  facet_grid(cols = vars(labels_Micro), scales = 'free_x', space = 'free_y') +
  geom_bar(aes(fill = taxa.family), stat = 'identity', position = 'fill', width = 1) +
  scale_fill_brewer(palette = 'Paired', name = "Taxa family") +
  scale_y_continuous(name = 'Relative Abundance (%)',
                     labels = scales::percent) + 
  theme_bw() +
  # Modify text settings
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y = element_text(color = 'black', size = 12),
        strip.text = element_text(face = 'bold', size = 10),  # Increase facet title size
        strip.background = element_blank()) +
  xlab('Samples')

ggsave(paste0(res_path, "relative_abundance_microbiome_barplot.png"), plot = p, width = 10, height = 6, units = "in", dpi = 600)
