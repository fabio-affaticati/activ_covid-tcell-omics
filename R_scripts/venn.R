library(rstudioapi)

# Get the path of the script's directory
script_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)

# Set working directory one level above the script's directory
setwd(dirname(script_dir))


# Define the directory path
res_path <- "results/R_scripts_plots/"
data_path <- "data/processed_data/chord_plot_data/"

# Check if the directory exists; if not, create it
if (!file.exists(res_path)) {
  dir.create(res_path, recursive = TRUE)
  cat("Directory created:", res_path, "\n")
} else {
  cat("Directory already exists:", res_path, "\n")
}


library(ggplot2)
library(ggvenn)
library(reshape2)
library(dplyr)



############### Venn plot
RNAclusters <- read.csv(paste0(data_path, 'RNA_clusters.csv'), sep = '\t', row.names = 1)
Microclusters <- read.csv(paste0(data_path, 'Micro_clusters.csv'), sep = '\t', row.names = 1)
Cytofclusters <- read.csv(paste0(data_path, 'Cytof_clusters.csv'), sep = '\t', row.names = 1)
TCRclusters <- read.csv(paste0(data_path, 'TCR_clusters.csv'), sep = '\t', row.names = 1)


venn_data <- list(
  RNAseq = RNAclusters$Donor, 
  Microbiome = Microclusters$Donor, 
  CyTOF = Cytofclusters$Donor,
  TCRseq = TCRclusters$Donor
)



venn_plot <- ggvenn(venn_data, fill_color = c("orange", "darkgreen", "lightblue", "darkred"),
       stroke_size = 0.5, set_name_size = 4)

ggsave(paste0(res_path,"vennplot.png"), scale = 1, plot = venn_plot, width = 8, height = 6)
