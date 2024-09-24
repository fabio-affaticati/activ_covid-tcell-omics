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
library(reshape2)
library(circlize)
library(dplyr)


############### RNA unimodal stability plot
png(paste0(res_path, "chord_diagram_RNAstability.png"), width = 8, height = 8, units = "in", res = 600)
chord = read.csv(paste0(data_path, 'RNA_chord.csv'), sep = '\t', row.names = 1)
circos.par(start.degree = 90)
chordDiagram(chord, big.gap = 30,
             order = unique(c(sort(chord$from, decreasing = TRUE), sort(chord$to))))
dev.off()
circos.clear()
#############################################



############################################
RNAclusters <- read.csv(paste0(data_path, 'RNA_clusters.csv'), sep = '\t', row.names = 1)
Microclusters <- read.csv(paste0(data_path, 'Micro_clusters.csv'), sep = '\t', row.names = 1)
Cytofclusters <- read.csv(paste0(data_path, 'Cytof_clusters.csv'), sep = '\t', row.names = 1)
TCRclusters <- read.csv(paste0(data_path, 'TCR_clusters.csv'), sep = '\t', row.names = 1)

merged_df <- merge(RNAclusters, Microclusters, by = "Donor", all = TRUE)
merged_df <- merge(merged_df, Cytofclusters, by = "Donor", all = TRUE)
merged_df <- merge(merged_df, TCRclusters, by = "Donor", all = TRUE)
merged_df[is.na(merged_df)] <- 0


# cytof vs rna
counts <- as.data.frame(table(merged_df$labels_RNA, merged_df$labels_Cytof))
names(counts) <- c("to", "from", "value")
counts <- counts %>%
  filter_all(all_vars(. != "0"))

png(paste0(res_path, "CytofRNA_chord_diagram.png"), width = 7, height = 7, units = "in", res = 600)
circos.par(start.degree = 90)
chordDiagram(counts, big.gap = 30, small.gap = 12,
             order = unique(c(sort(counts$to, decreasing = TRUE), sort(counts$from))))
dev.off()
circos.clear()


# microbiome vs rna
counts <- as.data.frame(table(merged_df$labels_RNA, merged_df$labels_Micro))
names(counts) <- c("to", "from", "value")
counts <- counts %>%
  filter_all(all_vars(. != "0"))
png(paste0(res_path, "MicroRNA_chord_diagram.png"), width = 7, height = 7, units = "in", res = 600)
circos.par(start.degree = 90)
chordDiagram(counts, big.gap = 30, small.gap = 5,
             order = unique(c(sort(counts$to, decreasing = TRUE), sort(counts$from))))
dev.off()
circos.clear()


# microbiome vs cytof
counts <- as.data.frame(table(merged_df$labels_Micro, merged_df$labels_Cytof))
names(counts) <- c("to", "from", "value")
counts <- counts %>%
  filter_all(all_vars(. != "0"))

png(paste0(res_path, "CytofMicro_chord_diagram.png"), width = 7, height = 7, units = "in", res = 600)
circos.par(start.degree = 90)
chordDiagram(counts, big.gap = 30, small.gap = 12,
             order = unique(c(sort(counts$from, decreasing = TRUE), sort(counts$to))))
dev.off()
circos.clear()
############################################







############### Multi chord plots
multi_chord = read.csv(paste0(data_path, 'multi_chord_cmr_mr.csv'), sep = '\t', row.names = 1)
multi_chord <- multi_chord %>%
  mutate(Micro.RNA = recode(Micro.RNA, 
                                  "M&R Upregulated metabolic sets" = "M&R Upregulated\nmetabolic sets",
                                  "M&R Ruminococcaceae & Christensenellaceae - driven" = "M&R\nRuminococcaceae &\nChristensenellaceae - driven"))

png(paste0(res_path, "Trimodal_MnR_chord_diagram.png"), width = 8, height = 8, units = "in", res = 600)
circos.par(start.degree = 90)
chordDiagram(multi_chord, big.gap = 30,
             order = unique(c(sort(multi_chord$Micro.RNA, decreasing = TRUE),
                              sort(multi_chord$Cytof.Micro.RNA))))
dev.off()
circos.clear()





multi_chord =  read.csv(paste0(data_path, 'multi_chord_cmr_cr.csv'), sep = '\t', row.names = 1)

png(paste0(res_path, "Trimodal_CnR_chord_diagram.png"), width = 8, height = 8, units = "in", res = 600)

circos.par(start.degree = 90)
chordDiagram(multi_chord, big.gap = 30,
             order = unique(c(sort(multi_chord$Cytof.RNA, decreasing = TRUE),
                              sort(multi_chord$Cytof.Micro.RNA))))
dev.off()
circos.clear()
