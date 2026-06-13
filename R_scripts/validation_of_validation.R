res_path <- "results/R_scripts_plots/"

library(tidyverse)
library(pheatmap)
library(grid)
library(gridExtra)
library(circlize)
library(pheatmap)
library(ggplot2)
library(ggpubr)

script_path_from_args <- function() {
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg) > 0) {
    return(normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = FALSE))
  }
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    active_path <- rstudioapi::getActiveDocumentContext()$path
    if (!is.null(active_path) && nzchar(active_path)) {
      return(normalizePath(active_path, mustWork = FALSE))
    }
  }
  NA_character_
}

script_path <- script_path_from_args()
script_dir <- if (!is.na(script_path)) dirname(script_path) else getwd()
repo_root <- if (basename(script_dir) == "R_scripts") dirname(script_dir) else script_dir
setwd(repo_root)
source(file.path(repo_root, "R_scripts", "cytometry_metafeatures_functions.R"))

set.seed(42)



### ----------------------------
### COLOR DEFINITIONS
### ----------------------------

cluster_colors <- c(
  "~Prevotella"      = "#e31a1c",
  "~Rumino&Christen" = "#FFA500",
  "~Bacteroidaceae"  = "#33a02c",
  "~Lachnospiraceae" = "#1f78b4"
)
ann_colors <- list(Cluster = cluster_colors)




### ----------------------------
### LOAD DATA
### ----------------------------

unstim_erc_cd4_base <- read.csv("data/5_per_file_MC_freq_CD4_ERC_neg.csv")

cd4_phenotype_summary <- read.csv("data/8_CD4_phenotype_summary.csv")

annot_metaclusters_cd4 <- read.csv("data/7_annotated_CD4_metaclusters.csv")

annot_metaclusters_cd4 <- annot_metaclusters_cd4 %>%
  rename(Metacluster = MC, Label = dictionary_label) %>%
  select(Metacluster, phenotype, Label) %>%
  mutate(
    phenotype = phenotype |>
      gsub("^lo ", "-",  x = _, fixed = TRUE) |>
      gsub("^hi ", "+",  x = _, fixed = TRUE) |>
      gsub("^lo",  "-",  x = _, fixed = TRUE) |>
      gsub("^hi",  "+",  x = _, fixed = TRUE)
  )

Meta <- read.csv("data/entero_clusters_cytometrycomparison.csv")

### ----------------------------
### PROCESS CD4
### ----------------------------

cd4_notjoined <- process_cytometry(unstim_erc_cd4_base) %>%
  filter(!str_detect(sample_id, "TYF")) %>%
  enforce_55_metaclusters() %>%
  remove_samples_percentage() %>%
  mutate(sample_time = paste(sample_id, Timepoint, sep = "_")) %>%
  filter(!sample_id %in% c("ERCVZV005.2", "ERCNICU06"))

cd4_clean <- inner_join(cd4_notjoined, Meta, by = "sample_id") %>% 
  select(-CMVIgG, -RunNumber) %>%
  rename(Enterotype = Cluster)




# --------------------------------------------------------------------
# 7. BOXPLOT — Check Metaclusters that were enriched in past analysis
# --------------------------------------------------------------------

cd4_phenotype_summary_long <- cd4_phenotype_summary %>%
  separate_rows(MCs, sep = ",") %>%       # split comma-separated values
  mutate(MCs = as.integer(MCs)) %>%
  rename(phenotype_mapped = phenotype)


cd4_clean_annot <- cd4_clean %>%
  inner_join(cd4_phenotype_summary_long, by = c("Metacluster" = "MCs")) %>%
  filter(!is.na(phenotype_mapped))

# Create folder if it doesn't exist
dir.create("results/R_scripts_plots/annotated_plots", showWarnings = FALSE)

etypes <- sort(unique(cd4_clean_annot$Enterotype))
comparisons <- combn(etypes, 2, simplify = FALSE)

cd4_sum <- cd4_clean_annot %>%
  group_by(sample_id, Enterotype, phenotype_mapped, Timepoint) %>%   # include Timepoint
  summarise(total_percentage = sum(percentage, na.rm = TRUE), .groups = "drop")


oldCytof <- read.csv("data/processed_data/Cytof_meta.csv", sep='\t', row.names = 1, check.names = FALSE)

oldEnterotypes <- read.csv("data/processed_data/chord_plot_data/Micro_clusters.csv", sep='\t', row.names = 1, check.names = FALSE)

oldEnterotypes <- oldEnterotypes %>%
  rename(sample_id = Donor) %>%
  rename(Enterotype = labels_Micro) %>%
  mutate(Enterotype = recode(Enterotype ,
                              "Bacteroidaceae driven" = "~Bacteroidaceae",
                              "Lachnospiraceae-driven" = "~Lachnospiraceae",
                              "Prevotella-driven\nHigh BMI" = "~Prevotella",
                              "Ruminococcaceae &\nChristensenellaceae\ndriven Low BMI" = "~Rumino&Christen"))

oldCytof_long <- oldCytof %>%
  pivot_longer(
    cols = -Donor,        # all columns except Donor
    names_to = "phenotype_mapped", 
    values_to = "total_percentage"
  ) %>%
  rename(sample_id = Donor) %>%
  mutate(phenotype_mapped = gsub("\\+", "_", phenotype_mapped)) %>%
  mutate(
    phenotype_mapped = gsub("Tna", "Tnaive", phenotype_mapped)
  ) %>%
  filter(phenotype_mapped %in% cd4_sum$phenotype_mapped)

oldCytof_long <- oldCytof_long %>%
  left_join(oldEnterotypes, by = "sample_id")



merged_df <- bind_rows(oldCytof_long, cd4_sum) %>%
  filter(!is.na(Enterotype)) %>%
  mutate(Timepoint = ifelse(is.na(Timepoint), "T1", Timepoint)) %>%  # fill NAs with T1
  filter(Timepoint == "T1")  


etypes <- sort(unique(merged_df$Enterotype))
comparisons <- combn(etypes, 2, simplify = FALSE)

# loop over phenotypes
for (p in unique(merged_df$phenotype_mapped)) {
  
  df_sub <- merged_df %>% filter(phenotype_mapped == p)
  
  p_plot <- ggplot(df_sub, aes(x = Enterotype, y = total_percentage, fill = Enterotype)) +
    geom_boxplot(outlier.alpha = 0, alpha = 0.6, width = 0.6, color = "black") +
    geom_jitter(aes(color = Enterotype), width = 0.15, size = 1.8, alpha = 0.8, show.legend = FALSE) +
    scale_fill_manual(values = cluster_colors) +
    scale_color_manual(values = cluster_colors) +
    ggtitle(paste("Phenotype:", p, "(Summed Percentages)")) +
    #stat_compare_means(
    #  comparisons = comparisons,
    #  method = "wilcox.test",
    #  label = "p.signif",
    #  hide.ns = TRUE,
    #  size = 4
    #) +
    theme_minimal(base_size = 12) +
    labs(
      x = "Enterotype Group",
      y = "Summed Percentage"
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      axis.text.y = element_text(size = 10),
      axis.title.x = element_blank(),
      axis.title.y = element_text(size = 12),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.7)
    )
  
  ggsave(
    filename = paste0("results/R_scripts_plots/annotated_plots/validation_summed_phenotype_", p, ".png"),
    plot = p_plot,
    width = 7,
    height = 6,
    dpi = 300
  )

}
