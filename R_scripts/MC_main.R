res_path <- "results/R_scripts_plots/"

library(tidyverse)
library(pheatmap)
library(grid)
library(gridExtra)
library(circlize)
library(pheatmap)
library(ggplot2)
library(ggpubr)

source("/Users/fabioaffaticati/Desktop/Desktop\ -\ Fabio’s\ MacBook\ Pro/Work/activ_covid-tcell-omics/R_scripts/cytometry_metafeatures_functions.R")
setwd("/Users/fabioaffaticati/Desktop/Desktop - Fabio’s MacBook Pro/Work/activ_covid-tcell-omics")

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
unstim_erc_cd8_base <- read.csv("data/5_per_file_MC_freq_CD8.csv")

cd4_phenotype_summary <- read.csv("data/8_CD4_phenotype_summary.csv")

annot_metaclusters_cd4 <- read.csv("data/7_annotated_CD4_metaclusters.csv")
annot_metaclusters_cd8 <- read.csv("data/7_annotated_CD8_metaclusters.csv")

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

annot_metaclusters_cd8 <- annot_metaclusters_cd8 %>%
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
### PROCESS CD4 & CD8
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

cd8_notjoined <- process_cytometry(unstim_erc_cd8_base) %>%
  filter(!str_detect(sample_id, "TYF")) %>%
  enforce_55_metaclusters() %>%
  remove_samples_percentage() %>%
  mutate(sample_time = paste(sample_id, Timepoint, sep = "_")) %>%
  filter(!sample_id %in% c("ERCVZV005.2", "ERCNICU06"))

cd8_clean <- inner_join(cd8_notjoined, Meta, by = "sample_id") %>% 
  select(-CMVIgG, -RunNumber) %>%
  rename(Enterotype = Cluster)

### ----------------------------
### ADD PHENOTYPES TO THE METACLUSTERS 
### ----------------------------

cd4_clean <- cd4_clean %>%
  left_join(annot_metaclusters_cd4, by = "Metacluster")
cd8_clean <- cd8_clean %>%
  left_join(annot_metaclusters_cd8, by = "Metacluster")



# -------------------------------------------------
# 1. Prepare SAMPLE-LEVEL matrix
# -------------------------------------------------

sample_meta_matrix <- cd4_clean %>%
  select(sample_time, Enterotype, percentage, phenotype) %>%
  pivot_wider(names_from = phenotype,
              values_from = percentage,
              values_fill = 0)

# extract matrix
mat <- as.matrix(sample_meta_matrix %>% select(-sample_time, -Enterotype))
rownames(mat) <- sample_meta_matrix$sample_time

# -------------------------------------------------
# 2. Feature-scaled version (z-score per feature)
# -------------------------------------------------

mat_scaled <- scale(mat)

# -------------------------------------------------
# 3. Annotation (Enterotype + Timepoint)
# -------------------------------------------------

row_annotation <- data.frame(
  Enterotype = sample_meta_matrix$Enterotype,
  Timepoint  = cd4_clean$Timepoint[match(sample_meta_matrix$sample_time, cd4_clean$sample_time)]
)
rownames(row_annotation) <- sample_meta_matrix$sample_time

# Define colors for Enterotype and Timepoint
timepoint_colors <- c(T1 = "#66c2a5", T2 = "#fc8d62")  # example: green/orange
ann_colors <- list(
  Enterotype = cluster_colors,
  Timepoint  = timepoint_colors
)

# -------------------------------------------------
# 4. HEATMAP — Scaled values (blue-white-red) with two annotations
# -------------------------------------------------

png(paste0(res_path, "heatmap_scaled_timepoint_cd4.png"), width = 2200, height = 2200, res = 180)
pheatmap(
  mat_scaled,
  annotation_row = row_annotation,
  annotation_colors = ann_colors,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  border_color = "gray",
  color = colorRampPalette(c("blue", "white", "red"))(100),
  main = "CD4 Metacluster Percentages (Z-score, Feature-level)"
)
dev.off()




# -------------------------------------------------
# 1. Prepare SAMPLE-LEVEL matrix
# -------------------------------------------------

sample_meta_matrix <- cd8_clean %>%
  select(sample_time, Enterotype, phenotype, percentage) %>%
  pivot_wider(names_from = phenotype,
              values_from = percentage,
              values_fill = 0)

# extract matrix
mat <- as.matrix(sample_meta_matrix %>% select(-sample_time, -Enterotype))
rownames(mat) <- sample_meta_matrix$sample_time

# -------------------------------------------------
# 2. Feature-scaled version (z-score per feature)
# -------------------------------------------------

mat_scaled <- scale(mat)

# -------------------------------------------------
# 3. Annotation (Enterotype + Timepoint)
# -------------------------------------------------

row_annotation <- data.frame(
  Enterotype = sample_meta_matrix$Enterotype,
  Timepoint  = cd8_clean$Timepoint[match(sample_meta_matrix$sample_time, cd8_clean$sample_time)]
)
rownames(row_annotation) <- sample_meta_matrix$sample_time

# Define colors for Enterotype and Timepoint
timepoint_colors <- c(T1 = "#66c2a5", T2 = "#fc8d62")  # example: green/orange
ann_colors <- list(
  Enterotype = cluster_colors,
  Timepoint  = timepoint_colors
)

# -------------------------------------------------
# 4. HEATMAP — Scaled values (blue-white-red) with two annotations
# -------------------------------------------------

png(paste0(res_path, "heatmap_scaled_timepoint_cd8.png"), width = 2200, height = 2200, res = 180)
pheatmap(
  mat_scaled,
  annotation_row = row_annotation,
  annotation_colors = ann_colors,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  border_color = "gray",
  color = colorRampPalette(c("blue", "white", "red"))(100),
  main = "CD8 Metacluster Percentages (Z-score, Feature-level)"
)
dev.off()



# -------------------------------------------------
# 5. Boxplots T1 vs T2 per Metacluster
# -------------------------------------------------

# Ensure Timepoint is a factor
cd8_notjoined$Timepoint <- factor(cd8_notjoined$Timepoint, levels = c("T1", "T2"))

cd8_notjoined <- cd8_notjoined %>%
  left_join(annot_metaclusters_cd8, by = "Metacluster")

# Boxplot with statistical comparison
p_timepoint_stats <- cd8_notjoined %>%
  ggplot(aes(x = Timepoint, y = percentage, fill = Timepoint)) +
  geom_boxplot(outlier.alpha = 0.2) +
  facet_wrap(~ Label, scales = "free_y") +
  scale_fill_manual(values = c(T1 = "#66c2a5", T2 = "#fc8d62")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Metacluster Percentage: T1 vs T2 (CD8) all data") +
  stat_compare_means(
    aes(group = Timepoint),
    method = "wilcox.test",           # or "wilcox.test"
    label = "p.format"
  )

# Save the plot
png(paste0(res_path, "boxplot_metaclusters_T1_vs_T2_cd8.png"), width = 2200, height = 2200, res = 150)
print(p_timepoint_stats)
dev.off()


# Ensure Timepoint is a factor
cd4_notjoined$Timepoint <- factor(cd4_notjoined$Timepoint, levels = c("T1", "T2"))

cd4_notjoined <- cd4_notjoined %>%
  left_join(annot_metaclusters_cd8, by = "Metacluster")



# Boxplot with statistical comparison
p_timepoint_stats <- cd4_notjoined %>%
  ggplot(aes(x = Timepoint, y = percentage, fill = Timepoint)) +
  geom_boxplot(outlier.alpha = 0.2) +
  facet_wrap(~ Label, scales = "free_y") +
  scale_fill_manual(values = c(T1 = "#66c2a5", T2 = "#fc8d62")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Metacluster Percentage: T1 vs T2 (CD4) all data") +
  stat_compare_means(
    aes(group = Timepoint),
    method = "wilcox.test",           # or "wilcox.test"
    label = "p.format"
  )

# Save the plot
png(paste0(res_path, "boxplot_metaclusters_T1_vs_T2_cd4.png"), width = 2200, height = 2200, res = 150)
print(p_timepoint_stats)
dev.off()


# -------------------------------------------------
# 6. BOXPLOT — Metacluster percentages per Enterotype
# -------------------------------------------------

p_box <- cd8_clean %>%
  ggplot(aes(x = Enterotype, y = percentage, fill = Enterotype)) +
  geom_boxplot(outlier.alpha = 0.2) +
  geom_point(size=.8, alpha=.7) +
  facet_wrap(~ Label, scales = "free_y") +
  scale_fill_manual(values = cluster_colors) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Metacluster Percentage by Enterotype (CD8)")

png(paste0(res_path, "boxplot_metacluster_enterotype_cd8.png"), width = 2200, height = 2200, res = 150)
print(p_box)
dev.off()

p_box <- cd4_clean %>%
  ggplot(aes(x = Enterotype, y = percentage, fill = Enterotype)) +
  geom_boxplot(outlier.alpha = 0.2) +
  geom_point(size=.8, alpha=.7) +
  facet_wrap(~ Label, scales = "free_y") +
  scale_fill_manual(values = cluster_colors) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Metacluster Percentage by Enterotype (CD4)")

png(paste0(res_path, "boxplot_metacluster_enterotype_cd4.png"), width = 2200, height = 2200, res = 150)
print(p_box)
dev.off()



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

# loop over timepoints
for (tp in unique(cd4_sum$Timepoint)) {
  
  df_tp <- cd4_sum %>% filter(Timepoint == tp)
  
  # loop over phenotypes
  for (p in unique(df_tp$phenotype_mapped)) {
    
    df_sub <- df_tp %>% filter(phenotype_mapped == p)
    
    p_plot <- ggplot(df_sub, aes(x = Enterotype, y = total_percentage, fill = Enterotype)) +
      geom_boxplot(outlier.alpha = 0, alpha = 0.6, width = 0.6, color = "black") +
      geom_jitter(aes(color = Enterotype), width = 0.15, size = 1.8, alpha = 0.8, show.legend = FALSE) +
      scale_fill_manual(values = cluster_colors) +
      scale_color_manual(values = cluster_colors) +
      ggtitle(paste("Phenotype:", p, "(Summed Percentages) -", tp)) +
      stat_compare_means(
        comparisons = comparisons,
        method = "wilcox.test",
        label = "p.signif",
        hide.ns = TRUE,
        size = 4
      ) +
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
      filename = paste0("results/R_scripts_plots/annotated_plots/", tp, "_summed_phenotype_", p, ".png"),
      plot = p_plot,
      width = 7,
      height = 6,
      dpi = 300
    )
    
  }
}
