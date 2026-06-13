# Define the directory path

res_path <- "results/R_scripts_plots/"

# Check if the directory exists; if not, create it
if (!file.exists(res_path)) {
  dir.create(res_path, recursive = TRUE)
  cat("Directory created:", res_path, "\n")
} else {
  cat("Directory already exists:", res_path, "\n")
}

setwd("/Users/fabioaffaticati/Desktop/Desktop\ -\ Fabio’s\ MacBook\ Pro/Work/activ_covid-tcell-omics")

library(tidyverse)
library(RColorBrewer)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(vegan)
library(gridExtra)
library(kernlab)
library(tidyverse)
library(cluster)   
library(mclust)    
library(readxl)
library(ggExtra)
library(gt)

set.seed(42)








extra_meta <- read_excel("data/Extra_microbiome_samples/benson/metadata.xlsx")

Meta = read.csv('metadata_validation.csv', sep = ',', row.names = 1)
Meta$sample_id <-gsub("-", ".", Meta$sample_id, fixed=TRUE)
colnames(Meta)[colnames(Meta) == "cluster"] <- "Cluster"

extra_meta <- extra_meta %>%
      rename(sample_id = Subjectnr)
extra_meta$sample_id <- str_replace_all(extra_meta$sample_id, " ", "")

extra_meta <- merge(extra_meta, Meta, by = "sample_id", all = TRUE)
sum(!is.na(extra_meta$CMVIgG))
extra_meta$CMVIgG <- str_replace_all(extra_meta$CMVIgG, "<", "")
extra_meta$CMVIgG <- str_replace_all(extra_meta$CMVIgG, ">", "")
extra_meta$CMVIgG[extra_meta$CMVIgG == "Unknown"] <- NA
sum(!is.na(extra_meta$CMVIgG))


Meta_encoded = read.csv('metadata_encoded_validation.csv', sep = ',', row.names = 1)
Meta_encoded$sample_id <-gsub("-", ".", Meta$sample_id, fixed=TRUE)


### pass relative abundances
Tab = read.csv('microbiome_abundances_validation.csv', sep = ',', row.names = 1)

### read or pass the taxa information
Tax = read.csv('microbiome_taxa_validation.csv', sep = ',', row.names = 1)

Tab$mapping <-rownames(Tab)


dat <- Tab %>%
  pivot_longer(-mapping, names_to = "sample_id", values_to = 'rel_ab')

dat <- dat %>%
  left_join(Tax, by='mapping')


### Not needed anymore after old taxonomy rework
#dat <- dat %>%
#  mutate(Family = ifelse(Genus == "Prevotella", "Prevotellaceae", Family))



family_order <- c('Other',
                  'Bacteroidaceae', 
                  'Lachnospiraceae', 
                  'Prevotellaceae', 
                  'Ruminococcaceae',
                  'Christensenellaceae')



shouldBecomeOther <- !(dat$Family %in% family_order)
dat$Family[shouldBecomeOther]<- "Other"

dat <- dat %>%
  mutate(Family = factor(Family, levels = family_order))

p <-  ggplot(dat, aes(x = sample_id, y = rel_ab)) + 
  geom_bar(aes(fill = Family), stat = 'identity', position = 'fill', width = 1) +
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
p
ggsave(paste0(res_path, "relative_abundance_microbiome_barplot_validation.png"), plot = p, width = 10, height = 6, units = "in", dpi = 600)




Tab$mapping <- NULL
# Step 1: Transform abundance table
Tab_hel <- decostand(Tab, method = "hellinger")

# Step 2: Distance matrix
bc_micro <- vegdist(t(Tab_hel), method = "bray")

# Step 3: PCoA
pcoa_result <- cmdscale(bc_micro, eig = TRUE, k = 2)
pcoa_points <- as.data.frame(pcoa_result$points)
colnames(pcoa_points) <- c("PCoA1", "PCoA2")
pcoa_points$SampleID <- rownames(pcoa_points)

# Step 4: Spectral clustering
sim_matrix <- exp(-as.matrix(bc_micro))
k <- 4
sc <- specc(as.kernelMatrix(sim_matrix), centers = k)

# Add cluster assignments
pcoa_points$cluster <- sc@.Data

# Add cluster assignments
pcoa_points$cluster <- as.factor(sc@.Data)
pcoa_points$cluster <- paste0("Cluster ", pcoa_points$cluster)

pcoa_points$cluster <- recode(pcoa_points$cluster,
                   "Cluster 2" = "~Prevotella",
                   "Cluster 3" = "~Lachnospiraceae",
                   "Cluster 1" = "~Rumino&Christen",
                   "Cluster 4" = "~Bacteroidaceae")


# Create metadata dataframe
meta <- data.frame(
  sample_id = colnames(Tab),         # samples from distance matrix
  cluster   = as.factor(sc@.Data)         # cluster assignments
)


p <- ggplot(pcoa_points, aes(x = PCoA1, y = PCoA2, color = cluster)) +
  geom_point(size = 3, alpha = 0.9) +
  labs(title = paste("PCoA (Bray-Curtis) with Spectral Clustering (k =", k, ")"),
       x = paste0("PCoA1 (", round(pcoa_result$eig[1] / sum(pcoa_result$eig) * 100, 1), "%)"),
       y = paste0("PCoA2 (", round(pcoa_result$eig[2] / sum(pcoa_result$eig) * 100, 1), "%)")) +
  theme_minimal() +
  scale_color_manual(
    values = c("~Prevotella"      = "#e31a1c",
               "~Rumino&Christen" = "#FFA500",
               "~Bacteroidaceae"  = "#33a02c",
               "~Lachnospiraceae" = "#1f78b4"),
    name = "Cluster"
  )
p
ggsave(paste0(res_path, "pcoa_spectralclustering_validation.png"), plot = p, width = 10, height = 6, units = "in", dpi = 600)





### pass relative abundances
Tab = read.csv('microbiome_abundances_validation.csv', sep = ',', row.names = 1)

Tab$mapping <-rownames(Tab)

 dat <- Tab %>%
  pivot_longer(-mapping, names_to = "sample_id", values_to = 'rel_ab')

dat <- dat %>%
  left_join(Tax, by='mapping')

# Join cluster assignments to your wide dataframe
dat <- dat %>%
  left_join(meta, by = c("sample_id" = "sample_id"))



dat <- dat %>%
  mutate(Family = ifelse(Genus == "Prevotella", "Prevotellaceae", Family))


family_order <- c('Other',
                  'Bacteroidaceae', 
                  'Lachnospiraceae', 
                  'Prevotellaceae', 
                  'Ruminococcaceae',
                  'Christensenellaceae')

shouldBecomeOther <- !(dat$Family %in% family_order)
dat$Family[shouldBecomeOther]<- "Other"

dat <- dat %>%
  mutate(Family = factor(Family, levels = family_order))

dat$cluster <- recode(dat$cluster,
                              "2" = "~Prevotella",
                              "3" = "~Lachnospiraceae",
                              "1" = "~Rumino&Christen",
                              "4" = "~Bacteroidaceae")

p <-  ggplot(dat, aes(x = sample_id, y = rel_ab)) + 
  facet_grid(cols = vars(cluster), scales = 'free_x', space = 'free_y') +
  geom_bar(aes(fill = Family), stat = 'identity', position = 'fill', width = 1) +
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
p
ggsave(paste0(res_path, "spectral_clustering_relative_abundance_microbiome_barplot_validation.png"), plot = p, width = 10, height = 6, units = "in", dpi = 600)








########################################################################
Meta = read.csv('metadata_validation.csv', sep = ',', row.names = 1)
Meta$sample_id <-gsub("-", ".", Meta$sample_id, fixed=TRUE)
colnames(Meta)[colnames(Meta) == "cluster"] <- "Cluster"
Meta_temp <- Meta
Meta_temp$Cluster <- recode(Meta_temp$Cluster,
                      "C 2" = "~Prevotella",
                      "C 3" = "~Lachnospiraceae",
                      "C 0" = "~Rumino&Christen",
                      "C 1" = "~Bacteroidaceae")
write_csv(Meta_temp, "data/entero_clusters_cytometrycomparison.csv")



dat <- dat %>%
  left_join(Meta, by = c("sample_id" = "sample_id"))

dat$Cluster <- recode(dat$Cluster,
                      "C 2" = "~Prevotella",
                      "C 3" = "~Lachnospiraceae",
                      "C 0" = "~Rumino&Christen",
                      "C 1" = "~Bacteroidaceae")


p <-  ggplot(dat, aes(x = sample_id, y = rel_ab)) + 
  facet_grid(cols = vars(Cluster), scales = 'free_x', space = 'free_y') +
  geom_bar(aes(fill = Family), stat = 'identity', position = 'fill', width = 1) +
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
p
ggsave(paste0(res_path, "snf_clustering_relative_abundance_microbiome_barplot_validation.png"), plot = p, width = 10, height = 6, units = "in", dpi = 600)



# Add cluster assignments
pcoa_points$Cluster <- as.factor(Meta$Cluster)
pcoa_points$Cluster <- recode(pcoa_points$Cluster,
                      "C 2" = "~Prevotella",
                      "C 3" = "~Lachnospiraceae",
                      "C 0" = "~Rumino&Christen",
                      "C 1" = "~Bacteroidaceae")

p <- ggplot(pcoa_points, aes(x = PCoA1, y = PCoA2, color = Cluster)) +
  geom_point(size = 3, alpha = 0.9) +
  labs(title = paste("PCoA (Bray-Curtis) with SNF + Spectral Clustering ( k =", k, ")"),
       x = paste0("PCoA1 (", round(pcoa_result$eig[1] / sum(pcoa_result$eig) * 100, 1), "%)"),
       y = paste0("PCoA2 (", round(pcoa_result$eig[2] / sum(pcoa_result$eig) * 100, 1), "%)")) +
  theme_minimal() +
  scale_color_manual(
    values = c("~Prevotella"      = "#e31a1c",
               "~Rumino&Christen" = "#FFA500",
               "~Bacteroidaceae"  = "#33a02c",
               "~Lachnospiraceae" = "#1f78b4"),
    name = "Cluster"
  )
p
ggsave(paste0(res_path, "pcoa_snf_spectralclustering_validation.png"), plot = p, width = 10, height = 6, units = "in", dpi = 600)






########################################################################
############### Eigentaxa encoded
Meta_encoded = read.csv('metadata_encoded_validation.csv', sep = ',', row.names = 1)
Meta_encoded$sample_id <-gsub("-", ".", Meta_encoded$sample_id, fixed=TRUE)

dat_encoded <- dat %>%
  left_join(Meta_encoded, by = c("sample_id" = "sample_id"))

dat_encoded$Cluster <- recode(dat_encoded$cluster_encoded,
                      "C 2" = "~Bacteroidaceae",
                      "C 3" = "~Lachnospiraceae",
                      "C 0" = "~Rumino&Christen",
                      "C 1" = "~Prevotella")


p <-  ggplot(dat_encoded, aes(x = sample_id, y = rel_ab)) + 
  facet_grid(cols = vars(Cluster), scales = 'free_x', space = 'free_y') +
  geom_bar(aes(fill = Family), stat = 'identity', position = 'fill', width = 1) +
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
p
ggsave(paste0(res_path, "snf_clustering_relative_abundance_encoded_microbiome_barplot_validation.png"), plot = p, width = 10, height = 6, units = "in", dpi = 600)


# Add cluster assignments
pcoa_points$Cluster_encoded <- as.factor(Meta_encoded$cluster_encoded)
pcoa_points$Cluster_encoded <- recode(pcoa_points$Cluster_encoded,
                                      "C 2" = "~Bacteroidaceae",
                                      "C 3" = "~Lachnospiraceae",
                                      "C 0" = "~Rumino&Christen",
                                      "C 1" = "~Prevotella")

p <- ggplot(pcoa_points, aes(x = PCoA1, y = PCoA2, color = Cluster)) +
  geom_point(size = 3, alpha = 0.9) +
  labs(title = paste("PCoA (Bray-Curtis) after encoding \n with SNF + Spectral Clustering ( k =", k, ")"),
       x = paste0("PCoA1 (", round(pcoa_result$eig[1] / sum(pcoa_result$eig) * 100, 1), "%)"),
       y = paste0("PCoA2 (", round(pcoa_result$eig[2] / sum(pcoa_result$eig) * 100, 1), "%)")) +
  theme_minimal() +
  scale_color_manual(
    values = c("~Prevotella"      = "#e31a1c",
               "~Rumino&Christen" = "#FFA500",
               "~Bacteroidaceae"  = "#33a02c",
               "~Lachnospiraceae" = "#1f78b4"),
    name = "Cluster"
  )
p
ggsave(paste0(res_path, "pcoa_snf_spectralclustering_encoded_validation.png"), plot = p, width = 10, height = 6, units = "in", dpi = 600)






























# --- Silhouette scores ---
# For clustering1
meta$cluster <- as.numeric(as.factor(meta$cluster))
sil1 <- silhouette(meta$cluster, bc_micro)
mean_sil1 <- mean(sil1[, 3])  # average silhouette width

# For clustering2
Meta$Cluster <- as.numeric(as.factor(Meta$Cluster))
sil2 <- silhouette(Meta$Cluster, bc_micro)
mean_sil2 <- mean(sil2[, 3])

# --- ARI score ---
ari <- adjustedRandIndex(meta$cluster, Meta$Cluster)

# Print results
cat("Average silhouette (clustering 1):", mean_sil1, "\n")
cat("Average silhouette (clustering 2):", mean_sil2, "\n")
cat("Adjusted Rand Index between clusterings:", ari, "\n")







dat_summary <- dat %>%
  filter(Family %in% c("Prevotellaceae",
                       "Ruminococcaceae",
                       'Christensenellaceae',
                       "Bacteroidaceae",
                       "Lachnospiraceae")) %>%
  group_by(sample_id, Family) %>%
  summarise(total_rel_ab = sum(rel_ab, na.rm = TRUE)) %>%
  ungroup()

# move rownames into a column
pcoa_points2 <- pcoa_points %>%
  rownames_to_column("sample_id")



# join cluster info onto df_summary
dat_summary1 <- dat_summary %>%
  left_join(pcoa_points2 %>% select(sample_id, cluster),
            by = "sample_id") %>%
  left_join(extra_meta, by = "sample_id")

bmis1 <- dat_summary1 %>%
  select(-Family, -total_rel_ab) %>%
  distinct() %>%
  mutate(BMI = round(as.numeric(BMI), 2))

boxplots <- ggplot(dat_summary1, aes(x = cluster, y = total_rel_ab, fill = cluster)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  facet_wrap(~ Family, scales = "free_y", nrow = 1) +  # all plots in one row
  labs(x = "Cluster", y = "Relative abundance (summed)", 
       title = "Relative Abundance per Family across Clusters") +
  theme_bw() +
  scale_fill_manual(
    values = c("~Prevotella"      = "#e31a1c",
               "~Rumino&Christen" = "#FFA500",
               "~Bacteroidaceae"  = "#33a02c",
               "~Lachnospiraceae" = "#1f78b4"),
    name = "Cluster"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  stat_compare_means(method = "kruskal.test", label = "p.format", size = 3)
boxplots
ggsave(paste0(res_path, "boxplots_spectralclustering_validation.png"), plot = boxplots, width = 12, height = 4, units = "in", dpi = 600)

boxplots_bmis <- ggplot(bmis1, aes(x = cluster, y = BMI, shape = cluster)) +
  geom_jitter(width = 0.2, alpha = 0.7, size = 2) +
  # per-cluster mean line (horizontal tick)
  stat_summary(fun = mean, geom = "point", 
               shape = 95, size = 15, color = "black") +
  labs(x = "Cluster", y = "BMI", 
       title = "BMI across Clusters") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"   # <- removes the legend
  ) +
  stat_compare_means(method = "wilcox.test", 
                     comparisons = combn(unique(bmis1$cluster), 2, simplify = FALSE),
                     label = "p.signif",
                     hide.ns = TRUE)
boxplots_bmis
ggsave(paste0(res_path, "boxplots_bmi_spectralclustering_validation.png"), plot = boxplots_bmis, width = 4, height = 5, units = "in", dpi = 600)




# join cluster info onto df_summary
dat_summary2 <- dat_summary %>%
  left_join(pcoa_points2 %>% select(sample_id, Cluster),
            by = "sample_id") %>%
  left_join(extra_meta %>% select(-Cluster), by = "sample_id")

bmis2 <- dat_summary2 %>%
  select(-Family, -total_rel_ab) %>%
  distinct() %>%
  mutate(BMI = round(as.numeric(BMI), 2))

dat_summary_summed <- dat_summary2 %>%
  mutate(Family = if_else(Family %in% c("Ruminococcaceae","Christensenellaceae"), "Rumi+Christ", Family))


boxplots <- ggplot(dat_summary_summed, aes(x = Cluster, y = total_rel_ab, fill = Cluster)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  facet_wrap(~ Family, scales = "free_y", nrow = 1) +
  labs(x = "Cluster", y = "Relative abundance (summed)", 
       title = "Relative Abundance per Family across Clusters (SNF)") +
  theme_bw() +
  scale_fill_manual(
    values = c("~Prevotella"      = "#e31a1c",
               "~Rumino&Christen" = "#FFA500",
               "~Bacteroidaceae"  = "#33a02c",
               "~Lachnospiraceae" = "#1f78b4"),
    name = "Cluster"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  stat_compare_means(method = "kruskal.test", label = "p.format", size = 3)

boxplots
ggsave(paste0(res_path, "boxplots_snf_spectralclustering_validation.png"), plot = boxplots, width = 12, height = 4, units = "in", dpi = 600)




bmis2 <- bmis2 %>%
  mutate(
    Gender = case_when(
      Gender == "m" ~ "Male",
      Gender == "f" ~ "Female",
      TRUE ~ Gender
    )
  )


boxplots_bmis <- ggplot(bmis2, aes(x = Cluster, y = BMI, shape = Cluster, color = Gender)) +
  geom_jitter(width = 0.2, alpha = 0.7, size = 2) +
  # per-cluster mean line (horizontal tick)
  stat_summary(fun = mean, geom = "point", 
               shape = 95, size = 15, color = "black") +
  labs(x = "Cluster", y = "BMI", 
       title = "BMI across Clusters (SNF)") +
  theme_bw() +
  guides(shape = "none") +   # <- hides shape legend only
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    #legend.position = "none"   # <- removes the legend
  ) +
  stat_compare_means(method = "wilcox.test", 
                     comparisons = combn(unique(bmis1$cluster), 2, simplify = FALSE),
                     label = "p.signif",
                     hide.ns = TRUE)
boxplots_bmis
ggsave(paste0(res_path, "boxplots_bmi_snf_spectralclustering_validation.png"), plot = boxplots_bmis, width = 4, height = 5, units = "in", dpi = 600)

bmis2$Age <- as.numeric(bmis2$Age)

boxplots_ages <- ggplot(bmis2, aes(x = Cluster, y = Age, shape = Cluster, color = Gender)) +
  geom_jitter(width = 0.2, alpha = 0.7, size = 2) +
  # per-cluster mean line (horizontal tick)
  stat_summary(fun = mean, geom = "point", 
               shape = 95, size = 15, color = "black") +
  labs(x = "Cluster", y = "Age", 
       title = "Age across Clusters (SNF)") +
  theme_bw() +
  guides(shape = "none") +   # <- hides shape legend only
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    #legend.position = "none"   # <- removes the legend
  ) +
  stat_compare_means(method = "wilcox.test", 
                     comparisons = combn(unique(bmis1$cluster), 2, simplify = FALSE),
                     label = "p.signif",
                     hide.ns = TRUE)
boxplots_ages
ggsave(paste0(res_path, "boxplots_age_snf_spectralclustering_validation.png"), plot = boxplots_ages, width = 4, height = 5, units = "in", dpi = 600)

bmis2$CMVIgG <- as.numeric(bmis2$CMVIgG)




plot_cmv <- ggplot(bmis2, aes(x = BMI, y = CMVIgG, color = Cluster)) +
  geom_point(alpha = 0.7, size = 2) +
  labs(
    x = "BMI",
    y = "CMVIgG (log10 scale)",
    title = "CMVIgG vs BMI across Clusters (SNF)"
  ) +
  theme_bw() +
  guides(shape = "none") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  scale_y_log10()

# ✅ Add marginal density plots on x and y axes
plot_cmv_marginals <- ggMarginal(
  plot_cmv,
  type = "density",   # can be "histogram", "boxplot", "violin", or "density"
  groupColour = TRUE, # color by Cluster
  groupFill = F
)

plot_cmv_marginals

ggsave(paste0(res_path, "plots_cmv_snf_spectralclustering_validation.png"), plot = plot_cmv_marginals, width = 6, height = 5, units = "in", dpi = 600)




###### eigentaxa encoded
dat_summary3 <- dat_summary %>%
  left_join(pcoa_points2 %>% select(sample_id, Cluster_encoded),
            by = "sample_id")

boxplots <- ggplot(dat_summary3, aes(x = Cluster_encoded, y = total_rel_ab, fill = Cluster_encoded)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  facet_wrap(~ Family, scales = "free_y", nrow = 1) +
  labs(x = "Cluster", y = "Relative abundance (summed)", 
       title = "Relative Abundance per Family across EncodedClusters (SNF)") +
  theme_bw() +
  scale_fill_manual(
    values = c("~Prevotella"      = "#e31a1c",
               "~Rumino&Christen" = "#FFA500",
               "~Bacteroidaceae"  = "#33a02c",
               "~Lachnospiraceae" = "#1f78b4"),
    name = "Cluster"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  stat_compare_means(method = "kruskal.test", label = "p.format", size = 3)

boxplots
ggsave(paste0(res_path, "boxplots_snf_spectralclustering_encoded_validation.png"), plot = boxplots, width = 12, height = 4, units = "in", dpi = 600)

















# apply to all clusters

# function to run one-vs-rest logistic regression
run_logreg <- function(cluster_label, data) {
  data <- data %>%
    mutate(y = ifelse(Cluster == cluster_label, 1, 0))
  
  fit <- glm(y ~ BMI + Age + CMVIgG, data = data, family = binomial)
  summary_fit <- summary(fit)
  
  # Helper function to safely extract coefficients
  safe_coef <- function(var, stat) {
    if (var %in% rownames(coef(summary_fit))) {
      coef(summary_fit)[var, stat]
    } else {
      NA
    }
  }
  
  data.frame(
    Cluster = cluster_label,
    Estimate_BMI = safe_coef("BMI", "Estimate"),
    Estimate_Age = safe_coef("Age", "Estimate"),
    Estimate_CMVIgG = safe_coef("CMVIgG", "Estimate"),
    p_value_BMI = safe_coef("BMI", "Pr(>|z|)"),
    p_value_Age = safe_coef("Age", "Pr(>|z|)"),
    p_value_CMVIgG = safe_coef("CMVIgG", "Pr(>|z|)")
  )
}



results <- lapply(levels(bmis2$Cluster), run_logreg, data = bmis2) %>%
  bind_rows() %>%
  mutate(p_adjust_BMI = p.adjust(p_value_BMI, method = "bonferroni")) %>% # Bonferroni correction
  mutate(p_adjust_Age = p.adjust(p_value_Age, method = "bonferroni")) %>% # Bonferroni correction
  mutate(p_adjust_CMVIgG = p.adjust(p_value_CMVIgG, method = "bonferroni")) # Bonferroni correction

print(results)



results_clean <- results %>%
  select(
    Cluster,
    Estimate_BMI, Estimate_Age, Estimate_CMVIgG,
    p_adjust_BMI, p_adjust_Age, p_adjust_CMVIgG
  ) %>%
  mutate(
    across(
      starts_with("Estimate_"),
      ~ ifelse(. > 0, "+", ifelse(. < 0, "-", NA))
    ),
    across(
      starts_with("p_adjust_"),
      ~ case_when(
        . <= 0.001 ~ "***",
        . <= 0.01  ~ "**",
        . <= 0.05  ~ "*",
        TRUE       ~ "ns"
      )
    )
  )


results_clean %>%
  gt() %>%
  tab_header(
    title = "Cluster-wise Logistic Regression Results",
    subtitle = "Signs of estimates (+/−) and significance levels (*, **, ***, ns)"
  ) %>%
  cols_label(
    Estimate_BMI = "BMI (Sign)",
    Estimate_Age = "Age (Sign)",
    Estimate_CMVIgG = "CMVIgG (Sign)",
    p_adjust_BMI = "BMI (Adj p)",
    p_adjust_Age = "Age (Adj p)",
    p_adjust_CMVIgG = "CMVIgG (Adj p)"
  ) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_body(columns = 2:4)
  ) %>%
  tab_options(
    table.font.size = 14,
    heading.background.color = "#f7f7f7",
    data_row.padding = px(4)
  )


write.csv(pcoa_points2, "metadata_for_entropy.csv")




# 1️⃣ Keep only numeric taxa columns
Tab_num <- Tab[, sapply(Tab, is.numeric)]

# 2️⃣ Calculate Shannon entropy
shannon <- diversity(t(Tab_num), index = "shannon")

# 3️⃣ Number of taxa present per sample
S <- specnumber(t(Tab_num))

# 4️⃣ Normalized Shannon entropy (0–1)
normalized_shannon <- shannon / log(S)

# 5️⃣ Combine into dataframe
normalized_shannon_df <- data.frame(
  Sample = rownames(t(Tab_num)),
  Normalized_Shannon = normalized_shannon,
  UnnormalizedShannon = shannon
)

# 6️⃣ Add Cluster info
normalized_shannon_df$Cluster <- pcoa_points2$Cluster

# 7️⃣ (Optional) If you later have more diversity metrics, pivot them to long format
# For now, we only have one metric — so we’ll keep it simple
div_long <- normalized_shannon_df %>%
  pivot_longer(
    cols = c(Normalized_Shannon), 
    names_to = "Metric", 
    values_to = "Value"
  )

# 8️⃣ Plot
# Ensure Cluster is a factor with the correct order
div_long$Cluster <- factor(
  div_long$Cluster,
  levels = c("~Lachnospiraceae", "~Rumino&Christen", "~Bacteroidaceae", "~Prevotella")
)

# Define pairwise comparisons
my_comparisons <- list(
  c("~Prevotella", "~Rumino&Christen"),
  c("~Prevotella", "~Bacteroidaceae"),
  c("~Prevotella", "~Lachnospiraceae"),
  c("~Rumino&Christen", "~Bacteroidaceae"),
  c("~Rumino&Christen", "~Lachnospiraceae"),
  c("~Bacteroidaceae", "~Lachnospiraceae")
)

# Compute dynamic lower limit
ymin <- min(div_long$Value, na.rm = TRUE) -0.1

# Plot
y_positions <- seq(.88, 1.15, length.out = length(my_comparisons))

# Plot
p <- ggplot(div_long, aes(x = Cluster, y = Value, fill = Cluster)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.6) +
  geom_jitter(width = 0.2, alpha = 0.6, size = 1.5) +
  #facet_wrap(~ Metric, scales = "free_y", nrow = 1) +
  labs(
    x = "Cluster",
    y = "Diversity Index",
    title = "Normalized Shannon Diversity across Clusters (SNF)"
  ) +
  theme_bw(base_size = 12) +
  scale_fill_manual(
    values = c("~Prevotella"      = "#e31a1c",
               "~Bacteroidaceae"  = "#33a02c",
               "~Rumino&Christen" = "#FFA500",
               "~Lachnospiraceae" = "#1f78b4"),
    name = "Cluster"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    panel.grid.minor = element_blank()
  ) +
  coord_cartesian(ylim = c(ymin, 1.15)) +
  stat_compare_means(
    method = "wilcox.test",
    comparisons = my_comparisons,
    label = "p.signif",
    hide.ns = T,
    label.y = y_positions
  )
p
ggsave(paste0(res_path, "plots_normshannon_snf_spectralclustering_validation.png"), plot = p, width = 6, height = 5, units = "in", dpi = 600)


