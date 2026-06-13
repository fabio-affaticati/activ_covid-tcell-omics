library(rstudioapi)

# Get the path of the script's directory
script_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)

# Set working directory one level above the script's directory
setwd(dirname(script_dir))

# Define the directory path
res_path  <- "results/R_scripts_plots/"
data_path <- "data/processed_data/"

# Check if the directory exists; if not, create it
if (!file.exists(res_path)) {
  dir.create(res_path, recursive = TRUE)
  cat("Directory created:", res_path, "\n")
} else {
  cat("Directory already exists:", res_path, "\n")
}

library(vegan)
library(gridExtra)
library(kernlab)
library(tidyverse)
library(compositions)   # for clr()

set.seed(42)

# =============================================================================
# LOAD DATA
# =============================================================================
micro_data <- read.csv(paste0(data_path, "microbial_unimodal_adonis.csv"), row.names = 1)
metadata   <- read.csv(paste0(data_path, "adonis_labels.csv"),             row.names = 1)
Tax        <- read.csv("data/microbiome_analysis/taxa.csv", sep = ",",     row.names = 1)

rownames(micro_data) <- metadata$Donor
metadata$labels      <- factor(metadata$labels)

# Custom colours — shared by both Bray-Curtis and Aitchison plots
custom_colors <- c(
  "Prevotella-driven\nHigh BMI"                            = "#e31a1c",
  "Bacteroidaceae driven"                                  = "#33a02c",
  "Ruminococcaceae &\nChristensenellaceae\ndriven Low BMI" = "#FFA500",
  "Lachnospiraceae-driven"                                 = "#1f78b4"
)

# =============================================================================
# SECTION 1: BRAY-CURTIS (original analysis)
# =============================================================================

bc_micro <- vegdist(micro_data)

# PERMDISP2
cat("\n--- PERMDISP2 (Bray-Curtis) ---\n")
print(anova(betadisper(bc_micro, metadata$labels)))

# PERMANOVA
permanova_bc <- adonis2(bc_micro ~ labels, data = metadata, by = "margin",
                        permutations = 1000, method = "braycurtis")
cat("\n--- PERMANOVA (Bray-Curtis) ---\n")
print(permanova_bc)
p_val_bc <- permanova_bc$`Pr(>F)`[1]
r2_bc    <- permanova_bc$R2[1]

metadata %>% dplyr::count(labels)

# PCoA
PCoA_bc <- vegan::betadisper(bc_micro, metadata$labels)
c_bc    <- cmdscale(bc_micro, eig = TRUE)

# Species correlations with PCoA axes
corr1_bc <- sapply(micro_data[], function(x) cor(PCoA_bc$vectors[, 1], x))
corr2_bc <- sapply(micro_data[], function(x) cor(PCoA_bc$vectors[, 2], x))

corr_matrix_bc              <- as.data.frame(cbind(corr1_bc, corr2_bc))
colnames(corr_matrix_bc)    <- c("corr1", "corr2")
corr_matrix_bc$x0           <- 0
corr_matrix_bc$y0           <- 0
corr_matrix_bc              <- corr_matrix_bc[
  order(abs(corr_matrix_bc$corr1), abs(corr_matrix_bc$corr2), decreasing = TRUE), ]
corr_matrix_bc$taxa.taxon_id <- rownames(corr_matrix_bc)
corr_matrix_bc <- corr_matrix_bc %>% left_join(Tax, by = "taxa.taxon_id")

to_plot_bc        <- as.data.frame(PCoA_bc$vectors[, 1], nm = c("PCoA1"))
to_plot_bc$PCoA2  <- PCoA_bc$vectors[, 2]
to_plot_bc$Cluster <- PCoA_bc$group

# --- Biplot (Bray-Curtis) ---
p_bc <- ggplot(to_plot_bc, aes(x = PCoA1, y = PCoA2, color = Cluster)) +
  geom_point(size = 4) +
  scale_color_manual(values = custom_colors) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
  labs(title = "") +
  labs(
    title = paste0("Bray-Curtis PCoA  |  PERMANOVA R²=", 
                   round(r2_bc, 3), "  p=<", round(p_val_bc, 5)),
    #x     = paste0("PCoA1 Exp Var = ", exp_var1_ait, "%"),
    #y     = paste0("PCoA2 Exp Var = ", exp_var2_ait, "%")
  ) +
  theme_minimal() +
  theme(
    plot.title   = element_text(size = 25),
    legend.title = element_text(size = 20),
    legend.text  = element_text(size = 16),
    axis.title   = element_text(size = 14)
  )


ggsave(paste0(res_path, "PCoAbiplot.png"), p_bc, scale = 1, width = 16, height = 12)

# --- Density-marginal PCoA (Bray-Curtis) ---
pos_eig_bc <- c_bc$eig[c_bc$eig > 0]
exp_var1_bc <- round(pos_eig_bc[1] / sum(pos_eig_bc), 2) * 100
exp_var2_bc <- round(pos_eig_bc[2] / sum(pos_eig_bc), 2) * 100

scat_bc <- ggplot(to_plot_bc, aes(x = PCoA1, y = PCoA2, color = Cluster)) +
  theme_bw() +
  geom_point() +
  scale_color_manual(values = custom_colors) +
  stat_ellipse(level = 0.95) +
  labs(
    x = paste0("PCoA1 Exp Var = ", exp_var1_bc, "%"),
    y = paste0("PCoA2 Exp Var = ", exp_var2_bc, "%")
  ) +
  geom_rug(aes(color = Cluster)) +
  theme(
    legend.position = "none",
    axis.ticks      = element_blank(),
    axis.text.x     = element_blank(),
    axis.text.y     = element_blank(),
    axis.title      = element_text(size = 16),
    axis.text       = element_text(size = 14),
    legend.text     = element_text(size = 12)
  )

xdensity_bc <- ggplot(to_plot_bc, aes(x = PCoA1, fill = Cluster)) +
  geom_density(alpha = .5) + theme_bw() +
  scale_fill_manual(values = custom_colors) +
  theme(legend.position = "none", axis.ticks = element_blank(),
        axis.text = element_blank(), axis.title = element_text(size = 16))

ydensity_bc <- ggplot(to_plot_bc, aes(y = PCoA2, fill = Cluster)) +
  geom_density(alpha = .5) + theme_bw() +
  scale_fill_manual(values = custom_colors) +
  theme(legend.position = "none", axis.ticks = element_blank(),
        axis.text = element_blank(), axis.title = element_text(size = 16))

blankPlot <- ggplot() + geom_blank(aes(1, 1)) +
  theme(plot.background = element_blank(),
        panel.grid = element_blank(), panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title = element_blank(), axis.text = element_blank(),
        axis.ticks = element_blank())

p_density_bc <- grid.arrange(xdensity_bc, blankPlot, scat_bc, ydensity_bc,
                             ncol = 2, nrow = 2,
                             widths = c(4, 1.4), heights = c(1.4, 4))
ggsave(paste0(res_path, "PCoA_density.png"),
       p_density_bc, scale = 1, width = 20, height = 12)

# =============================================================================
# SECTION 2: CLR + AITCHISON DISTANCE (sensitivity analysis)
# =============================================================================

# --- CLR transformation ---
micro_clr_mat <- as.data.frame(clr(micro_data))

# Aitchison distance = Euclidean distance on CLR-transformed data
aitchison_dist <- dist(micro_clr_mat, method = "euclidean")

# --- PERMDISP2 (Aitchison) ---
cat("\n--- PERMDISP2 (Aitchison) ---\n")
bd_aitchison <- betadisper(aitchison_dist, metadata$labels)
print(anova(bd_aitchison))

# --- PERMANOVA (Aitchison) ---
cat("\n--- PERMANOVA (Aitchison) ---\n")
permanova_aitchison <- adonis2(
  aitchison_dist ~ labels,
  data         = metadata,
  by           = "margin",
  permutations = 1000
)
print(permanova_aitchison)
p_val_aitchison <- permanova_aitchison$`Pr(>F)`[1]
r2_aitchison    <- permanova_aitchison$R2[1]

# --- PCoA via cmdscale (Aitchison) ---
cmd_ait      <- cmdscale(aitchison_dist, k = 3, eig = TRUE) # k=3 to get the 3rd dimension
cmd_coords   <- as.data.frame(cmd_ait$points[, 1:3])
colnames(cmd_coords) <- c("PCoA1", "PCoA2", "PCoA3")
cmd_coords$Cluster   <- metadata$labels

# --- NEW: Scree Plot (Variance Explained) ---
pos_eig_ait <- cmd_ait$eig[cmd_ait$eig > 0]
variances   <- pos_eig_ait / sum(pos_eig_ait) * 100

png(paste0(res_path, "Aitchison_ScreePlot.png"), width = 800, height = 600)
barplot(variances[1:10], 
        names.arg = paste0("PC", 1:10), 
        main = "Aitchison Variance Explained per Axis",
        ylab = "% Variance", col = "steelblue")
dev.off()

exp_var1_ait <- round(variances[1], 1)
exp_var2_ait <- round(variances[2], 1)
exp_var3_ait <- round(variances[3], 1)


# --- Species correlations (using first 2 axes for your standard biplot) ---
corr1_ait <- sapply(micro_clr_mat, function(x) cor(cmd_coords$PCoA1, x))
corr2_ait <- sapply(micro_clr_mat, function(x) cor(cmd_coords$PCoA2, x))

corr_matrix_ait              <- as.data.frame(cbind(corr1_ait, corr2_ait))
colnames(corr_matrix_ait)    <- c("corr1", "corr2")
corr_matrix_ait$x0           <- 0
corr_matrix_ait$y0           <- 0
corr_matrix_ait              <- corr_matrix_ait[
  order(abs(corr_matrix_ait$corr1), abs(corr_matrix_ait$corr2), decreasing = TRUE), ]
corr_matrix_ait$taxa.taxon_id <- rownames(corr_matrix_ait)
corr_matrix_ait <- corr_matrix_ait %>% left_join(Tax, by = "taxa.taxon_id")

# Scale arrows
scale_factor <- max(abs(range(cmd_coords$PCoA1)), abs(range(cmd_coords$PCoA2))) * 0.8
corr_matrix_ait$corr1_sc <- corr_matrix_ait$corr1 * scale_factor
corr_matrix_ait$corr2_sc <- corr_matrix_ait$corr2 * scale_factor

# --- Biplot (Aitchison) ---
p_ait <- ggplot(cmd_coords, aes(x = PCoA1, y = PCoA2, color = Cluster)) +
  geom_point(size = 4) +
  scale_color_manual(values = custom_colors) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
  coord_equal() +
  labs(
    title = paste0("Aitchison distance PCoA  |  PERMANOVA R²=", 
                   round(r2_aitchison, 3), "  p=<", round(p_val_aitchison, 5)),
    #x     = paste0("PCoA1 Exp Var = ", exp_var1_ait, "%"),
    #y     = paste0("PCoA2 Exp Var = ", exp_var2_ait, "%")
  ) +
  theme_minimal() +
  theme(
    plot.title   = element_text(size = 25),
    legend.title = element_text(size = 20),
    legend.text  = element_text(size = 15),
    axis.title   = element_text(size = 20)
  )

ggsave(paste0(res_path, "PCoAbiplot_Aitchison.png"),
       p_ait, scale = 1, width = 16, height = 12)



# =============================================================================
# SUMMARY TABLE — paste into reviewer response
# =============================================================================
cat("\n========================================\n")
cat("  SENSITIVITY ANALYSIS SUMMARY\n")
cat("========================================\n")
cat(sprintf("  Bray-Curtis  PERMANOVA:  R² = %.4f,  p = %.4f\n", r2_bc,         p_val_bc))
cat(sprintf("  Aitchison    PERMANOVA:  R² = %.4f,  p = %.4f\n", r2_aitchison,  p_val_aitchison))
cat("========================================\n")