library(rstudioapi)

# Get the path of the script's directory
script_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)

# Set working directory one level above the script's directory
setwd(dirname(script_dir))


# Define the directory path
res_path <- "results/R_scripts_plots/"
data_path <- "data/processed_data/bubble_plot_data/"

# Check if the directory exists; if not, create it
if (!file.exists(res_path)) {
  dir.create(res_path, recursive = TRUE)
  cat("Directory created:", res_path, "\n")
} else {
  cat("Directory already exists:", res_path, "\n")
}


# Function to standardize columns
standardize <- function(x) {
  (x - mean(x)) / sd(x)
}


library(reshape2)
library(ggplot2)


############# TCR unimodal Bubble plot

# Define the matrix
tcrmodules_coefficients = as.matrix(read.csv(paste0(data_path,'TCR_bubbles_coefficients.csv'), sep = '\t', row.names = 1))
tcrmodules_pvalues = as.matrix(read.csv(paste0(data_path,'TCR_bubbles_pvalues.csv'), sep = '\t', row.names = 1))

# Convert the matrix to a dataframe
tcrmodules_coefficients <- as.data.frame(tcrmodules_coefficients)
# Set the row names
rownames(tcrmodules_coefficients) <- c("TCR C1", "TCR C2", "TCR C3", "TCR C4")
rownames(tcrmodules_pvalues) <- c("TCR C1", "TCR C2", "TCR C3", "TCR C4")


# Apply standardization to all columns except the first (assuming it's an ID or non-numeric)
standard_values <- as.data.frame(apply(tcrmodules_coefficients[, 1:6], 2, standardize))

tcrmodules_coefficients$Cluster <- rownames(tcrmodules_coefficients)
data_melted <- melt(tcrmodules_coefficients, id.vars = "Cluster")
standard_values$Cluster <- rownames(standard_values)
standard_values <- melt(standard_values, id.vars = "Cluster")

colnames(data_melted)[colnames(data_melted) == "value"] <- "log2FC"
data_melted$standard_log2fc <- standard_values$value

tcrmodules_pvalues <- as.data.frame(tcrmodules_pvalues)
tcrmodules_pvalues$Cluster <- rownames(tcrmodules_pvalues)
pvalues_melted <- melt(tcrmodules_pvalues, id.vars = "Cluster")
colnames(pvalues_melted)[colnames(pvalues_melted) == "value"] <- "pvalue_annotation"

# Merge the melted p-values with the melted coefficients dataframe
data_melted <- merge(data_melted, pvalues_melted, by = c("Cluster", "variable"))



p <- ggplot(data_melted, aes(x=Cluster, y = variable, color = standard_log2fc, size = abs(log2FC))) + 
  geom_point()+ 
  scale_size_continuous(range = c(2, 10)) +  # Set range of point sizes
  scale_colour_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, 
                         limit = c(min(data_melted$standard_log2fc), max(data_melted$standard_log2fc)), 
                         name = "Scaled log2 Fold Change") +
  cowplot::theme_cowplot() + 
  theme_minimal(base_size = 12) +  # `base_size` to make the plot smaller
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Clusters", y = "Features", title = "") +
  geom_text(aes(label = pvalue_annotation), color = 'black', vjust = -0.5, size = 6, angle = -45)  # Add p-value annotations

ggsave(paste0(res_path,'TCR_unimodal_bubble.png'), p, scale = 1, width = 8, height = 5, dpi = 600)





############# C&R Bubble plot
gsea_bacteria = read.csv(paste0(data_path,'bubbleplot_C&R 1.csv'), sep = '\t', row.names = 1)
gsea_bacteria$Cluster <- "C&R1"
aux = read.csv(paste0(data_path,'bubbleplot_C&R 2.csv'), sep = '\t', row.names = 1)
aux$Cluster <- "C&R2"
gsea_bacteria <- bind_rows(gsea_bacteria, aux)
aux = read.csv(paste0(data_path,'bubbleplot_C&R 3.csv'), sep = '\t', row.names = 1)
aux$Cluster <- "C&R3"
gsea_bacteria <- bind_rows(gsea_bacteria, aux)
aux = read.csv(paste0(data_path,'bubbleplot_C&R 4.csv'), sep = '\t', row.names = 1)
aux$Cluster <- "C&R4"
gsea_bacteria <- bind_rows(gsea_bacteria, aux)
aux = read.csv(paste0(data_path,'bubbleplot_C&R 5.csv'), sep = '\t', row.names = 1)
aux$Cluster <- "C&R5"
gsea_bacteria <- bind_rows(gsea_bacteria, aux)

p <- ggplot(gsea_bacteria, aes(y = Term, x = Cluster, color = NES, size = -log10(FDR.q.val) )) +
  geom_point() +
  scale_size_continuous(range = c(1, 5)) +  # Adjust size range as needed
  scale_colour_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +  # Adjust color gradient
  labs(
    x = "Cluster",
    y = "Bacterial Family",
    title = ""
  ) +
  geom_text(data = subset(gsea_bacteria, FDR.q.val < 0.25),  # Annotate points with FDR < 0.25
            aes(label = round(FDR.q.val, 2)),                # Format FDR to 2 decimal places
            hjust = 0.4, vjust = -2,                 # Adjust text position
            size = 5) +
  theme_minimal()  # Adjust theme as needed

ggsave(paste0(res_path,'CnR_GSEAbacteria_bubble.png'), p, scale = 1, width = 8, height = 5, dpi = 600)




############# Microbiome unimodal Bubble plot
gsea_bacteria = read.csv(paste0(data_path,'bubbleplot_Prevotella-driven_High_BMI.csv'), sep = '\t', row.names = 1)
gsea_bacteria$Cluster <- "Prevotella-driven\nHigh BMI"
aux = read.csv(paste0(data_path,'bubbleplot_Bacteroidaceae_driven.csv'), sep = '\t', row.names = 1)
aux$Cluster <- "Bacteroidaceae-driven"
gsea_bacteria <- bind_rows(gsea_bacteria, aux)
aux = read.csv(paste0(data_path,'bubbleplot_Ruminococcaceae_&_Christensenellaceae_driven_Low_BMI.csv'), sep = '\t', row.names = 1)
aux$Cluster <- "Ruminococcaceae &\nChristensenellaceae driven\nLow BMI"
gsea_bacteria <- bind_rows(gsea_bacteria, aux)
aux = read.csv(paste0(data_path,'bubbleplot_Lachnospiraceae-driven.csv'), sep = '\t', row.names = 1)
aux$Cluster <- "Lachnospiraceae-driven"
gsea_bacteria <- bind_rows(gsea_bacteria, aux)
gsea_bacteria$annot <- ifelse(gsea_bacteria$FDR.q.val > 0, round(gsea_bacteria$FDR.q.val, 4), "< 0.0001")
gsea_bacteria$size <- gsea_bacteria$FDR.q.val
gsea_bacteria$size[gsea_bacteria$size == 0] <- 0.0001


p <- ggplot(gsea_bacteria, aes(y = Term, x = Cluster, color = NES, size = -log10(size) )) +
  geom_point() +
  scale_size_continuous(range = c(1, 5)) +  # Adjust size range as needed
  scale_colour_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +  # Adjust color gradient
  labs(
    x = "Cluster",
    y = "Bacterial Family",
    title = ""
  ) +
  geom_text(data = subset(gsea_bacteria, FDR.q.val < 0.25),  # Annotate points with FDR < 0.25
            aes(label = annot),                # Format FDR to 2 decimal places
            hjust = 0.4, vjust = -2,                 # Adjust text position
            size = 5) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(paste0(res_path,'uniM_GSEAbacteria_bubble.png'), p, scale = 1, width = 8, height = 5, dpi = 600)



############# M&R Bubble plot
gsea_bacteria = read.csv(paste0(data_path,'bubbleplot_M&R Upregulated metabolic sets.csv'), sep = '\t', row.names = 1)
gsea_bacteria$Cluster <- "M&R Upregulated\nmetabolic sets"
aux = read.csv(paste0(data_path,'bubbleplot_M&R CMV seropositive.csv'), sep = '\t', row.names = 1)
aux$Cluster <- "M&R CMV seropositive"
gsea_bacteria <- bind_rows(gsea_bacteria, aux)
aux = read.csv(paste0(data_path,'bubbleplot_M&R Higher BMI (Prevotella).csv'), sep = '\t', row.names = 1)
aux$Cluster <- "M&R Higher BMI\n(Prevotella)"
gsea_bacteria <- bind_rows(gsea_bacteria, aux)
aux = read.csv(paste0(data_path,'bubbleplot_M&R Ruminococcaceae & Christensenellaceae - driven.csv'), sep = '\t', row.names = 1)
aux$Cluster <- "M&R\nRuminococcaceae &\nChristensenellaceae - driven "
gsea_bacteria <- bind_rows(gsea_bacteria, aux)
gsea_bacteria$annot <- ifelse(gsea_bacteria$FDR.q.val > 0, round(gsea_bacteria$FDR.q.val, 4), "< 0.0001")
gsea_bacteria$size <- gsea_bacteria$FDR.q.val
gsea_bacteria$size[gsea_bacteria$size == 0] <- 0.0001


p <- ggplot(gsea_bacteria, aes(y = Term, x = Cluster, color = NES, size = -log10(size) )) +
  geom_point() +
  scale_size_continuous(range = c(1, 5)) +  # Adjust size range as needed
  scale_colour_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +  # Adjust color gradient
  labs(
    x = "Cluster",
    y = "Bacterial Family",
    title = ""
  ) +
  geom_text(data = subset(gsea_bacteria, FDR.q.val < 0.25),  # Annotate points with FDR < 0.25
            aes(label = annot),                # Format FDR to 2 decimal places
            hjust = 0.4, vjust = -2,                 # Adjust text position
            size = 5) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(paste0(res_path,'MnR_GSEAbacteria_bubble.png'), p, scale = 1, width = 8, height = 5, dpi = 600)

