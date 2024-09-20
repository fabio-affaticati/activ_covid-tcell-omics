library(circlize)
library(dplyr)

library(UpSetR)
library(ggplot2)

library(ggvenn)
library(reshape2)

workingDir = "/Users/fabioaffaticati/Desktop/Work/activ_covid-and-omics"
setwd(workingDir)

RNAclusters <- read.csv('RNA_clusters.csv', sep = '\t', row.names = 1)
Microclusters <- read.csv('Micro_clusters.csv', sep = '\t', row.names = 1)
Cytofclusters <- read.csv('Cytof_clusters.csv', sep = '\t', row.names = 1)
TCRclusters <- read.csv('TCR_clusters.csv', sep = '\t', row.names = 1)

merged_df <- merge(RNAclusters, Microclusters, by = "Donor", all = TRUE)
merged_df <- merge(merged_df, Cytofclusters, by = "Donor", all = TRUE)
merged_df <- merge(merged_df, TCRclusters, by = "Donor", all = TRUE)
merged_df[is.na(merged_df)] <- 0

# cytof vs rna
counts <- as.data.frame(table(merged_df$labels_RNA, merged_df$labels_Cytof))
names(counts) <- c("to", "from", "value")
counts$to <- paste("RNA", counts$to, sep = " ")
counts$from <- paste("Cytof", counts$from, sep = " ")
circos.par(start.degree = 90)
chordDiagram(counts, big.gap = 20, small.gap = 7,
             order = c("RNA C 4", "RNA C 3", "RNA C 2",
                       "RNA C 1", "RNA C 0", "Cytof C 0", "Cytof C 1",
                       "Cytof C 2", "Cytof C 3", "Cytof C 4", "Cytof C 5")
)

png("CytofRNA_circle_plot.png")
dev.off()

# microbiome vs rna
counts <- as.data.frame(table(merged_df$labels_RNA, merged_df$labels_Micro))
names(counts) <- c("to", "from", "value")
counts$to <- paste("RNA", counts$to, sep = " ")
counts$from <- paste("Micro", counts$from, sep = " ")
circos.par(start.degree = 90)
chordDiagram(counts, big.gap = 20, small.gap = 7,
             order = c("RNA C 4", "RNA C 3", "RNA C 2",
                       "RNA C 1", "RNA C 0", "Micro C 0", "Micro C 1",
                       "Micro C 2", "Micro C 3")
)

png("MicroRNA_circle_plot.png")
dev.off()






merged_df[merged_df != 0] <- 1
names(merged_df) <- c("Donor", "RNAseq", "Microbiome", "CyTOF", "TCRseq")
merged_df[, 2:ncol(merged_df)] <- lapply(merged_df[, 2:ncol(merged_df)] , as.numeric)
upset_plot <- upset(merged_df, sets = c("RNAseq", "Microbiome", "CyTOF", 'TCRseq'), order.by = 'freq', group.by = 'degree', show.numbers = TRUE,text.scale = 3, point.size = 5,
      sets.bar.color = c("orange", "darkgreen", "lightblue", "darkred"))
upset_plot



venn_data <- list(
  RNAseq = RNAclusters$Donor, 
  Microbiome = Microclusters$Donor, 
  CyTOF = Cytofclusters$Donor,
  TCRseq = TCRclusters$Donor
)



venn_plot <- ggvenn(venn_data, fill_color = c("orange", "darkgreen", "lightblue", "darkred"),
       stroke_size = 0.5, set_name_size = 4)
# Save the plot
ggsave("results/vennplot.png", scale = 1, plot = venn_plot, width = 8, height = 6)




workingDir = "/Users/fabioaffaticati/Desktop/Work/activ_covid-and-omics/results"
setwd(workingDir) 


chord = read.csv('chord.csv', sep = '\t', row.names = 1)
circos.par(start.degree = 90)
chordDiagram(chord, big.gap = 30, order = c("EC 3", "EC 2", "EC 1", "EC 0",
                                            "C 0", "C 1", "C 2", "C 3") 
              )


circos.clear()
cytof_vs_microbiome = read.csv('cytof_vs_microbiome.csv', sep = '\t', row.names = 1)
cytof_vs_microbiome <- cytof_vs_microbiome %>%
                          mutate(to = recode(to, 
                                                  "M_C 0" = "Prevotella-driven\nHigh BMI",
                                                  "M_C 1" = "Bacteroidaceae-driven",
                                                  "M_C 2" = "Ruminococcaceae &\nChristensenellaceae-driven\nLow BMI",
                                                  "M_C 3" = "Lachnospiraceae-driven"))
cytof_vs_microbiome <- cytof_vs_microbiome %>%
  mutate(from = recode(from, 
                     "C_C 0" = "LAG3 cluster",
                     "C_C 1" = "Monocyte-driven",
                     "C_C 2" = "CMV seropositive\nsenescent",
                     "C_C 3" = "FAS cluster",
                     "C_C 4" = "Naive group",
                     "C_C 5" = "Elderly group",))




circos.par(start.degree = 90)
par(cex = 2, mar = c(0, 0, 0, 0))
chordDiagram(cytof_vs_microbiome, big.gap = 30,
             order = c("Elderly group","Naive group","FAS cluster", "CMV seropositive\nsenescent",  "Monocyte-driven", "LAG3 cluster",
                       "Prevotella-driven\nHigh BMI", "Bacteroidaceae-driven", "Ruminococcaceae &\nChristensenellaceae-driven\nLow BMI", "Lachnospiraceae-driven")
)
annotationTrackHeight = c(0.1, 0.01)
title(main = "Microbiome (left) and Cytof (right) cluster 'conservation'")

png("Microbiome_vs_Cytof.png")
dev.off()














randindexdata = read.csv('randindexpermu.csv', sep = '\t', row.names = 1)


# from https://rdrr.io/github/fruciano/GeometricMorphometricsMix/src/R/adjRand_test.R
adjRand_test=function(A, B, perm=1000) {
  if (length(A)!=length(B)) { stop("A and B should have the same length") }
  # Make sure that the two groups of partitions have the same length
  
  ARIo=mclust::adjustedRandIndex(A, B)
  # Observed adjusted Rand index
  print(ARIo)
  
  Aperm=lapply(seq_len(perm), function(X) sample(A, length(A), replace=FALSE))
  Bperm=lapply(seq_len(perm), function(X) sample(B, length(B), replace=FALSE))
  # Generate permuted samples
  
  ARIperm=unlist(lapply(seq_len(perm),
                        function(i) mclust::adjustedRandIndex(Aperm[[i]], Bperm[[i]])))
  # Compute adjusted Rand index for the permuted samples
  
  m=mean(ARIperm); v=var(ARIperm)
  # compute mean and variance of the permuted samples
  
  NARI=(ARIperm-m)/sqrt(v)
  # compute NARI (normalized ARI)
  
  NARIo=(ARIo-m)/sqrt(v)
  # compute observed NARI
  
  p_value=(length(which(NARI>NARIo))+1)/(perm+1)
  # Compute p value as proportion of permuted NARI larger than the observed
  print(p_value)

  return(ARIperm)
}

set.seed(42)
values <- adjRand_test(randindexdata$to, randindexdata$from)
# Convert the array to a data frame and calculate the frequency of each value
df <- as.data.frame(table(values))

# Rename the columns for clarity
colnames(df) <- c("Value", "Frequency")
# Convert Value to numeric
df$Value <- as.numeric(as.character(df$Value))
df$AboveThreshold <- ifelse(df$Value > 0.0163, "Above", "Below")

# Create the bar plot
bar_plot <- ggplot(df, aes(x = Value, y = Frequency, fill = AboveThreshold)) +
  geom_bar(stat = "identity") +   # Use 'identity' to plot actual values
  scale_fill_manual(values = c("Below" = "gray", "Above" = "red")) +
  labs(title = "", x = "Adjusted Rand Index", y = "Frequency") +
  theme_bw() +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
  geom_vline(xintercept = 0.0163, linetype = "dashed", color = "red", linewidth = 1) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 1) +
  geom_vline(xintercept = 0.0163, linetype = "dashed", color = "red", linewidth = 1) +
  geom_hline(yintercept = 0, linetype = "solid", color = "black", linewidth = 1) +
  annotate("text", x = 0.05, y = 25, label = "Observed Adjusted Rand index = 0.0163\nP-value = 0.1309", angle = 0, vjust = 1.5, color = "red", size = 5) +
  theme(
    legend.position = "none",
    axis.title.x = element_text(size = 20),  # Increase x-axis title size
    axis.title.y = element_text(size = 14),  # Increase y-axis title size
    axis.text.x = element_text(size = 14),   # Increase x-axis text size
    axis.text.y = element_text(size = 12),    # Increase y-axis text size
    panel.grid.major = element_blank(),      # Remove major gridlines
    panel.grid.minor = element_blank(),      # Remove minor gridlines
  )

ggsave('randindexpermutations.png', bar_plot, scale = 1, width = 10, height = 8)
bar_plot





















chord = read.csv('cytof_chord.csv', sep = '\t', row.names = 1)

chordDiagram(chord, big.gap = 20, annotationTrack = "grid", preAllocateTracks = 1,
             order = c("Unstable", "Stable C 5", "Stable C 4", "Stable C 3",
                       "Stable C 2", "Stable C 1", "Stable C 0",
                       "C 0", "C 1", "C 2", "C 3", "C 4", "C 5")
)
circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  sector.name = get.cell.meta.data("sector.index")
  circos.text(mean(xlim), ylim[1] + .1, sector.name, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
  circos.axis(h = "top", labels.cex = 0.5, major.tick.percentage = 0.2, sector.index = sector.name, track.index = 2)
}, bg.border = NA)



chord = read.csv('RNA_chord.csv', sep = '\t', row.names = 1)
circos.par(start.degree = 90)
chordDiagram(chord, big.gap = 30,
            order = c("Unstable", "Stable C 4", "Stable C 3", "Stable C 2",
                       "Stable C 1", "Stable C 0", "C 0",
                       "C 1", "C 2", "C 3", "C 4")
             )







chord = read.csv('CytofRNA_chord.csv', sep = '\t', row.names = 1)
circos.par(start.degree = 90)
chordDiagram(chord, big.gap = 30,
             order = c("Unstable", "Stable C 4", "Stable C 3", "Stable C 2",
                       "Stable C 1", "Stable C 0", "C 0",
                       "C 1", "C 2", "C 3", "C 4")
)







chord = read.csv('CD4+CD8_chord.csv', sep = '\t', row.names = 1)
circos.par(start.degree = 90)

chordDiagram(chord[,c(1,2,3)], big.gap = 30,
             
### set visibility after correcting for multiple testing
#             link.visible = chord$p.value < (0.05/4),
             link.visible = chord$p.value < 0.05,
#             scale = T
             order = c("CD4_C 3", "CD4_C 2", "CD4_C 1", "CD4_C 0",
                       "CD8_C 0", "CD8_C 1", "CD8_C 2", "CD8_C 3")
)




multi_chord = read.csv('multi_chord_cmr_mr.csv', sep = '\t', row.names = 1)
multi_chord <- multi_chord %>%
  mutate(Cytof.Micro.RNA = recode(Cytof.Micro.RNA, 
                     "CMR_C 0" = "Trimodal - C 1",
                     "CMR_C 1" = "Trimodal - C 2",
                     "CMR_C 2" = "Trimodal - C 3",
                     "CMR_C 3" = "Trimodal - C 4"))
multi_chord <- multi_chord %>%
  mutate(Micro.RNA = recode(Micro.RNA, 
                                  "C 0" = "M&R Upregulated\nmetabolic sets",
                                  "C 1" = "M&R CMV seropositive",
                                  "C 2" = "M&R Higher BMI (Prevotella)",
                                  "C 3" = "M&R \nRuminococcaceae &\nChristensenellaceae - driven"))



circos.par(start.degree = 90)
par(cex = 1.5, mar = c(0, 0, 0, 0))
chordDiagram(multi_chord, big.gap = 30, order = c("M&R \nRuminococcaceae &\nChristensenellaceae - driven",
                                                  "M&R Higher BMI (Prevotella)", "M&R CMV seropositive", "M&R Upregulated\nmetabolic sets",
                                                  "Trimodal - C 1", "Trimodal - C 2", "Trimodal - C 3", "Trimodal - C 4") 
)





multi_chord = read.csv('multi_chord_cmr_cm.csv', sep = '\t', row.names = 1)
chordDiagram(multi_chord, big.gap = 30, order = c("C 3", "C 2", "C 1", "C 0",
                                                  "CMR_C 0", "CMR_C 1", "CMR_C 2", "CMR_C 3") 
)


multi_chord = read.csv('multi_chord_cmr_cr.csv', sep = '\t', row.names = 1)
multi_chord <- multi_chord %>%
  mutate(Cytof.Micro.RNA = recode(Cytof.Micro.RNA, 
                     "CMR_C 0" = "Trimodal - C 1",
                     "CMR_C 1" = "Trimodal - C 2",
                     "CMR_C 2" = "Trimodal - C 3",
                     "CMR_C 3" = "Trimodal - C 4"))
multi_chord <- multi_chord %>%
  mutate(Cytof.RNA = recode(Cytof.RNA, 
                            "C 0" = "C&R 1",
                            "C 1" = "C&R 2",
                            "C 2" = "C&R 3",
                            "C 3" = "C&R 4",
                            "C 4" = "C&R 5",
                            ))

chordDiagram(multi_chord, big.gap = 30, order = c("C&R 5", "C&R 4", "C&R 3", "C&R 2", "C&R 1",
                                                  "Trimodal - C 1", "Trimodal - C 2", "Trimodal - C 3", "Trimodal - C 4")
)





############# Bubble plots
  
# Define the matrix
tcrmodules_coefficients = as.matrix(read.csv('TCR_bubbles_coefficients.csv', sep = '\t', row.names = 1))
tcrmodules_pvalues = as.matrix(read.csv('TCR_bubbles_pvalues.csv', sep = '\t', row.names = 1))


# Convert the matrix to a dataframe
tcrmodules_coefficients <- as.data.frame(tcrmodules_coefficients)
# Set the row names
rownames(tcrmodules_coefficients) <- c("TCR C1", "TCR C2", "TCR C3", "TCR C4")
rownames(tcrmodules_pvalues) <- c("TCR C1", "TCR C2", "TCR C3", "TCR C4")

# Function to standardize columns
standardize <- function(x) {
  (x - mean(x)) / sd(x)
}

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

p
ggsave('bubble.png', p, scale = 1, width = 8, height = 5, dpi = 600)





### TODO change directory for the proper analysis
gsea_bacteria = read.csv('bubbleplot_C 0.csv', sep = '\t', row.names = 1)
gsea_bacteria$Cluster <- "C&R1"
aux = read.csv('bubbleplot_C 1.csv', sep = '\t', row.names = 1)
aux$Cluster <- "C&R2"
gsea_bacteria <- bind_rows(gsea_bacteria, aux)
aux = read.csv('bubbleplot_C 2.csv', sep = '\t', row.names = 1)
aux$Cluster <- "C&R3"
gsea_bacteria <- bind_rows(gsea_bacteria, aux)
aux = read.csv('bubbleplot_C 3.csv', sep = '\t', row.names = 1)
aux$Cluster <- "C&R4"
gsea_bacteria <- bind_rows(gsea_bacteria, aux)
aux = read.csv('bubbleplot_C 4.csv', sep = '\t', row.names = 1)
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
p
ggsave('CnR_GSEAbacteria.png', p, scale = 1, width = 8, height = 5, dpi = 600)




### TODO change directory for the proper analysis
gsea_bacteria = read.csv('bubbleplot_C 0.csv', sep = '\t', row.names = 1)
gsea_bacteria$Cluster <- "Prevotella-driven\nHigh BMI"
aux = read.csv('bubbleplot_C 1.csv', sep = '\t', row.names = 1)
aux$Cluster <- "Bacteroidaceae-driven"
gsea_bacteria <- bind_rows(gsea_bacteria, aux)
aux = read.csv('bubbleplot_C 2.csv', sep = '\t', row.names = 1)
aux$Cluster <- "Ruminococcaceae &\nChristensenellaceae driven\nLow BMI"
gsea_bacteria <- bind_rows(gsea_bacteria, aux)
aux = read.csv('bubbleplot_C 3.csv', sep = '\t', row.names = 1)
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
p
ggsave('uniM_GSEAbacteria.png', p, scale = 1, width = 8, height = 5, dpi = 600)



### TODO change directory for the proper analysis
gsea_bacteria = read.csv('bubbleplot_C 0.csv', sep = '\t', row.names = 1)
gsea_bacteria$Cluster <- "M&R Upregulated\nmetabolic sets"
aux = read.csv('bubbleplot_C 1.csv', sep = '\t', row.names = 1)
aux$Cluster <- "M&R CMV seropositive"
gsea_bacteria <- bind_rows(gsea_bacteria, aux)
aux = read.csv('bubbleplot_C 2.csv', sep = '\t', row.names = 1)
aux$Cluster <- "M&R Higher BMI\n(Prevotella)"
gsea_bacteria <- bind_rows(gsea_bacteria, aux)
aux = read.csv('bubbleplot_C 3.csv', sep = '\t', row.names = 1)
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
p
ggsave('MnR_GSEAbacteria.png', p, scale = 1, width = 8, height = 5, dpi = 600)
