library(rstudioapi)

# Get the path of the script's directory
script_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)

# Set working directory one level above the script's directory
setwd(dirname(script_dir))


# Define the directory path
res_path <- "results/R_scripts_plots/"
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

set.seed(42)



micro_data <- read.csv(paste0(data_path,"microbial_unimodal_adonis.csv"), row.names = 1)
metadata <- read.csv(paste0(data_path,"adonis_labels.csv"), row.names = 1)
Tax = read.csv('data/microbiome_analysis/taxa.csv', sep = ',', row.names = 1)



bc_micro <- vegdist(micro_data)
rownames(micro_data) <- metadata$Donor
metadata$labels <- factor(metadata$labels)

### It is important to note that the application of PERMANOVA assumes homogeneous
# group dispersions (variances). This can be tested with the PERMDISP2 method 
# (Anderson 2006) by using the same assay and distance method than in PERMANOVA.
anova(betadisper(bc_micro, metadata$labels))


permanova <- adonis2(bc_micro~labels, data = metadata, by = 'margin',
                     permutations = 1000, method='braycurtis')
permanova
str(permanova)
p_val <- permanova$`Pr(>F)`[1]


metadata %>% dplyr::count(labels)


PCoA <- vegan::betadisper(bc_micro, metadata$labels)
c <- cmdscale(bc_micro, eig = TRUE)


plot(PCoA, xlab=paste0('MDS1 Exp Var=', round(c$eig[1]/sum(abs(c$eig[c$eig])), 2)*100, '%'), xaxt='n',
           ylab=paste0('MDS2 Exp Var=', round(c$eig[2]/sum(abs(c$eig[c$eig])), 2)*100, '%'), yaxt='n')
boxplot(PCoA)

to_plot = as.data.frame(PCoA$vectors[,1], nm = c('PCoA1'))
to_plot$PCoA2 = PCoA$vectors[,2]
to_plot$Cluster = PCoA$group

# Define custom colors
custom_colors <- c("Prevotella-driven\nHigh BMI" = "#e31a1c",
                   "Bacteroidaceae driven" = "#33a02c",
                   "Ruminococcaceae &\nChristensenellaceae\ndriven Low BMI" = "#FFA500",
                   "Lachnospiraceae-driven" = "#1f78b4")


corr1 <- sapply(micro_data[], function(x) cor(to_plot$PCoA1, x))
corr2 <- sapply(micro_data[], function(x) cor(to_plot$PCoA2, x))

corr_matrix <- as.data.frame(cbind(corr1,corr2))
corr_matrix$x0 <- rep(0, nrow(corr_matrix))
corr_matrix$y0 <- rep(0, nrow(corr_matrix))

corr_matrix <- corr_matrix[order(abs(corr_matrix$corr1), abs(corr_matrix$corr2), decreasing = TRUE), ]
corr_matrix$taxa.taxon_id <- rownames(corr_matrix)


corr_matrix <- corr_matrix %>%
  left_join(Tax, by='taxa.taxon_id')



p <- ggplot(to_plot, aes(x=PCoA1, y=PCoA2, color=Cluster)) +
  geom_point(size = 4) +
  scale_color_manual(values = custom_colors) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
  labs(title = "") +
  theme_minimal() +
  theme(
    legend.title = element_text(size = 20),    
    legend.text = element_text(size = 16),      
    axis.title = element_text(size = 14)  
    )


for (i in seq(1:5)) {
  temp <- corr_matrix[corr_matrix$corr1>0,]
  temp <- temp[temp$corr2>0,]
  p <- p + geom_segment(
      x = temp$x0[i], y = temp$y0[i],
      xend = temp$corr1[i], yend = temp$corr2[i],
      arrow = arrow(length = unit(0.2, "cm")),
      color = "gray",  # Arrow color
      size = 0.5,      # Arrow size
      alpha = 0.5      # Arrow transparency
    ) + 
    annotate("text", label = temp$taxa.taxon_name[i],#temp$taxa.family[i],
             x = temp$corr1[i], y = temp$corr2[i],
             size = 6, colour = 'black')
}
for (i in seq(1:5)) {
  temp <- corr_matrix[corr_matrix$corr1<0,]
  temp <- temp[temp$corr2<0,]
  p <- p + geom_segment(
    x = temp$x0[i], y = temp$y0[i],
    xend = temp$corr1[i], yend = temp$corr2[i],
    arrow = arrow(length = unit(0.2, "cm")),
    color = "gray",  # Arrow color
    size = 0.5,      # Arrow size
    alpha = 0.5      # Arrow transparency
  ) + 
    annotate("text", label = temp$taxa.taxon_name[i],#temp$taxa.family[i],
             x = temp$corr1[i], y = temp$corr2[i],
             size = 6, colour = 'black')
}
for (i in seq(1:5)) {
  temp <- corr_matrix[corr_matrix$corr1>0,]
  temp <- temp[temp$corr2<0,]
  p <- p + geom_segment(
    x = temp$x0[i], y = temp$y0[i],
    xend = temp$corr1[i], yend = temp$corr2[i],
    arrow = arrow(length = unit(0.2, "cm")),
    color = "gray",  # Arrow color
    size = 0.5,      # Arrow size
    alpha = 0.5      # Arrow transparency
  ) + 
    annotate("text", label = temp$taxa.taxon_name[i],#temp$taxa.family[i],
             x = temp$corr1[i], y = temp$corr2[i],
             size = 6, colour = 'black')
}
for (i in seq(1:5)) {
  temp <- corr_matrix[corr_matrix$corr1<0,]
  temp <- temp[temp$corr2>0,]
  p <- p + geom_segment(
    x = temp$x0[i], y = temp$y0[i],
    xend = temp$corr1[i], yend = temp$corr2[i],
    arrow = arrow(length = unit(0.2, "cm")),
    color = "gray",  # Arrow color
    size = 0.5,      # Arrow size
    alpha = 0.5      # Arrow transparency
  ) + 
    annotate("text", label = temp$taxa.taxon_name[i],#temp$taxa.family[i],
             x = temp$corr1[i], y = temp$corr2[i],
             size = 6,colour = 'black')
}
ggsave(paste0(res_path,"PCoAbiplot.png"),p, scale = 1, width = 20, height = 12)






sc <- specc(PCoA$vectors, centers=3, iterations=1000)


### Plot the spectral clustering
scat <- ggplot(to_plot, aes(x=PCoA1, y=PCoA2, color=Cluster)) + 
  theme_bw() +
  geom_point() + 
  scale_color_manual(values = custom_colors) +
  stat_ellipse(level = 0.95) + 
  labs(x=paste0('PCoA1 Exp Var=', round(c$eig[1]/sum(c$eig[c$eig > 0]), 2)*100, '%'),
       y=paste0('PCoA2 Exp Var=', round(c$eig[2]/sum(c$eig[c$eig > 0]), 2)*100, '%'),
       xaxt='n', yaxt='n',) + 
  geom_rug(aes(color = Cluster)) + 
  theme(legend.position = "none",
        axis.ticks = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title = element_text(size = 16),  # Increase size of axis titles
        axis.text = element_text(size = 14),    # Increase size of axis text
        legend.text = element_text(size = 12))  # Increase legend text size)



xdensity <- ggplot(to_plot, aes(x=PCoA1, fill=Cluster)) + 
  geom_density(alpha=.5) + 
  theme_bw() +
  scale_fill_manual(values = custom_colors) +
  theme(legend.position = "none",
        axis.ticks = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title = element_text(size = 16),  # Increase size of axis titles
        axis.text = element_text(size = 14)    # Increase size of axis text
  )
ydensity <- ggplot(to_plot, aes(y=PCoA2, fill=Cluster)) + 
  geom_density(alpha=.5) + 
  theme_bw() +
  scale_fill_manual(values = custom_colors) +
  theme(legend.position = "none",
        axis.ticks = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title = element_text(size = 16),  # Increase size of axis titles
        axis.text = element_text(size = 14)    # Increase size of axis text
        )

blankPlot <- ggplot()+geom_blank(aes(1,1))+
  theme(plot.background = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(), 
        axis.text.y = element_blank(),
        axis.ticks = element_blank()
  )


p <- grid.arrange(xdensity, blankPlot, scat, ydensity, 
             ncol=2, nrow=2, widths=c(4, 1.4), heights=c(1.4, 4))

ggsave(paste0(res_path,"PCoA_density.png"),p, scale = 1, width = 20, height = 12)

