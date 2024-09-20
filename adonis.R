library(tidyverse)
library(vegan)
library(ape)
library(gridExtra)
library(kernlab)
library(Rtsne)

set.seed(42)

workingDir = "/Users/fabioaffaticati/Desktop/Work/activ_covid-and-omics/"
setwd(workingDir) 

micro_data <- read.csv("microbial_unimodal_adonis.csv",row.names = 1)
metadata <- read.csv("adonis_labels.csv",row.names = 1)

bc_micro <- vegdist(micro_data)

#colnames(micro_data) <- metadata$Donor
rownames(micro_data) <- metadata$Donor


metadata$labels <- factor(metadata$labels)
#micro_data <- micro_data %>% as.dist

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

to_plot <- to_plot %>%
  mutate(Cluster = recode(Cluster, 
                          "C 0" = "Prevotella-driven\nHigh BMI",
                          "C 1" = "Bacteroidaceae-driven",
                          "C 2" = "Ruminococcaceae &\nChristensenellaceae-driven\nLow BMI",
                          "C 3" = "Lachnospiraceae-driven"))
# Define custom colors
custom_colors <- c("Prevotella-driven\nHigh BMI" = "#e31a1c",
                   "Bacteroidaceae-driven" = "#33a02c",
                   "Ruminococcaceae &\nChristensenellaceae-driven\nLow BMI" = "#FFA500",
                   "Lachnospiraceae-driven" = "#1f78b4")


corr1 <- sapply(micro_data[], function(x) cor(to_plot$PCoA1, x))
corr2 <- sapply(micro_data[], function(x) cor(to_plot$PCoA2, x))

corr_matrix <- as.data.frame(cbind(corr1,corr2))
corr_matrix$x0 <- rep(0, nrow(corr_matrix))
corr_matrix$y0 <- rep(0, nrow(corr_matrix))

corr_matrix <- corr_matrix[order(abs(corr_matrix$corr1), abs(corr_matrix$corr2), decreasing = TRUE), ]
corr_matrix$taxa.taxon_id <- rownames(corr_matrix)

Tax = read.csv('data/microbiome_analysis/taxa.csv', sep = ',', row.names = 1)

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
    legend.title = element_text(size = 20),    # Increase legend title size
    legend.text = element_text(size = 16),      # Increase legend text size
    axis.title = element_text(size = 14)  
    )
p






#for (i in seq(1:nrow(corr_matrix))) {
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
ggsave(filename = "results/PCoAbiplot.png",p, scale = 1, width = 20, height = 12)
p





sc <- specc(PCoA$vectors, centers=3, iterations=1000)

#sc
#centers(sc)
#size(sc)
#withinss(sc)


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
ggsave(filename = "results/PCoAcircleplot.png",p, scale = 1, width = 20, height = 12)
p



names(metadata)[names(metadata) == "Donor"] <- 'sample_id'
metadata$Cluster <- plot.spectral$Cluster
### pass relative abundances
Tab = read.csv('results/relative_abundances_forplot.csv', sep = '\t', row.names = 1)


### read or pass the taxa information
Tax = read.csv('data/microbiome_analysis/taxa.csv', sep = ',', row.names = 1)



Tab$taxa.taxon_id <-rownames(Tab)


dat <- Tab %>%
  pivot_longer(-taxa.taxon_id, names_to = "sample_id", values_to = 'rel_ab')

dat <- dat %>%
  left_join(Tax, by='taxa.taxon_id')

dat$sample_id <-gsub(".", " ", dat$sample_id, fixed=TRUE)


dat <- dat %>%
  left_join(metadata, by='sample_id')

family_order <- c('Other',
                  'Bacteroidaceae', 
                  'Lachnospiraceae', 
                  'Prevotellaceae', 
                  'Ruminococcaceae',
                  #'Rikenellaceae',
                  'Christensenellaceae')

shouldBecomeOther <- !(dat$taxa.family %in% family_order)
dat$taxa.family[shouldBecomeOther]<- "Other"

dat <- dat %>%
  mutate(taxa.family = factor(taxa.family, levels = family_order))


dat %>%
  ggplot(aes(x = sample_id, y = rel_ab)) + 
  facet_grid(cols = vars(Cluster), scales = 'free_x', space = 'free_y') +
  geom_bar(aes(fill = taxa.family), stat = 'identity', position = 'fill', width = 1) +
  scale_fill_brewer(palette = 'Paired') +
  scale_y_continuous(name = 'Relative abundance',
                     labels = scales::percent) + 
  theme(axis.text.x = element_text(angle = 90, size = 6),
        axis.title.y = element_text(color = 'black'),
        strip.text = element_text(face = 'bold'),
        strip.background = element_blank())









##### RNA COVID
workingDir = "/Users/fabioaffaticati/Desktop/Work/activ_covid-and-omics/results/"
setwd(workingDir) 


countMat <- read.csv("normalized_counts_renamed.csv",row.names = 1, sep='\t')
meta.data <- read.csv("patients_raw_counts.csv", sep='\t' ,row.names = 1)

countMat = countMat[!duplicated(countMat$genename),]
rownames(countMat) <- countMat$genename
countMat <- subset(countMat, select = -genename)

rownames(meta.data) <- meta.data$Subjectnr
head(meta.data)


dist_RNA <- dist(countMat)


meta.data$labels <- factor(meta.data$labels)

anova(betadisper(dist_RNA, meta.data$labels))


permanova <- adonis2(dist_RNA~labels, data = meta.data, by = 'margin',
                     permutations = 100, method='euclidean')
permanova
str(permanova)
p_val <- permanova$`Pr(>F)`[1]


metadata %>% dplyr::count(labels)


PCoA <- vegan::betadisper(dist_RNA, meta.data$labels)
c <- cmdscale(dist_RNA, eig = TRUE)


plot(PCoA, xlab=paste0('MDS1 Exp Var=', round(c$eig[1]/sum(c$eig[c$eig > 0]), 2)*100, '%'), xaxt='n',
     ylab=paste0('MDS2 Exp Var=', round(c$eig[2]/sum(c$eig[c$eig > 0]), 2)*100, '%'), yaxt='n')
boxplot(PCoA)









pairwise_p <- numeric()

########## C 0 vs C 1
meta_C0_C1 <- metadata %>% filter(labels == "C 0"| labels == "C 1")
micro_data_C0_C1 <- micro_data[rownames(micro_data) %in% meta_C0_C1$Donor,]
permanova_C0_C1 <- adonis2(micro_data_C0_C1~labels, data = meta_C0_C1, permutations = 1000)
permanova_C0_C1
pairwise_p["permanova_C0_C1"] <- permanova_C0_C1$`Pr(>F)`[1]

########## C 0 vs C 2
meta_C0_C1 <- metadata %>% filter(labels == "C 0"| labels == "C 2")
micro_data_C0_C1 <- micro_data[rownames(micro_data) %in% meta_C0_C1$Donor,]
permanova_C0_C1 <- adonis2(micro_data_C0_C1~labels, data = meta_C0_C1, permutations = 1000)
permanova_C0_C1
pairwise_p["permanova_C0_C2"] <- permanova_C0_C1$`Pr(>F)`[1]

########## C 0 vs C 3
meta_C0_C1 <- metadata %>% filter(labels == "C 0"| labels == "C 3")
micro_data_C0_C1 <- micro_data[rownames(micro_data) %in% meta_C0_C1$Donor,]
permanova_C0_C1 <- adonis2(micro_data_C0_C1~labels, data = meta_C0_C1, permutations = 1000)
permanova_C0_C1
pairwise_p["permanova_C0_C3"] <- permanova_C0_C1$`Pr(>F)`[1]

########## C 1 vs C 3
meta_C0_C1 <- metadata %>% filter(labels == "C 1"| labels == "C 3")
micro_data_C0_C1 <- micro_data[rownames(micro_data) %in% meta_C0_C1$Donor,]
permanova_C0_C1 <- adonis2(micro_data_C0_C1~labels, data = meta_C0_C1, permutations = 1000)
permanova_C0_C1
pairwise_p["permanova_C1_C3"] <- permanova_C0_C1$`Pr(>F)`[1]

########## C 1 vs C 2
meta_C0_C1 <- metadata %>% filter(labels == "C 1"| labels == "C 2")
micro_data_C0_C1 <- micro_data[rownames(micro_data) %in% meta_C0_C1$Donor,]
permanova_C0_C1 <- adonis2(micro_data_C0_C1~labels, data = meta_C0_C1, permutations = 1000)
permanova_C0_C1
pairwise_p["permanova_C1_C2"] <- permanova_C0_C1$`Pr(>F)`[1]



p.adjust(pairwise_p, method = "BH")