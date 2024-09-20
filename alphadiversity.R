library(reshape2)
library(ggplot2)
library(dplyr)
library(ggpubr)


workingDir = "/Users/fabioaffaticati/Desktop/Work/activ_covid-and-omics"
setwd(workingDir) 

alpha = read.csv('results/alphadiversitydata.csv', sep = '\t', row.names = 1)

alpha$Cluster <- as.factor(alpha$Cluster)


p <- ggboxplot(alpha, x = "Cluster", y = "Shannon.entropy",
          color = "Cluster", palette = "jco",
                add = "jitter", shape = "Cluster")

p <- p + geom_pwc(
  method = "wilcox.test", label = "p.adj.format",
) +  stat_compare_means(method = "kruskal.test",label.y=9)   

ggadjust_pvalue(
  p, p.adjust.method = "bonferroni",
  label = "{p.adj.signif}", hide.ns = F
)                   




p <- ggboxplot(alpha, x = "Cluster", y = "Simpson.index",
               color = "Cluster", palette = "jco",
               add = "jitter", shape = "Cluster")

p <- p + geom_pwc(
  method = "wilcox.test", label = "p.adj.format",
) +  stat_compare_means(method = "kruskal.test",label.y=1.4)   

ggadjust_pvalue(
  p, p.adjust.method = "bonferroni",
  label = "{p.adj.signif}", hide.ns = F
)  




p <- ggboxplot(alpha, x = "Cluster", y = "Chao1.richness",
               color = "Cluster", palette = "jco",
               add = "jitter", shape = "Cluster")

p <- p + geom_pwc(
  method = "wilcox.test", label = "p.adj.format",
) +  stat_compare_means(method = "kruskal.test",label.y=330)   

ggadjust_pvalue(
  p, p.adjust.method = "bonferroni",
  label = "{p.adj.signif}", hide.ns = F
)

