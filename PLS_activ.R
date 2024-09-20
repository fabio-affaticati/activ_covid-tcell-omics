library(data.table)
library(tidyr)
library(mixOmics)
library(reticulate)

workingDir = "/Users/fabioaffaticati/Desktop/multiomics_tests/"
setwd(workingDir) 

set.seed(42)

metadata <- read.csv("metadata.csv", row.names=1, sep = '\t')
rna <- read.csv("RNAseq.csv", row.names=1, sep = '\t')
tcr <- read.csv("TCRseq.csv", row.names=1, sep = '\t')
microbiome <- read.csv("Microbiome.csv", row.names=1, sep = '\t')
cytof <- read.csv("Cytof.csv", row.names=1, sep = '\t')
cytof_meta <- read.csv("Cytof_meta.csv", row.names=1, sep = '\t')



rna_long <- pivot_longer(rna,
             cols = names(rna)[-which(names(rna) == 'Donor')],
             names_to = "feature",
             values_to = "value")

tcr_long <- pivot_longer(tcr,
                         cols = names(tcr)[-which(names(tcr) == 'Donor')],
                         names_to = "feature",
                         values_to = "value")

microbiome_long <- pivot_longer(microbiome,
                         cols = names(microbiome)[-which(names(microbiome) == 'Donor')],
                         names_to = "feature",
                         values_to = "value")

cytof_long <- pivot_longer(cytof,
                         cols = names(cytof)[-which(names(cytof) == 'Donor')],
                         names_to = "feature",
                         values_to = "value")
rna_long$view <- "RNA"
tcr_long$view <- "TCR"
microbiome_long$view <- "Micro"
cytof_long$view <- "Cytof"


# Extract common IDs using intersect
common_ids <- intersect(rna$Donor, cytof$Donor)
common_ids <- intersect(common_ids, metadata$Subjectnr)

# Keep rows present in both data frames
rna <- rna[rna$Donor %in% common_ids, ]
cytof <- cytof[cytof$Donor %in% common_ids, ]
cytof_meta <- cytof_meta[cytof_meta$Donor %in% common_ids, ]

metadata <- metadata[metadata$Subjectnr %in% common_ids, ]


pca.gene <- pca(rna[, -which(names(rna) == 'Donor')], ncomp = 10, center = TRUE, scale = TRUE)
pca.cytof <- pca(cytof_meta[, -which(names(cytof_meta) == 'Donor')], ncomp = 10, center = TRUE, scale = TRUE)

plot(pca.gene)
plot(pca.cytof)

plotIndiv(pca.gene, comp = c(1, 2), 
          group = metadata$Group, 
          ind.names = metadata$Subjectnr, 
          legend = TRUE, title = 'BTMs, PCA comp 1 - 2')

plotIndiv(pca.cytof, comp = c(1, 2), 
          group = metadata$Group, 
          ind.names = metadata$Subjectnr, 
          legend = TRUE, title = 'Cytof, PCA comp 1 - 2')



# set range of test values for number of variables to use from X dataframe
list.keepX <- c(seq(20, 50, 5))
# set range of test values for number of variables to use from Y dataframe
#list.keepY <- c(seq(10, 40, 5))
list.keepY <- c(seq(5, 10, 1))


tuned.spls <- tune.spls(rna[, -which(names(rna) == 'Donor')], cytof_meta[, -which(names(cytof_meta) == 'Donor')],
                             ncomp = 2,
                             test.keepX = list.keepX,
                             test.keepY = list.keepY,
                             nrepeat = 1, folds = 10, # use 10 folds
                             mode = 'canonical', measure = 'cor') 
plot(tuned.spls)         # use the correlation measure for tuning
# extract optimal number of variables for X dataframe


optimal.keepX <- tuned.spls$choice.keepX 

# extract optimal number of variables for Y datafram
optimal.keepY <- tuned.spls$choice.keepY

optimal.ncomp <-  length(optimal.keepX) # extract optimal number of components



pls.result <- spls(rna[, -which(names(rna) == 'Donor')], cytof_meta[, -which(names(cytof_meta) == 'Donor')], 
                   mode = 'canonical', ncomp = 5) # run the method


perf.spls <- perf(pls.result, validation = 'Mfold',
                        folds = 10, nrepeat = 5) 

plot(perf.spls , criterion = 'Q2.total')


# use all tuned values from above
pls.result <- spls(rna[, -which(names(rna) == 'Donor')], cytof_meta[, -which(names(cytof_meta) == 'Donor')], ncomp = optimal.ncomp, 
                         #keepX = optimal.keepX,
                         keepX = 50,
                         #keepY = optimal.keepY,
                         keepY = 22,
                         mode = "canonical") # explanitory approach being used, 
# hence use regression mode


plotIndiv(pls.result, 
          group = metadata$Group, 
          ind.names = metadata$Subjectnr,
          legend = TRUE)   # plot the samples
plotVar(pls.result, cex = c(4,4), )     # plot the variables

plotIndiv(pls.result, ind.names = FALSE, 
          rep.space = "X-variate", # plot in X-variate subspace
          group = metadata$Group,  # colour by time group
          pch = as.factor(metadata$Group), 
          col.per.group = color.mixo(1:3), 
          legend = TRUE, legend.title = 'Group', legend.title.pch = 'Group')

plotIndiv(pls.result, ind.names = FALSE,
          rep.space = "Y-variate", # plot in Y-variate subspace
          group = metadata$Group,  # colour by time group
          pch = as.factor(metadata$Group), 
          col.per.group = color.mixo(1:3), 
          legend = TRUE, legend.title = 'Group', legend.title.pch = 'Group')


plotIndiv(pls.result, ind.names = FALSE, 
          rep.space = "XY-variate", # plot in averaged subspace
          group = metadata$Group,  # colour by time group
          pch = as.factor(metadata$Group), 
          col.per.group = color.mixo(1:3), 
          legend = TRUE, legend.title = 'Group', legend.title.pch = 'Group')



plotArrow(pls.result, ind.names = FALSE, 
          rep.space = "XY-variate", # plot in averaged subspace
          group = metadata$Group,  # colour by time group
          col.per.group = color.mixo(1:3), 
          legend = TRUE, legend.title = 'Group', legend.title.pch = 'Group')



# form new perf() object which utilises the final model
perf.spls <- perf(pls.result, 
                        folds = 5, nrepeat = 10, # use repeated cross-validation
                        validation = "Mfold", 
                        dist = "max.dist",  # use max.dist measure
                        progressBar = FALSE)

# plot the stability of each feature for the first two components, 
# 'h' type refers to histogram
par(mfrow=c(1,2)) 
plot(perf.spls$features$stability.X[[1]], type = 'h',
     ylab = 'Stability',
     xlab = 'Features',
     main = '(a) Comp 1', las =2,
     xlim = c(0, 150))
plot(perf.spls$features$stability.X$comp2, type = 'h',
     ylab = 'Stability',
     xlab = 'Features',
     main = '(b) Comp 2', las =2,
     xlim = c(0, 300))



color.edge <- color.GreenRed(50)  # set the colours of the connecting lines

# X11() # To open a new window for Rstudio
net <- network(pls.result, comp = 1:2,
        cutoff = 0.5, # only show connections with a correlation above 0.7
        shape.node = c("rectangle", "circle"),
        #shape.node = c("none", "none"),
        color.node = c("cyan", "pink"),
        color.edge = color.edge,
        #interactive = T,
        #save = 'png', # save as a png to the current working directory
        name.save = 'sPLS Network Plot')
library(igraph)
write_graph(net$gR, file = "network.gml", format = "gml")

cim(pls.result, comp = 1:1, xlab = "Cytof", ylab = "BTMs", keysize=c(1,1),keysize.label = c(1), margins = c(20,15))




pls.result <- pls(rna[, -which(names(rna) == 'Donor')], cytof_meta[, -which(names(cytof_meta) == 'Donor')], ncomp = 2, 
                   mode = "canonical") # explanitory approach being used,

cim(pls.result, comp = 1:2, xlab = "Cytof", ylab = "BTMs", keysize=c(.5,.5), keysize.label = c(.5), margins = c(20,10),
    scale = TRUE)

png("annotation_highres.png", width = 8000, height = 12000, res = 300)

# Create the CIM again, but this time it will be saved as a high-resolution image
cim(pls.result, comp = 1:2, xlab = "Cytof", ylab = "BTMs", keysize=c(.5,.5),keysize.label = c(.5), margins = c(20,10),
    scale = TRUE)

# Close the graphics device
dev.off()


