library(WGCNA)
library(tidyverse)
library(dplyr)


options(stringsAsFactors = FALSE)
enableWGCNAThreads()

setwd('~/genotypes129/results/16S/WGCNA')

#===============================================================================
#
#  Read the ASV table and plot the sample tree
#
#===============================================================================

bac.ps = readRDS("~/genotypes129/intermediate_data/16S/bac.abund.vst.norm.rmBatch.RDS")
bac.rs = subset_samples(bac.ps, Compartment=='Rhizosphere') 
#bac.rt = subset_samples(bac.ps, Compartment=='Root') 
datExpr.bac = t(as.data.frame(otu_table(bac.rs)))
# build sample clusters
sampleTree = hclust(dist(datExpr.bac), method = "average")
# plot sample tree
pdf(file = "./Rhizo/1-n-sampleClustering.pdf", width = 80, height = 9)
par(cex = 0.8)
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="",
     cex.lab = 1.5,cex.axis = 1.5, cex.main = 2)
dev.off()

#===============================================================================
#
#  Choose soft threshold parameter
#
#===============================================================================

# Choose a set of soft threshold parameters
powers = c(c(1:10), seq(from = 12, to=30, by=2))
sft = pickSoftThreshold(datExpr.bac, powerVector = powers, verbose = 5, networkType = 'signed', 
                        corOptions = list(method='spearman')) 
# Scale-free topology fit index as a function of the soft-thresholding power
pdf(file = "./Rhizo/2-n-sft.pdf", width = 9, height = 5)
par(mfrow = c(1,2))
cex1 = 0.9
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
	xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
	main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
labels=powers,cex=cex1,col="red")
# this line corresponds to using an R^2 cut-off of h
abline(h=0.8,col="red") 
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
	xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
	main = paste("Mean connectivity")) 
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()

#===============================================================================
#
#  Turn data expression into topological overlap matrix
#
#===============================================================================

# Turn data expression into topological overlap matrix
power = sft$powerEstimate   
TOM = TOMsimilarityFromExpr(datExpr.bac, power = power, TOMType="signed")
#saveRDS(TOM, file = "./Rhizo/TOM_power6.RDS")
dissTOM = 1-TOM
# Plot gene tree
geneTree = hclust(as.dist(dissTOM), method = "average")
pdf(file = "./Rhizo/3-cluster-rhizo.pdf", width = 12, height = 9)
plot(geneTree, xlab="", sub="", main = "clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04)
dev.off()

#===============================================================================
#
#  Construct modules
#
#===============================================================================

# Module identification using dynamic tree cut
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 4, # The higher the value, the more smaller clusters will be produced. 
                            pamRespectsDendro = FALSE, minClusterSize = 5)
table(dynamicMods)
length(table(dynamicMods)) 
# Convert numeric labels into colors
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
# Plot the dendrogram and colors underneath
pdf(file = "4-module_tree.pdf", width = 8, height = 6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",dendroLabels = FALSE,
                    hang = 0.03,addGuide = TRUE, guideHang = 0.05,main = "Gene dendrogram and module colors")
dev.off()

#===============================================================================
#
#                            Merge modules
#
#===============================================================================

# Merge close modules
MEDissThres=0.25
merge = mergeCloseModules(datExpr.bac, dynamicColors, cutHeight = MEDissThres, verbose = 3) 
mergedColors.otu = merge$colors  
table(mergedColors.otu)
mergedMEs = merge$newMEs  
# Plot merged module tree
pdf(file = "./5-merged_Module_Tree.pdf", width = 20, height = 9)  
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors.otu), 
                    c("Dynamic Tree Cut", "Merged dynamic"), dendroLabels = FALSE, 
                    hang = 0.03, addGuide = TRUE, guideHang = 0.05)  
dev.off()

write.table(merge$newMEs,file="./newMEs.txt")

#===============================================================================
#
#  Plot the heatmap of module eigen-OTUs and samples
#
#===============================================================================

library("pheatmap")

# Heatmap of new module eigen-OTUs and samples
pdf(file="./Rhizo/newMEs-rhizo-signed.pdf",heigh=80,width=20)
row.names(merge$newMEs) = rownames(datExpr.bac)
pheatmap(merge$newMEs,cluster_col=T,cluster_row=T, breaks = c(seq(-0.5, 0.5, 0.01)),
         show_rownames=T,show_colnames=T,fontsize=10)
dev.off()

write.table(merge(merge$newMEs, sample_data(bac.rs), by="row.names"), 
            "./Rhizo/ME_merged_rhizo_signed.txt", 
            sep = "\t", quote = F, row.names = F)

#===============================================================================
#
#  Correlation between ASV modules and phenotypic traits
#
#===============================================================================

# Define numbers of ASVs and samples
nMicrobes = ncol(datExpr.bac)
nSamples = nrow(datExpr.bac)
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr.bac, mergedColors.otu)$eigengenes
MEs.bac = orderMEs(MEs0)   # equal to merge$newMEs  !!!
rownames(datExpr.bac) = gsub('S', '', rownames(datExpr.bac))
rownames(MEs.bac) = gsub('S', '', rownames(MEs.bac))

#### !!!! sample names should be consistent in eigen OTUs and traits !!!!
# read biomass and nitrogen data
biomass = read.csv('~/genotypes129/data/Summary_Nitrogen.csv', sep = ',')
biomass$ID = gsub('-', '_', biomass$ID)
biomass.scaled = apply(biomass[, 2:4], 2, scale)
rownames(biomass.scaled) = biomass$ID
biomass.scaled = as.data.frame(biomass.scaled)
# Calculate spearman correlation coefficients between module eigen-OTUs and traits
moduleTraitCor = cor(MEs.bac[,1:31], biomass.scaled[rownames(MEs.bac), ], use="p", method='spearman')
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
# check significant modules
rownames(moduleTraitCor)[moduleTraitPvalue[, 2] < 0.05 & moduleTraitCor[, 2]>0.1]

# read environmental data 
env.data = read.csv('~/genotypes129/data/97 lines environmental data.csv', sep = ',')
env.data.scaled = apply(env.data[, 2:ncol(env.data)], 2, scale)
rownames(env.data.scaled) = env.data$Genotype
# correlation between MEs AND environmental data
MEs.bac$Treatment = gsub('_.*', '', rownames(MEs.bac))
MEs.bac$Genotype = sapply(strsplit(rownames(MEs.bac), "_"), function(x) x[2])
# calculate mean value for each genotype under each soil treatment
MEs.bac.mean.ck = aggregate(.~Genotype, MEs.bac[MEs.bac$Treatment == 'CK', 1:32], mean)
# Calculate spearman correlation coefficients
moduleTraitCor = cor(MEs.bac.mean.ck[, 2:32], env.data.scaled[match(MEs.bac.mean.ck$Genotype, rownames(env.data.scaled)), c(3,41,119)], use = "p", method = 'spearman')
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nrow(env.data.scaled))

#===============================================================================
#
#  Plot heatmap of module-traits relationship
#
#===============================================================================

sizeGrWindow(10,6)
# display correlations and their p-values
textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 2), ")", sep = "")
dim(textMatrix) = dim(moduleTraitCor)
rownames(textMatrix) = colnames(MEs.bac.mean.ck)[2:32]
colnames(textMatrix) = colnames(env.data.scaled)[c(3,41,119)]

pdf("./Rhizo/rhizo_CK_bac_env_corr.pdf", width = 10, height = 20)
par(mar = c(15, 12, 5, 5))
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = colnames(moduleTraitCor),
               yLabels = rownames(moduleTraitCor),
               ySymbols = rownames(moduleTraitCor),
               xSymbols = colnames(moduleTraitCor),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 1,
               zlim = c(-1,1), 
               main = paste("bacteria module and traits relationships"))
dev.off()


#===============================================================================
#
#  Export of networks to external software
#
#===============================================================================

# Export the ASV list of new modules 
# Select modules
modules = gsub('ME', '', rownames(moduleTraitCor)[moduleTraitPvalue[, 2]<0.05 & moduleTraitCor[, 2]>0.1])
# Select module ASVs
microbes = colnames(datExpr.bac)
inModule = is.finite(match(mergedColors.otu, modules))
modOTUs = microbes[inModule]
modGenus = data.frame(tax_table(bac.rs))[modOTUs, 'Genus']
# Select the corresponding Topological Overlap
modTOM = TOM[inModule, inModule]
dimnames(modTOM) = list(modOTUs, modOTUs)
# Export the network into edge and node list files Cytoscape can read
cyt = exportNetworkToCytoscape(modTOM,
                               edgeFile = paste("CytoscapeInput-edges-", paste(modules, collapse="-"), ".txt", sep=""),
                               nodeFile = paste("CytoscapeInput-nodes-", paste(modules, collapse="-"), ".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0.0001,
                               nodeNames = modOTUs,
                               altNodeNames = modGenus,
                               nodeAttr = mergedColors.otu[inModule])

############ output all modules to node file ############

nodes = data.frame(cbind(microbes, mergedColors.otu, tax_table(bac.rs)[microbes, ]))
write.table(nodes, './Rhizo/bac_rhizosphere_modules_tab.txt', row.names = F, quote = F, sep = '\t')

#===============================================================================
#
#    plot genus distribution of each module
#
#===============================================================================

mods.nodes = read.table('./Rhizo/CytoscapeInput-nodes-royalblue-salmon-cyan-magenta-paleturquoise-white-lightgreen-pink-lightyellow-darkmagenta-saddlebrown-brown.txt', header = T)
# calculate the frequency of ASVs
nodes = mods.nodes %>% group_by(nodeAttr, altName) %>% summarise(n = n())
nodes.pl = nodes[nodes$n>1 & nodes$altName != 'uncultured' & !is.na(nodes$altName), ]
nodes.pl = rbind(nodes.pl, data.frame(nodeAttr = modules,
                                      altName = rep('other', 12),
                                      n = c(4, 6, 4, 13, 1, 3, 11, 11, 6, 3, 7, 37)))
nodes.pl$altName = factor(nodes.pl$altName, levels = c(names(table(nodes.pl[nodes.pl$altName != 'other', ]$altName)), "other"))
# Stacked barplot with multiple groups
library(randomcoloR)
palette = distinctColorPalette(22)
fig = ggplot(data=nodes.pl, aes(x=nodeAttr, y=n, fill=altName)) +
    geom_bar(stat="identity", position = position_stack(reverse = TRUE), width = 0.7) +
    scale_fill_manual(values = c(palette, 'grey')) +
    mytheme +
    theme(axis.text.x = element_text(angle =90, vjust=1, hjust=1)) +
    coord_flip() +
    ylab('number of ASVs') +
    xlab('Module')
ggsave(fig, file = './Rhizo/genus_number.pdf', height = 12, width = 25, units = 'cm')
# output selected modules ASV table
mods.tab = data.frame(otu_table(bac.rs)[mods.nodes$nodeName, ])
mods.tab = cbind(mods.tab, Module = mods.nodes$nodeAttr)
write.table(mods.tab, './Rhizo/rhizosphere_12_modules_ASV_table_normalized.txt', row.names = F, quote = F, sep = '\t')






