library(ggplot2)
library(RColorBrewer)
library(tidyverse)
library(DESeq2)
library(vegan)
library(phyloseq)
library(dplyr)
library(reshape2)
library(biomformat)
library(qiime2R)
library(FSA)
library(lme4)
library(rcompanion)
library(viridis)

mytheme = theme_bw() +
  theme(plot.margin = unit(c(1,1,1,1), "cm"), plot.title = element_text(size = 10),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        axis.text = element_text(colour = "black", size = 8)) 


setwd("~/genotypes129")


#-----------------------------------------------------------------------#
#                 import to phyloseq                                    #
#-----------------------------------------------------------------------#

# read asv biom table 
biom_table = read_biom("./data/ITS/feature-table.biom")
rm(biom_table)
feature.table = as.matrix(biom_data(biom_table))
dim(feature.table)                                             
# convert feature ID to ASV 
featureID2asv = cbind("FeatureID"=rownames(feature.table), "ASV"=paste("fASV", seq(1, nrow(feature.table)), sep = ""))
write.table(featureID2asv, "./intermediate_data/ITS/pre-analysis/fungi_feature_ID_to_ASV_file.txt", row.names = F, sep = "\t", col.names = T)
asv.table = feature.table
rownames(asv.table) = featureID2asv[,2]

# read taxonomy 
taxonomy = read.table("./data/ITS/taxonomy.tsv", sep = "\t", header = T)
taxa.table = taxonomy %>%
  separate(Taxon, c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = ";") %>%
  mutate(Kingdom = gsub("k__", "", Kingdom),
         Phylum = gsub("p__", "", Phylum),
         Class = gsub("c__", "", Class),
         Order = gsub("o__", "", Order),
         Family = gsub("f__", "", Family),
         Genus = gsub("g__", "", Genus),
         Species = gsub("s__", "", Species))
table(rownames(feature.table) == taxa.table$Feature.ID)
rownames(taxa.table) = featureID2asv[, 2]
taxa.table = taxa.table[, -1]
taxa.table = as.matrix(taxa.table)

# read metadata file
md = read.csv("./data/16S/samples_metadata_new.csv", header = T)

# create phyloseq object
OTU.tab = otu_table(asv.table, taxa_are_rows = TRUE)
TAX.tab = tax_table(taxa.table)
md = md[match(colnames(OTU.tab), md$Sample.ID), ]
rownames(md) = md$Sample.ID
sample.tab = sample_data(md)
table(sample.tab$Sample.ID == colnames(OTU.tab))
ps = phyloseq(OTU.tab, TAX.tab, sample.tab)
# read phylogenetic tree file
rooted.tree = qza_to_phyloseq(tree = "./data/ITS/rooted-tree.qza")
table(rooted.tree$tip.label %in% rownames(feature.table))
# change tip label to ASV number
new_tiplabel = featureID2asv[match(rooted.tree$tip.label, featureID2asv[,1]), 1]
table(rooted.tree$tip.label == new_tiplabel)
new_tiplabel = featureID2asv[match(rooted.tree$tip.label, featureID2asv[,1]), 2]
rooted.tree$tip.label = new_tiplabel
ps.tree = merge_phyloseq(ps, phy_tree(rooted.tree))
ps.tree                                             
saveRDS(ps.tree, "./intermediate_data/ITS/pre-analysis/raw_fungi_asv_table_phyloseq.RDS")


table(tax_table(ps.tree)[, "Phylum"])  # remove unidentified asv at phylum level
ps.tree.prefiltered = subset_taxa(ps.tree, Phylum!=' unidentified')

# function for filter taxa those express at least K reads in M*3 samples
expressKinM = function(K, M, otuTab){
  keep.taxa = c()
  for (i in 1:(ncol(otuTab)-1)) {
    print(i)
    if(sum(table(otuTab[otuTab[, i]>=K, ]$rep)>=3) >= M){
      keep.taxa = c(keep.taxa, colnames(otuTab)[i])
    }
  }
  return(keep.taxa)
}
# filter lowly-expressed taxa
tmp.otu = as.data.frame(cbind(t(otu_table(ps.tree.prefiltered)), "rep"=as.factor(rep(1:1056, each=3))))
keep.tax = expressKinM(10, 2, tmp.otu)  
ps.tree.filtered.abund = subset_taxa(ps.tree.prefiltered, taxa_names(ps.tree.prefiltered)%in%keep.tax) 
range(sample_sums(ps.tree.filtered.abund))

#--------------------------------------------------------------------#
#                  alpha diversity                                   #
#--------------------------------------------------------------------#

setwd("~/genotypes129/results/ITS/alpha_div")
# plot rarefaction curve
rarecurve(t(otu_table(ps.tree.filtered.abund)), step=100, label = F, sample = c(1000,8000))
# rarefaction 
ps.rare = rarefy_even_depth(ps.tree.filtered.abund, 1000, rngseed = 1625458, replace = FALSE)
# generate a data.frame with alpha diversity measures
adiv = data.frame(
  "Observed" = phyloseq::estimate_richness(ps.rare, measures = "Observed"),
  "Shannon" = phyloseq::estimate_richness(ps.rare, measures = "Shannon"))
adiv = cbind(adiv, sample_data(ps.rare))
write.table(adiv, file = "~/genotypes129/results/ITS/alpha_div/rarefied1000reads/alpha_diversity_table_1000.txt", sep = "\t", quote = F)

# box plot of alpha diversity
alphaDivPlot = function(alphaDivDF, plotFator1, plotFator2, levelsOrder, divIndex){
  print(divIndex)
  alphaDivDF$CompTre = paste(alphaDivDF[,plotFator1], alphaDivDF[,plotFator2], sep = "_")
  #f <- reformulate(divIndex, response="CompTre")
  f = as.formula(paste(divIndex, "CompTre", sep="~"))
  print(f)
  # Kruskal-Wallis test 
  krus.res = kruskal.test(f, data = alphaDivDF) 
  dunn.res = dunnTest(f, data = alphaDivDF, method = "bh")
  dunn.res$res$Compartment = gsub("_.*$", "", dunn.res$res$Comparison)
  #comOrd = c("Soil", "Rhizosphere", "Root")
  dunn.res$res = dunn.res$res[order(match(dunn.res$res$Compartment, levelsOrder)), ]
  write.csv(dunn.res$res, paste(divIndex, "dunnTest_results.csv", sep = "_"),  quote = F, row.names = F)
  # convert adjusted p value to letters
  label.dunn = cldList(P.adj ~ Comparison, data = dunn.res$res, threshold = 0.05)[, 1:2]
  label.dunn[, plotFator1] = gsub("_.*$", "", label.dunn$Group)
  label.dunn[, plotFator2] = sub("^[^_]*_", "", label.dunn$Group)
  label.dunn[, plotFator1] = factor(label.dunn[, plotFator1], levels = levelsOrder)
  print(label.dunn)
  alphaDivDF[, plotFator1] = factor(alphaDivDF[, plotFator1], levels = levelsOrder)
  #print(head(alphaDivDF))
  # box plot
  temp = ggplot(data= alphaDivDF, aes_string(x = plotFator1, y = divIndex, color=plotFator2)) +
    geom_boxplot(outlier.colour = NA) +
    geom_point(position=position_jitterdodge(jitter.width = 0.1, dodge.width = 0.75), 
               aes_string(color=plotFator2), alpha=0.6, size=1) +
    scale_color_viridis(discrete=TRUE) +
    geom_text(data = label.dunn, aes(y=rep(max(alphaDivDF[, divIndex])*1.1, nrow(label.dunn)), 
                                     label=Letter), position = position_dodge(width = 1)) +
    ylim(NA, max(alphaDivDF[, divIndex])*1.2) +
    labs(
      x = plotFator1,
      y = divIndex,
      title = paste(paste("Kruskal-Wallis chi-squared", signif(krus.res$statistic, 5), sep = "="), 
                    paste("df", krus.res$parameter, sep = "="),
                    paste("p value", signif(krus.res$p.value,3), sep = "="), sep = " ")
    ) + mytheme
  
  ggsave(temp, file = paste(divIndex, ".pdf", sep = ""), width = 16, height = 10, units = "cm")
}

alphaDivPlot(adiv, "Compartment", "Treatment", c("Soil", "Rhizosphere", "Root"), "Shannon")


#---------------------------------------------------------------------#
#                     beta diversity                                  #
#---------------------------------------------------------------------#

# DESeq2 Normalization 
bac.deseq = phyloseq_to_deseq2(ps.tree.filtered.abund, ~Block+Compartment)
bac.deseq = estimateSizeFactors(bac.deseq, type="poscount")
bac.vst = varianceStabilizingTransformation(bac.deseq, blind=F)   
# REMOVE batch effects
bac.vst$batch =as.factor(paste(bac.vst$Harvest.batch, bac.vst$Treatment.batch, sep = "_"))
bac.vst$plotBatch = as.factor(paste(bac.vst$Block, bac.vst$Plot, sep = "_"))
mm <- model.matrix(~Compartment, colData(bac.vst))
mat <- limma::removeBatchEffect(assay(bac.vst), batch=bac.vst$plotBatch, 
                                batch2 = bac.vst$batch, design=mm)
assay(bac.vst) <- mat
bac.vst.norm = ps.tree.filtered.abund
otu_table(bac.vst.norm) = otu_table(assay(bac.vst), taxa_are_rows = T)
# save normalized ASV file after removing batch effects
saveRDS(bac.vst.norm, "./ps.abund.vst.norm.rmBatch.RDS")


