################################################################################
####                        figure 5a plot                                 #####
################################################################################

library(ggtree)
library(phyloseq)
library(ggplot2)
#library(ggtreeExtra)

# read data 
bac.asv = readRDS("~/genotypes129/intermediate_data/16S/ps-filtered-abund-815asv3140samples.RDS")
plot.asv = read.csv("67 ASVs_circus plot.csv", header = T, sep = "\t")
plot.asv = as.data.frame(plot.asv)
# extract samples and taxa
bac.asv = subset_samples(bac.asv, Compartment=='Root'&Treatment=='LN')
bac.asv1 = subset_taxa(bac.asv, taxa_names(bac.asv)%in%gsub('b', '', plot.asv$ASV))
# plot phylogenetic tree 
p0 <- ggtree(bac.asv1, layout="circular", branch.length='none', size=1, aes(color=Family)) +
  geom_tiplab(size=1, show.legend = FALSE, hjust = -0.5)
ggsave(p0, file = "./67asv-tree-draft1.pdf", width = 20, height = 20, units = "cm")
# first column must have same name with tree tip
plot.asv$ASV = gsub('b', '', plot.asv$ASV)
rownames(plot.asv) = plot.asv$ASV
# add tippoint with size proportion to average relative abundance
p = p0 %<+% plot.asv +
  geom_tippoint(aes(size=ave_RA*10000), alpha=0.9, color="grey") +
  scale_size(range = c(1,3))
# add heatmap of heritability for family, genus and ASV levels
p1 = gheatmap(p, plot.asv[, c(12,13,6)], offset = 2.5, color=NULL, 
         colnames_position="top", low = "white",
         high = "darkgreen", width = 0.15,
         colnames_angle=90, colnames_offset_y = 2, 
         hjust=0, font.size=1) +
  scale_fill_gradient(low = "white", high = "darkgreen", limits=c(0,0.6), na.value = "white")
# add barplot of sum of variance explained by SNPs for each ASV
p2 = p1 + geom_fruit_list(geom_fruit(geom=geom_bar,
                                    mapping=aes(y=ASV, x=Sum_R2),
                                    color="brown",
                                    offset =0.43, 
                                    orientation="y", 
                                    stat="identity",
                                    width=0.15,
                                    pwidth = 0.15,
                                    grid.params = list(vline=T, linetype=2),
                                    axis.params = list(axis="x", nbreak=4, vjust=1))) 
# save file
ggsave(p2, file = "./67asv-tree-draft7.pdf", width = 20, height = 20, units = "cm")
# negative prediction ability means not predictable
plot.asv[which(plot.asv[, 10] < 0, arr.ind = T),10] = NA
# plot heatmap of prediction ability
p3 = gheatmap(p, plot.asv[, c(8:10)], offset = 10, color=NULL, 
              colnames_position="top", low = "white",
              high = "blue", width = 0.15,
              colnames_angle=90, colnames_offset_y = 2, 
              hjust=0, font.size=1, legend_title = "Prediction ability") +
  scale_fill_gradient(low = "white", high = "blue", limits=c(0,0.5), na.value = "white")
# save file
ggsave(p3, file = "./67asv-tree-draft6.pdf", width = 20, height = 20, units = "cm")




