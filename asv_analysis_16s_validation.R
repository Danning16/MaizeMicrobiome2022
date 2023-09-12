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

setwd("~/genotypes129")

#==============================================================================#
#
#                 validation experiment 16S data analysis                      #
#
#==============================================================================#

# import to phyloseq
# read asv table file in biom format
biom_table = read_biom("~/genotypes129/data/new_validation/data/clean/feature-table.biom")
feature.table = as.matrix(biom_data(biom_table))
dim(feature.table)    
# create bacteria feature ID to ASV file
featureID2asv = cbind("FeatureID"=rownames(feature.table), "ASV"=paste("ASV", seq(1, nrow(feature.table)), sep = ""))
write.table(featureID2asv, "~/genotypes129/data/new_validation/data/clean/feature_ID_to_ASV_file.txt", row.names = F, sep = "\t", col.names = T)
asv.table = feature.table
rownames(asv.table) = featureID2asv[,2]
# read taxonomy file
taxonomy = read.table("~/genotypes129/data/new_validation/data/clean/taxonomy.tsv", sep = "\t", header = T)
taxa.table = taxonomy %>%
    separate(Taxon, c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = "; ") %>%
    mutate(Kingdom = gsub("d__", "", Kingdom),
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
taxa.table[is.na(taxa.table)] = 'unidentified'
# read metadata file
md = read.csv("~/genotypes129/data/validation/metadata.csv", header = T, sep = ',')
md = md[order(md$repID, md$Treatment, md$Compartment), ]
asv.table = asv.table[, order(colnames(asv.table))]
colnames(asv.table) = paste(md$repID, md$Treatment, md$Compartment, sep = '_')
md$sampleID = colnames(asv.table)
# create phyloseq object
OTU.tab = otu_table(asv.table, taxa_are_rows = TRUE)
TAX.tab = tax_table(taxa.table)
rownames(md) = md$sampleID
sample.tab = sample_data(md)
ps = phyloseq(OTU.tab, TAX.tab, sample.tab)
saveRDS(ps, '~/genotypes129/data/new_validation/data/clean/asv_ps_raw.RDS')
# ps = readRDS('~/genotypes129/data/new_validation/data/clean/asv_ps_raw.RDS')
# keep only bacteria and remove Mitochondria and Chloroplast ASVs 
ps.prefiltered = subset_taxa(ps, Kingdom == 'Bacteria' & 
                                 Family != 'Mitochondria' & 
                                 Order != 'Chloroplast')
# check if every sample has >1000 reads
range(colSums(otu_table(ps.prefiltered)) )
range(rowSums(otu_table(ps.prefiltered)) )

sum(taxa_sums(ps.prefiltered)) / sum(taxa_sums(ps))    # =0.57

# ra > 0.05% at least in 5% samples
ps.RA = transform_sample_counts(ps.prefiltered, function(x) x/sum(x))
ps.ra.abund = filter_taxa(ps.RA, function(x) sum(x > 0.0005) >= 0.05*nsamples(ps.RA), TRUE)  # 1040 taxa
ps.abund = subset_taxa(ps.prefiltered, taxa_names(ps.prefiltered)%in%taxa_names(ps.ra.abund))
range(colSums(otu_table(ps.abund)) )  # all samples >1000 reads
sum(taxa_sums(ps.abund)) / sum(taxa_sums(ps.prefiltered))   # = 0.88
saveRDS(ps.abund, '~/genotypes129/results/validation/intermediate_data/clean/ps.abund.root.RDS')

#===============================================================================
#
#                alpha diversity
#
#===============================================================================
# choose rarefaction threshold
rarecurve(t(data.frame(otu_table(ps.abund))), step=100, label = F, sample = c(1000,10000))
# rarefaction at 8000 reads
ps.rare = rarefy_even_depth(ps.abund, 8000, rngseed = 123555, replace = FALSE)
# generate a data frame with alpha-diversity measures
adiv = data.frame(
    "Observed" = phyloseq::estimate_richness(ps.rare, measures = "Observed"),
    "Shannon" = phyloseq::estimate_richness(ps.rare, measures = "Shannon"))
adiv = cbind(adiv, sample_data(ps.rare))
write.table(adiv, file = "~/genotypes129/results/validation/intermediate_data/clean/alpha_diversity_table_8000.txt", sep = "\t", quote = F)

## BOX plot of alpha diversity
library(FSA)
library(rcompanion)
# do Dunn test 
adiv$Genotype = gsub('-', '_', adiv$Genotype)
adiv$CompTre = paste(adiv$Compartment, adiv$Treatment, sep = "_")
adiv$CompTre = as.factor(adiv$CompTre)
kruskal.test(Shannon~CompTre, data = adiv) 
dunn.res = dunnTest(Shannon~CompTre, data = adiv, method = "bh")
label.dunn = cldList(P.adj ~ Comparison, data = dunn.res$res, threshold = 0.05)[, 1:2]   # genotype does not affect alpha-div
label.dunn$Compartment = gsub("_.*$", "", label.dunn$Group)
label.dunn$Treatment = sub("^[^_]*_", "", label.dunn$Group)
# define my plot theme
mytheme = theme_bw() +
    theme(plot.margin = unit(c(1,1,1,1), "cm"), plot.title = element_text(size = 10),
          panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
          axis.text = element_text(colour = "black", size = 10)) 

adiv$Genotype = factor(adiv$Genotype, levels = c('B73', 'F_0598', 'F7', 'D_0170', 'Bulk_soil'))
pal = viridisLite::viridis(4)

fig5a = ggplot(data = adiv, aes(x = Compartment, y = Shannon, fill=Treatment)) +
    geom_boxplot(outlier.colour = NA) +
    geom_point(position=position_jitterdodge(jitter.width = 0.1, dodge.width = 0.75), 
               aes(fill=Treatment, color=Genotype), alpha=0.6, size=2, stroke=0.5) +
    scale_fill_manual(values = c("brown", "orange")) +
    scale_color_manual(values = c(pal, 'black')) + 
    geom_text(data = label.dunn, aes(y=rep(max(adiv[, "Shannon"])*1.1, nrow(label.dunn)), 
                                     label=Letter), position = position_dodge(width = 1)) +
    ylim(NA, max(adiv[, "Shannon"])*1.2) +
    mytheme

ggsave(fig5a, file = "~/genotypes129/results/validation/v4/clean/Shannon_rarefied8000.pdf", width = 15, height = 10, units = "cm")

#===============================================================================
#
#                     beta diversity
#
#===============================================================================

# DESeq2 Normalization 
bac.deseq = phyloseq_to_deseq2(ps.abund, ~Compartment+Treatment)
bac.deseq = estimateSizeFactors(bac.deseq, type="poscount")
bac.vst = varianceStabilizingTransformation(bac.deseq, blind=F)   
bac.vst.norm = ps.abund
otu_table(bac.vst.norm) = otu_table(assay(bac.vst), taxa_are_rows = T)
# Round negative values up to zeroes, to enable Bray-Curtis calculations
bac.vst.norm = transformSampleCounts(bac.vst.norm,function(x) ifelse(x<0,0,x)) 
# do PCoA analysis
PCoA.bray.vst = ordinate(bac.vst.norm, method = "PCoA", distance = "bray")
pcoa.df = cbind.data.frame(PCoA.bray.vst$vectors[, 1:5], sample_data(bac.vst.norm))
pct = PCoA.bray.vst$values$Relative_eig[1:2]
percentage = paste( c("PCoA1", "PCoA2"), " (", paste( as.character(round(pct*100,2)), "%", ")", sep=""), sep = "" )
pcoa.df$Genotype = factor(pcoa.df$Genotype, levels = c('B73', 'F-0598', 'F7', 'D-0170', 'Bulk-soil'))
# PCoA plot
fig5b = ggplot(data= pcoa.df, aes(x = Axis.1, y = Axis.2, fill = Treatment, shape = Compartment)) +
    geom_point(aes(color=Genotype), alpha=1, size=2, stroke=1) +
    scale_shape_manual(values = c(21:25)) +
    stat_ellipse(linetype = 2) +
    scale_fill_manual(values = c("brown", "orange")) +
    scale_color_manual(values = c(pal, 'black')) + 
    ggtitle("PCoA plot ") + 
    mytheme +
    xlab(percentage[1]) + ylab(percentage[2]) +
    coord_fixed() 

## PERMANOVA test ##
#test if treatment affect the beta diversity
bray.dist.part = phyloseq::distance(bac.vst.norm, method = "bray") 
permanova = adonis2(bray.dist.part~Compartment+Treatment*Genotype, 
                    data = as(sample_data(bac.vst.norm), "data.frame"), 
                    permutations = 3999, parallel = 18)

# test for genotype effects 
bac.vst.norm.rt = subset_samples(bac.vst.norm, Compartment == 'Root')
bray.dist.rt = phyloseq::distance(bac.vst.norm.rt, method = "bray") 
# define permutation blocks
perm = how(nperm = 3999)
md = sample_data(bac.vst.norm.rt)
setBlocks(perm) = with(md, Treatment)
permanova2 = adonis2(bray.dist.rt~Treatment*Genotype, 
                     data = as(sample_data(bac.vst.norm.rt), "data.frame"), 
                     permutations = perm, parallel = 18)

#===============================================================================
#
#  identify genus significantly different between WT and mutants         
#
#===============================================================================
library(rcompanion)
library(multcomp)
library(viridis)

# group ASVs to genus level
bac.ge = tax_glom(ps.abund, taxrank = "Genus")
bac.ge = subset_taxa(bac.ge, !Genus %in% c("uncultured", 'unidentified') )   # 235 taxa
# extract only root samples
bac.ge = subset_samples(bac.ge, Compartment == 'Root')
bac.ge.ra = transformSampleCounts(bac.ge, function(x){x/sum(x)})
# only keep genera with RA > 1% in at least 1 samples
bac.ge.abund = filter_taxa(bac.ge.ra, function(x) sum(x > 0.01) >= 2, TRUE)   # 34 taxa
bac.ge.otutab = data.frame(otu_table(bac.ge.abund))
rownames(bac.ge.otutab) = tax_table(bac.ge.abund)[, 'Genus']
# create data frame for box plot
bac.ge.df = cbind.data.frame(t(bac.ge.otutab), sample_data(bac.ge.abund))
pl.df = pivot_longer(bac.ge.df, cols = rownames(bac.ge.otutab), names_to = c('Genus'))
pl.df$Genotype = as.factor(pl.df$Genotype)
pl.df$Genotype = gsub('-', '_', pl.df$Genotype)

# do wilcoxon test for each genus between wild type and its mutant under each nitrogen level
wilcoxon.res = c()
for (treat in c('HN', 'LN')) {
    for (comp in c('Root')) {
        for (g in unique(pl.df$Genus)) {
            df = pl.df[pl.df$Genus == g & pl.df$Treatment == treat & pl.df$Compartment == comp, ]
            formula1 = as.formula(paste('value', 'Genotype', sep = '~'))
            # for the first pair
            if(sum(df[df$Genotype %in% c('B73', 'D_0170'), 'value']) != 0){
                w.res1 = wilcox.test(formula1, data = df[df$Genotype %in% c('B73', 'D_0170'), ]) 
                print(w.res1$p.value)
                if(w.res1$p.value < 0.05){
                    wilcoxon.res = rbind(wilcoxon.res, cbind(Compartment = comp, Treatment = treat, Genus = g, Genotype = 'B73', pval = w.res1$p.value))
                }
            }
            # for the other pair
            if(sum(df[df$Genotype %in% c('F7', 'F_0598'), 'value']) != 0){
                w.res2 = wilcox.test(formula1, data = df[df$Genotype %in% c('F7', 'F_0598'), ]) 
                print(w.res2$p.value)
                if(w.res2$p.value < 0.05){
                    wilcoxon.res = rbind(wilcoxon.res, cbind(Compartment = comp, Treatment = treat, Genus = g, Genotype = 'F7', pval = w.res2$p.value))
                }
            }
        }
    }
}
# no significant genus between WT and its mutant under high nitrogen condition
wilcoxon.res = data.frame(wilcoxon.res)
wilcoxon.res$pval = as.numeric(wilcoxon.res$pval)
wilcoxon.res$label = cut(wilcoxon.res$pval, breaks = c(-Inf, 0.001, 0.01, 0.05, Inf), 
                         labels = c("***", "**", "*", "n.s."), right = FALSE)
write.table(wilcoxon.res, '~/genotypes129/results/validation/intermediate_data/clean/Genus_RA_wilcoxonTest_results.txt', row.names = F, sep = '\t', quote = F)

## BOX plot Relative abundance of genus significantly diff between WT and mutant  ##
mytheme = theme_bw() +
    theme(plot.margin = unit(c(1,1,1,1), "cm"), plot.title = element_text(size = 10),
          panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
          axis.text = element_text(colour = "black", size = 8)) 

fig.df = pl.df[pl.df$Genus %in% unique(wilcoxon.res$Genus) & pl.df$Compartment == 'Root' & pl.df$Treatment == 'LN', ]
fig.df$Genotype = factor(fig.df$Genotype, levels = c('B73', 'D_0170', 'F7', 'F_0598'))
fig.df$Genus = factor(fig.df$Genus, levels = unique(wilcoxon.res$Genus))

fig5c = ggplot(data= fig.df, aes(x = Genus, y = value, color=Genotype)) +
    geom_boxplot(outlier.colour = NA) +
    geom_point(position=position_jitterdodge(jitter.width = 0.1, dodge.width = 0.75), 
               aes(color=Genotype), alpha=0.6, size=2, stroke=0) +
    scale_color_viridis(discrete=TRUE) +
    geom_text(data = wilcoxon.res,
              aes(x = Genus, y = rep(0.3, 3), label = label, color = Genotype), inherit.aes = F,
              position = position_dodge(width = 0.75)) +
    # ylim(NA, max(adiv[, "Shannon"])*1.2) +
    labs(
        y = "Relative abundance",
        title = "Wilcoxon rank-sum test"
    ) +
    mytheme +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))


ggsave(fig5c, file = "~/genotypes129/results/validation/v4/clean/root_LN_genus_RA_boxplot_new.pdf", 
       width = 15, height = 10, units = "cm")

fig5 = ggarrange(fig5a, fig5b, fig5c, nrow = 1, align = 'h')
ggsave(fig5, file = "~/genotypes129/results/validation/v4/clean/fig5-1.pdf", 
       width = 35, height = 10, units = "cm")


