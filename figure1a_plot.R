library(ggplot2)
library(dplyr)
library(phyloseq)
library(vegan)


#------------------------------------------------------------------------------#
#                                                                              #
#                          Figure 1a plot                                      #             
#                                                                              #   
#------------------------------------------------------------------------------#

mytheme = theme_bw() +
    theme(plot.margin = unit(c(1,1,1,1), "cm"), 
        plot.title = element_text(size = 10),
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.text = element_text(colour = "black", size = 10)) 


## bacterial constrained ordination
bac.vst.norm = readRDS("~/genotypes129/intermediate_data/16S/bac.abund(2489taxa).vst.norm.rmBatch.RDS")
bac.vst.norm = transformSampleCounts(bac.vst.norm,function(x) ifelse(x<0,0,x)) 
cap.bray.bac.vst = ordinate(bac.vst.norm, method = "CAP", distance = "bray", 
                            formula = ~Compartment+Compartment:Treatment)
# extract scores
score.sites = scores(cap.bray.bac.vst, display="sites",choice=c(1,2))
df.bac = cbind.data.frame(score.sites, as.data.frame(sample_data(bac.vst.norm)))
# extract percentage of PCs
screeplot(cap.bray.bac.vst, type="line")
ev = eigenvals(cap.bray.bac.vst) %>% summary() 
percentage.bac = paste( c("CAP1", "CAP2"), "(", paste( as.character(round(ev[2,c(1,2)]*100,2)), "%", ")", sep=""), sep = "" )
# plot
temp = ggplot(data= df.bac, aes(x = CAP1, y = CAP2, shape = Compartment)) +
  geom_point(aes(fill=Treatment), color='black', alpha=0.5, size=3, stroke=0.1) +
  stat_ellipse(linetype = 2, type = "norm", level = 0.9) +
  scale_fill_manual(values = c('#D3D3D3', '#F5C76D', '#90EE90','#ADD8E6')) +
  scale_shape_manual(values=c(21,22,24))+
  ggtitle(" constrained db-RDA plot of bacteria") + 
  mytheme +
  xlab(percentage.bac[1]) + ylab(percentage.bac[2]) +
  scale_x_continuous(limits=c(-0.7,0.7), breaks = seq(-1, 1, 0.5)) +
  scale_y_continuous(limits=c(-0.7,0.7), breaks = seq(-1, 1, 0.5)) +
  coord_fixed() 
# save pdf file
setwd("~/genotypes129/results/")
ggsave(temp, file = "./plots/16S/capscale-all-samples-rm-batch(2489taxa).pdf", 
       width = 20, height = 20, units = "cm")

## fungal constrained ordination
# read data and perform onstrained ordination using capscale
f.vst.norm = readRDS("~/genotypes129/intermediate_data/ITS/beta_div/fungi.abund.vst.norm.rmBatch.RDS")
f.vst.norm = transformSampleCounts(f.vst.norm,function(x) ifelse(x<0, 0, x)) 
cap.bray.f.vst = ordinate(f.vst.norm, method = "CAP", distance = "bray", 
                            formula = ~Compartment+Compartment:Treatment)
score.sites = scores(cap.bray.f.vst, display="sites",choice=c(1,2))
df = cbind.data.frame(score.sites, as.data.frame(sample_data(f.vst.norm)))
screeplot(cap.bray.f.vst, type="line")
ev = eigenvals(cap.bray.f.vst) %>% summary() 
percentage <- paste( c("CAP1", "CAP2"), "(", paste( as.character(round(ev[2,c(1,2)]*100,2)), "%", ")", sep=""), sep = "" )
# plot 
temp1 = ggplot(data= df, aes(x = CAP1, y = CAP2, shape = Compartment)) +
  geom_point(aes(fill=Treatment), color='black', alpha=0.5, size=3, stroke=0.1) +
  stat_ellipse(linetype = 2, type = "norm", level = 0.9) +
  scale_fill_manual(values = c('#D3D3D3', '#F5C76D', '#90EE90','#ADD8E6')) +
  scale_shape_manual(values=c(21,22,24))+
  ggtitle(" constrained db-RDA plot of fungi") + 
  mytheme +
  xlab(percentage[1]) + ylab(percentage[2]) +
  scale_x_continuous(limits=c(-2.5,2), breaks = seq(-2, 2, 1)) +
  scale_y_continuous(limits=c(-2.5, NA), breaks = seq(-2, 2, 1)) +
  coord_fixed() 
# save in pdf file
ggsave(temp1, file = "./plots/ITS/capscale-all-samples-rm-batch(248taxa)-fungi.pdf", 
       width = 20, height = 20, units = "cm")

library(patchwork)
p = temp + temp1
ggsave(p, file = "./plots/constrained-ordination-all-samples-bacteria-fungi.pdf", 
       width = 30, height = 15, units = "cm")

## statistic test
library(RVAideMemoire)

bray.dist <- phyloseq::distance(f.vst.norm, method = "bray") 
# PERMANOVA test 
# !!! when test for genotype, soil samples are excluded
permanova = adonis2(as.matrix(bray.dist)~Compartment+
                      Treatment, #+Plant.Name, 
                    data = as(sample_data(f.vst.norm), "data.frame"), 
                    permutations = 3999, parallel = 6, by="margin")


