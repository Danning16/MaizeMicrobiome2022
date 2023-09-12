#-------------------------------------------------------------#
#                   Figure 1b plot                            #
#-------------------------------------------------------------#

library(multcomp)
library(FSA)
library(dplyr)

# read data
dat = read.csv("~/genotypes129/data/heritability_family_genus.csv", sep = "\t")
dat$compk = paste(dat$Compartment, dat$Kingdom, sep="_")
# remove rows without heritability
dat = dat[!is.na(dat$heritability), ]
# plot barplot for each taxon level
dat = dat[dat$Classification=="Genus", ]
dat$Treatment = as.factor(dat$Treatment)

###  !!! use non-parametric method that do not reqiure assumptions    !!!  ###

k.test = kruskal.test(H2~Treatment, dat[dat$Compartment=='Root.bac', ])
# create letter significance level 
label.dunn.all = c()
for (comp in c("Rhizosphere_Bacteria", "Root_Bacteria", "Rhizosphere_Fungi", "Root_Fungi")) {
    dunn.res = dunnTest(heritability~Treatment, data = dat[dat$compk==comp, ], method = "bh")
    label.dunn = cldList(P.adj ~ Comparison, data = dunn.res$res, threshold = 0.05)[, 1:2]
    label.dunn$compk = comp
    label.dunn$Treatment = label.dunn$Group
    label.dunn.all = rbind(label.dunn.all, label.dunn)
}
# set the order of display
label.dunn.all$compk = factor(label.dunn.all$compk, levels = c("Root_Bacteria","Rhizosphere_Bacteria", "Root_Fungi", "Rhizosphere_Fungi"))
label.dunn.all$Treatment = factor(label.dunn.all$Treatment, levels = c("CK","D", "LN", "LP", "Across"))
dat$compk= factor(dat$compk, levels = c("Root_Bacteria","Rhizosphere_Bacteria", "Root_Fungi", "Rhizosphere_Fungi"))
dat$Treatment = factor(dat$Treatment, levels = c("CK","D", "LN", "LP", "Across"))
# plot
mytheme = theme_bw() +
    theme(plot.margin = unit(c(1,1,1,1), "cm"), 
          plot.title = element_text(size = 10),
          panel.border = element_blank(), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          axis.line = element_line(colour = "black"),
          axis.text = element_text(colour = "black", size = 10)) 
temp = ggplot(data= dat, aes(x = compk, y = heritability, fill=Treatment)) +
    geom_boxplot(outlier.fill = 'white', color="black", outlier.shape = 21, outlier.stroke = 0.2,
                 position=position_dodge(width=0.75), width=0.5) +
    scale_fill_manual(values = c( '#D3D3D3', '#F5C76D', '#90EE90','#ADD8E6', 'white')) +
    geom_text(data = label.dunn.all, aes(y=rep(max(dat[, "heritability"])*1.1, nrow(label.dunn.all)), 
                                        label=Letter), position = position_dodge(width = 1)) +
    ylim(0, max(dat[, "heritability"])*1.2) +
    labs(
        x = "Compartment",
        y = "Heritability",
        title = "genus heritability"
    ) + mytheme

ggsave(temp, file = "~/genotypes129/results/plots/h2_genus.pdf", width = 16, height = 10, units = "cm")

#-----------------------------------------------------------------------------#
##                                    Pie plot                               ##
#-----------------------------------------------------------------------------#

# classify heritability into diff groups
dat1 = dat %>% group_by(compk, Treatment) %>% summarise(count = sum(heritability <=0.1)) %>% mutate(type = '0_0.1') 
dat1 = rbind(dat1, dat %>% group_by(compk, Treatment) %>% 
                 summarise(count = sum(heritability <=0.3 & heritability>0.1)) %>%
                 mutate(type = '0.1_0.3'))
dat1 = rbind(dat1, dat %>% group_by(compk, Treatment) %>% 
                 summarise(count = sum(heritability <=0.5 & heritability>0.3)) %>%
                 mutate(type = '0.3_0.5'))
dat1 = rbind(dat1, dat %>% group_by(compk, Treatment) %>% 
                 summarise(count = sum(heritability > 0.5)) %>%
                 mutate(type = '0.5'))
# set the order of x-axis labels
dat1$comb = paste(dat1$compk, dat1$Treatment)
dat2 = dat1 %>% group_by(comb) %>% mutate(num = sum(count))
dat2$comb = factor(dat2$comb, levels = c("Root_Bacteria CK", "Root_Bacteria D",
                                         "Root_Bacteria LN", "Root_Bacteria LP",
                                         "Root_Bacteria Across", "Rhizosphere_Bacteria CK",
                                         "Rhizosphere_Bacteria D", "Rhizosphere_Bacteria LN",
                                         "Rhizosphere_Bacteria LP", "Rhizosphere_Bacteria Across",
                                         "Root_Fungi CK", "Root_Fungi D", "Root_Fungi LN",
                                         "Root_Fungi LP", "Root_Fungi Across", 
                                         "Rhizosphere_Fungi CK", "Rhizosphere_Fungi D",
                                         "Rhizosphere_Fungi LN","Rhizosphere_Fungi LP",
                                         "Rhizosphere_Fungi Across"))

# Basic pie chart
g1 = ggplot(dat2, aes(x="", y=count, fill=type)) +
    geom_bar(position = "fill", stat="identity") +
    facet_grid(.~ comb) + 
    coord_polar("y", start=0) +
    theme_void() +
    scale_fill_manual(values = c('tan', 'tan1','tan3', 'tan4'))
ggsave(g1, file = "~/genotypes129/results/plots/h2_genus_pie.pdf", width = 16, height = 10, units = "cm")





