library(ggplot2)
library(tidyverse)

#------------------------------------------------------#
#            figure 5b                                 #
#------------------------------------------------------#

# read heterosis analysis results
bac.rs = read.table("./16S/heterosis/bac.rs.heterosis.sig.res.txt", 
                    sep = "\t", header = T)
bac.root = read.table("./16S/heterosis/bac.root.heterosis.significant.res.txt", 
                      sep = "\t", header = T)
# combine root and rhizosphere data
bac = rbind.data.frame(cbind(bac.rs, Comp="Rhizosphere"), cbind(bac.root, Comp="Root"))
# remove ASVs without family level 
bac = bac[bac$Family != " uncultured" & !is.na(bac$Family), ]
# remove families which has <5 ASVs
retained.fa = names(table(bac$Family))[which(table(bac$Family) >= 5)]
bac.df = bac[bac$Family %in% retained.fa, ]
bac.df$Paternal = factor(bac.df$Paternal, levels = rev(c("TZSTRI119", "TZMI722", "TZMI721", "H99", 
                                                     "CML108", "Oh43", "H84", "A554", "W64A", "Mo17")))
# add a "row" column which will be the y position in the plot
bac.df = bac.df %>% group_by(Paternal) 
bac.df = mutate(bac.df, row = group_indices())
# add a "col" column which will be the x position in the plot
bac.df = bac.df %>% group_by(Family) 
bac.df = mutate(bac.df, col = group_indices())
# barble plot
p = ggplot(bac.df, aes(x=factor(col), y=factor(row), color=Phylum)) +
  geom_point(alpha=0.5, size=3) +    # plot as points
  scale_x_discrete(breaks=1:length(unique(bac.df$Family)), 
                   labels=sort(unique(bac.df$Family)), 
                   position='top', guide = guide_axis(angle = 90)) +   # set the labels on the X axis
  scale_y_discrete(breaks=1:length(unique(bac.df$Paternal)), 
                   labels=levels(bac.df$Paternal)) +                 # set the labels on the Y axis
  theme_bw() +
  theme(axis.line = element_blank(),            # disable axis lines
        axis.title = element_blank(),           # disable axis titles
        panel.border = element_blank(),         # disable panel border
        panel.grid.major.x = element_blank(),   # disable lines in grid on X-axis
        panel.grid.minor.x = element_blank())   # disable lines in grid on X-axis
ggsave(p, filename = "./16S/plots/heterosis/heterosis-family-bac.pdf", height = 10, width = 10, units = "cm")

# do the same for fungal data
f.rs = read.table("./ITS/heterosis/fungi.rhizo.heterosis.sig.results.txt", sep = "\t", header = T)
f.root = read.table("./ITS/heterosis/fungi.root.heterosis.sig.results.txt", sep = "\t", header = T)

fungi = rbind.data.frame(cbind(f.rs, Comp="Rhizosphere"), cbind(f.root, Comp="Root"))
fungi = fungi[fungi$Family != " uncultured" & !is.na(fungi$Family), ]
retained.fa = names(table(fungi$Family))[which(table(fungi$Family) >= 5)]
fungi.df = fungi[fungi$Family %in% retained.fa, ]
fungi.df$Paternal = factor(fungi.df$Paternal, levels = rev(c("TZSTRI119", "TZMI722", "TZMI721", "H99", 
                                                         "CML108", "Oh43", "H84", "A554", "W64A", "Mo17")))
# add a "row" column which will be the y position in the plot
fungi.df <- fungi.df %>% group_by(Paternal) 
fungi.df = mutate(fungi.df, row = group_indices())
# add a "col" column which will be the x position in the plot
fungi.df <- fungi.df %>% group_by(Family) 
fungi.df = mutate(fungi.df, col = group_indices())

p = ggplot(fungi.df, aes(x=factor(col), y=factor(row), color=Phylum)) +
  geom_point(alpha=1, size=3) +    # plot as points
  #geom_text(aes(label=value, x=col + 0.25), alpha=1.0, size=3) +   # display the value next to the "balloons"
  #scale_size_area(max_size = 5) +
  scale_x_discrete(breaks=1:length(unique(fungi.df$Family)), labels=sort(unique(fungi.df$Family)), position='top', guide = guide_axis(angle = 90)) +   # set the labels on the X axis
  scale_y_discrete(breaks=1:length(unique(fungi.df$Paternal)), labels=levels(fungi.df$Paternal)) +                 # set the labels on the Y axis
  theme_bw() +
  theme(axis.line = element_blank(),            # disable axis lines
        axis.title = element_blank(),           # disable axis titles
        panel.border = element_blank(),         # disable panel border
        panel.grid.major.x = element_blank(),   # disable lines in grid on X-axis
        panel.grid.minor.x = element_blank())   # disable lines in grid on X-axis
ggsave(p, filename = "./ITS/plots/heterosis/heterosis-family-fungi.pdf", height = 10, width = 10, units = "cm")


####     bar plot of heterosis results      ####

# combine the data frame
het.df = rbind.data.frame(bac, fungi)
het.df = het.df[!is.na(het.df$Kingdom), ]
het.df$Paternal = factor(het.df$Paternal, levels = c("TZSTRI119", "TZMI722", "TZMI721", "H99", 
                                                             "CML108", "Oh43", "H84", "A554", "W64A", "Mo17"))
# add a "row" column which will be the y position in the plot
het.df <- het.df %>% group_by(Paternal) 
het.df = mutate(het.df, row = group_indices())
# count the number of heterotic ASVs in root and rhizosphere
het.df.bac = het.df[het.df$Kingdom=="Bacteria", ]
het.df.bac = het.df.bac %>% group_by(Comp, Paternal) %>% summarise(., Freq=n()) 
het.df.bac$Comp = factor(het.df.bac$Comp, levels = c("Root", "Rhizosphere"))

het.df.f = het.df[het.df$Kingdom=="Fungi", ]
het.df.f = het.df.f %>% group_by(Comp, Paternal) %>% summarise(., Freq=n()) 
het.df.f$Comp = factor(het.df.f$Comp, levels = c("Root", "Rhizosphere"))
# bar plot for bacteria and fungi respectively
p = ggplot(het.df.f, aes(x=Paternal, y=Freq, fill=Comp)) +
    geom_bar(stat = "identity", alpha=1, width = 0.3) +     
    scale_fill_manual(values = c( "chocolate1", "coral4"))+
    mytheme +
    ylim(c(0, 40))
ggsave(p, filename = "./ITS/plots/heterosis/heterosis-compartment-fungi.pdf", 
       height = 9, width = 15, units = "cm")

mytheme = theme_bw() +
    theme(plot.margin = unit(c(1,1,1,1), "cm"), plot.title = element_text(size = 10),
          panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
          axis.text = element_text(colour = "black", size = 8)) 



