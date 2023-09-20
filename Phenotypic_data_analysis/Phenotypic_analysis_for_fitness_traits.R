###----------------------------------------------------------
#       Analyzing the phenotypic data for fitness traits
###----------------------------------------------------------


## to run the code please first cuctomize the following file directory

# setwd("~/code_for_the_manuscript/Phenotypic_data_analysis/")
# datapwd <- "~/code_for_the_manuscript/Input_data/"


library(asreml)

phdata <- read.table(paste0(datapwd,"PHdata_fitness.txt"),header=TRUE)
dim(phdata)
head(phdata)
rownames(phdata) <- phdata$Sample.ID

traitname <- colnames(phdata)[7:9]
nt <- length(traitname)

treatmentname <- levels(phdata$Treatment)
ntr <- nlevels(phdata$Treatment)

ngen <- nlevels(phdata$Genotype)
gen_name <- levels(phdata$Genotype)

nrep <- nlevels(as.factor(phdata$Block))


### estimating the heritability within each treatment (repeatability)

repetability <- matrix(0,nt,ntr)
rownames(repetability) <- traitname
colnames(repetability) <- c(treatmentname)

for (i in 1:nt)
{
  for (j in 1:(ntr))
  {
  
      subphdata <- phdata[which(phdata$Treatment==treatmentname[j]),]
      subphdata$Block <- factor(subphdata$Block)
      subphdata$Plot <- factor(subphdata$Plot)
      subphdata$Genotype <- factor(subphdata$Genotype)
      
      fitmm_r <- asreml(fixed = as.formula(paste(traitname[i],"~ 1")),
                        random = ~ Genotype + Block + at(Block):Plot,
                        data = subphdata)
      
      vc <- summary(fitmm_r)$varcomp
   
      var_G <- vc["Genotype","component"]
      var_e <- vc["units!R","component"]
      
      repetability[i,j] <- var_G/(var_G+var_e/nrep)
  }
}

write.table(repetability,"Result_heritability_fitness.txt",quote=FALSE)



### outlier test

res_outlier <- NULL

for (i in 1:nt)
{
  for (j in 1:ntr)
  {
    subphdata <- phdata[which(phdata$Treatment==treatmentname[j]),]
    subphdata$Block <- factor(subphdata$Block)
    subphdata$Plot <- factor(subphdata$Plot)
    subphdata$Germplasm <- factor(subphdata$Germplasm)
    fitmm <- asreml(fixed = as.formula(paste(traitname[i],"~ 1+ Germplasm")),
                    random = ~ Genotype + Block + at(Block):Plot,
                    data = subphdata)
    resd <- fitmm$residuals
    MAD <- median(abs(resd-median(resd)))
    s_r <- 1.4826*MAD
    resd_m <- resd/s_r
    pval <- 2*pnorm(-abs(resd_m))
    adjpval <- p.adjust(pval,method="holm")
    res_table <- data.frame(subphdata[,c("Genotype","Germplasm","Treatment","Block","Plot")],
                            subphdata[,traitname[i]],resd_m,pval,adjpval)
    colnames(res_table)[6:9] <- c("Trait_value","Res_std","P","P_adj")
    
    pos <- which(adjpval < 0.05)
    if (length(pos) > 0)
    {
      if (is.null(res_outlier))
      {
        res_outlier <- data.frame(res_table[pos,],Trait=rep(traitname[i],length(pos)))
      }else{
        add_outlier <- data.frame(res_table[pos,],Trait=rep(traitname[i],length(pos)))
        res_outlier <- rbind(res_outlier,add_outlier)
      }
    }
    
    write.table(res_table,paste0("Result_details_outlier_test_",traitname[i],"_",treatmentname[j],".txt"),quote=FALSE)
  } 
}

write.table(res_outlier,"Result_outliers_all_trait_and_treatment.txt",quote=FALSE)


# re-estimate the repeatability after remove all outliers

ot <- read.table("Result_outliers_all_trait_and_treatment.txt",header=TRUE)
str(ot)
head(ot)
two <- levels(ot$Trait)

repetability <- matrix(0,nt,ntr)
rownames(repetability) <- traitname
colnames(repetability) <- c(treatmentname)

for (i in 1:nt)
{
  if (traitname[i] %in% two)
  {
    samid <- rownames(ot)[which(ot$Trait==traitname[i])] 
    pos <- which(rownames(phdata)%in%samid)
    phdata_o <- phdata[-pos,]
  }else{
    phdata_o <- phdata
  }
  
  for (j in 1:ntr)
  {
      subphdata <- phdata_o[which(phdata_o$Treatment==treatmentname[j]),]
      subphdata$Block <- factor(subphdata$Block)
      subphdata$Plot <- factor(subphdata$Plot)
      subphdata$Genotype <- factor(subphdata$Genotype)
      
      fitmm_r <- asreml(fixed = as.formula(paste(traitname[i],"~ 1")),
                        random = ~ Genotype + Block + at(Block):Plot,
                        data = subphdata)
      
      vc <- summary(fitmm_r)$varcomp
  
      var_G <- vc["Genotype","component"]
      var_e <- vc["units!R","component"]
      
      repetability[i,j] <- var_G/(var_G+var_e/nrep)
  }
}

write.table(repetability,"Result_heritability_outlier_removed.txt",quote=FALSE)


### The BLUEs

for (i in 1:nt)
{
  if (traitname[i] %in% two)
  {
    samid <- rownames(ot)[which(ot$Trait==traitname[i])] 
    pos <- which(rownames(phdata)%in%samid)
    phdata_o <- phdata[-pos,]
  }else{
    phdata_o <- phdata
  }
  
  blue_table <- matrix(0,ngen,ntr)
  rownames(blue_table) <- gen_name
  colnames(blue_table) <- treatmentname
  
  for (j in 1:ntr)
  {
      subphdata <- phdata_o[which(phdata_o$Treatment==treatmentname[j]),]
      subphdata$Block <- factor(subphdata$Block)
      subphdata$Plot <- factor(subphdata$Plot)
      subphdata$Genotype <- factor(subphdata$Genotype)
      
      fitmm_f <- asreml(fixed = as.formula(paste(traitname[i],"~ Genotype")),
                        random = ~ Block + at(Block):Plot,
                        data = subphdata)
      blue <- predict.asreml(fitmm_f,classify="Genotype")$pvals
      rownames(blue) <- as.character(blue$Genotype)
      
      blue_table[,j] <- blue[gen_name,"predicted.value"] 
  }
  write.table(blue_table,paste0("BLUE_",traitname[i],".txt"),quote=FALSE)
}




