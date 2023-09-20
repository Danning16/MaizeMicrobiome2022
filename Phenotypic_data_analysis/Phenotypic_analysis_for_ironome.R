## to run the code please first cuctomize the following file directory

#setwd("~/Phenotypic_data_analysis/")
#datapwd <- "~/Input_data/"

library(asreml)


###-------------------------------------------------
#       Analyzing the phenotypic data of ionome
###-------------------------------------------------

phdata <- read.table(paste0(datapwd,"PHdata_ionome.txt"),header=TRUE)
dim(phdata)
head(phdata)
rownames(phdata) <- phdata$Sample.ID

str(phdata)

traitname <- colnames(phdata)[-1:-6]
nt <- length(traitname)

treatmentname <- levels(phdata$Treatment)
ntr <- nlevels(phdata$Treatment)

ngen <- nlevels(phdata$Genotype)
gen_name <- levels(phdata$Genotype)

nrep <- nlevels(as.factor(phdata$Block))



### estimating the heritability within each treatment (repeatability)

repetability <- matrix(0,nt,ntr)
rownames(repetability) <- traitname
colnames(repetability) <- treatmentname

for (i in 1:nt)
{
  for (j in 1:ntr)
  {
      subphdata <- phdata[which(phdata$Treatment==treatmentname[j]),]
      subphdata$Block <- factor(subphdata$Block)
      subphdata$Plot <- factor(subphdata$Plot)
      subphdata$Genotype <- factor(subphdata$Genotype)
      nbk <- nlevels(subphdata$Block)
      
      fitmm_r <- asreml(fixed = as.formula(paste(traitname[i],"~ 1")),
                        random = ~ Genotype + Block + at(Block):Plot,
                        data = subphdata)
      
      vc <- summary(fitmm_r)$varcomp
      write.table(vc,paste0("varcomp_",traitname[i],"_",treatmentname[j],".txt"),quote=FALSE)
      
      var_G <- vc["Genotype","component"]
      var_e <- vc["units!R","component"]
      
      repetability[i,j] <- var_G/(var_G+var_e/nrep)
  }
}

write.table(repetability,paste0("Result_heritability_Ironome.txt"),quote=FALSE)



### The BLUEs

for (i in 1:nt)
{
  phdata_o <- phdata
  
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


# re-arrange the tables

for (i in 1:ntr)
{
  res <- NULL
  for (j in 1:nt)
  {
    tempres <- read.table(paste0("BLUE_",traitname[j],".txt"),header=TRUE)
    if (is.null(res))
    {
      res <- tempres[,treatmentname[i],drop=FALSE]
    }else{
      res <- cbind(res,tempres[,treatmentname[i],drop=FALSE])
    }
  }
  colnames(res) <- traitname
  write.table(res,paste0("BLUE_Ironome_",treatmentname[i],".txt"),quote=FALSE)
}