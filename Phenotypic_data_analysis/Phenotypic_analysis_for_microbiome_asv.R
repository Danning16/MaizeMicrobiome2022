## to run the code please first cuctomize the following file directory

#setwd("~/Phenotypic_data_analysis/")
#datapwd <- "~/Input_data/"

library(asreml)

###--------------------------------
#  analyzing the fungi asv data
###--------------------------------


## the root sample

FUdata_r <- read.table(paste0(datapwd,"Fungi_data_root.txt"),header=TRUE)

notu <- ncol(FUdata_r)-6
otuname <- colnames(FUdata_r)[-1:-6]

FUdata_r$Germplasm <- factor(FUdata_r$Germplasm)
FUdata_r$Treatment <- factor(FUdata_r$Treatment)
FUdata_r$Compartment <- factor(FUdata_r$Compartment)
FUdata_r$Block <- factor(FUdata_r$Block)
FUdata_r$Plot <- factor(FUdata_r$Plot)

treatmentname <- levels(FUdata_r$Treatment)
ntr <- nlevels(FUdata_r$Treatment)

ngen <- nlevels(FUdata_r$Genotype)
gen_name <- levels(FUdata_r$Genotype)

nrep <- nlevels(as.factor(FUdata_r$Block))


# estimating the heritability

repetability <- matrix(0,notu,ntr)
rownames(repetability) <- otuname
colnames(repetability) <- ctreatmentname

for (i in 1:notu)
{
  for (j in 1:ntr)
  {
      subphdata <- FUdata_r[which(FUdata_r$Treatment==treatmentname[j]),]
      subphdata$Block <- factor(subphdata$Block)
      subphdata$Plot <- factor(subphdata$Plot)
      fitmm_r <- asreml(fixed = as.formula(paste(otuname[i],"~ 1")),
                        random = ~ Genotype + Block + at(Block):Plot,
                        data = subphdata)
      
      vc <- summary(fitmm_r)$varcomp
   
      var_G <- vc["Genotype","component"]
      var_e <- vc["units!R","component"]
      
      repetability[i,j] <- var_G/(var_G+var_e/nrep)
  }
}

write.table(repetability,"Result_heritability_fungi_asv_root.txt",quote=FALSE)


# BLUPs

for (j in 1:ntr)
{
  blue_mat <- matrix(0,ngen,notu)
  rownames(blue_mat) <- gen_name
  colnames(blue_mat) <- otuname
  
  subphdata <- FUdata_r[which(FUdata_r$Treatment==treatmentname[j]),]
  subphdata$Block <- factor(subphdata$Block)
  subphdata$Plot <- factor(subphdata$Plot)
  
  for (i in 1:notu)
  {

    fitmm_r <- asreml(fixed = as.formula(paste(otuname[i],"~ 1")),
                      random = ~ Genotype + Block + at(Block):Plot,
                      data = subphdata)
 
    blue <- predict(fitmm_r,classify="Genotype")$pvals
    rownames(blue) <- as.character(blue$Genotype)
    blue_mat[,i] <- blue[gen_name,"predicted.value"]
  }

    write.table(blue_mat,paste0("BLUE_fungi_asv_root_",treatmentname[j],".txt"),quote=FALSE)

}



## the rhizoshere sample

FUdata_s <- read.table(paste0(datapwd,"Fungi_data_rhizosphere.txt"),header=TRUE)

notu <- ncol(FUdata_s)-6
otuname <- colnames(FUdata_s)[-1:-6]

FUdata_s$Germplasm <- factor(FUdata_s$Germplasm)
FUdata_s$Treatment <- factor(FUdata_s$Treatment)
FUdata_s$Compartment <- factor(FUdata_s$Compartment)
FUdata_s$Block <- factor(FUdata_s$Block)
FUdata_s$Plot <- factor(FUdata_s$Plot)

treatmentname <- levels(FUdata_s$Treatment)
ntr <- nlevels(FUdata_s$Treatment)

ngen <- nlevels(FUdata_s$Genotype)
gen_name <- levels(FUdata_s$Genotype)

nrep <- nlevels(as.factor(FUdata_s$Block))


# estimating the heritability

repetability <- matrix(0,notu,ntr)
rownames(repetability) <- otuname
colnames(repetability) <- treatmentname

for (i in 1:notu)
{
  for (j in 1:ntr)
  {
      subphdata <- FUdata_s[which(FUdata_s$Treatment==treatmentname[j]),]
      subphdata$Block <- factor(subphdata$Block)
      subphdata$Plot <- factor(subphdata$Plot)
      subphdata$Genotype <- factor(subphdata$Genotype)
      fitmm_r <- asreml(fixed = as.formula(paste(otuname[i],"~ 1")),
                        random = ~ Genotype + Block + at(Block):Plot,
                        data = subphdata)
      
      vc <- summary(fitmm_r)$varcomp
    
      var_G <- vc["Genotype","component"]
      var_e <- vc["units!R","component"]
      
      repetability[i,j] <- var_G/(var_G+var_e/nrep)
  }
}

write.table(repetability,"Result_heritability_fungi_asv_rhizosphere.txt",quote=FALSE)


# BLUPs

for (j in 1:ntr)
{
  blue_mat <- matrix(0,ngen,notu)
  rownames(blue_mat) <- gen_name
  colnames(blue_mat) <- otuname
  
  subphdata <- FUdata_s[which(FUdata_s$Treatment==treatmentname[j]),]
  subphdata$Block <- factor(subphdata$Block)
  subphdata$Plot <- factor(subphdata$Plot)
  
  for (i in 1:notu)
  {
    fitmm_s <- asreml(fixed = as.formula(paste(otuname[i],"~ 1")),
                        random = ~ Genotype + Block + at(Block):Plot,
                        data = subphdata)

    
    blue <- predict(fitmm_s,classify="Genotype")$pvals
    rownames(blue) <- as.character(blue$Genotype)
    blue_mat[,i] <- blue[gen_name,"predicted.value"]
  }
  write.table(blue_mat,paste0("BLUE_fungi_asv_rhizosphere_",treatmentname[j],".txt"),quote=FALSE)
}



###--------------------------------
# analyzing the bacteria asv data
###--------------------------------

## the root sample

BAdata_r <- read.table(paste0(datapwd,"Bacteria_data_root.txt"),header=TRUE)

notu <- ncol(BAdata_r)-6
otuname <- colnames(BAdata_r)[-1:-6]

BAdata_r$Germplasm <- factor(BAdata_r$Germplasm)
BAdata_r$Treatment <- factor(BAdata_r$Treatment)
BAdata_r$Compartment <- factor(BAdata_r$Compartment)
BAdata_r$Block <- factor(BAdata_r$Block)
BAdata_r$Plot <- factor(BAdata_r$Plot)

treatmentname <- levels(BAdata_r$Treatment)
ntr <- nlevels(BAdata_r$Treatment)

ngen <- nlevels(BAdata_r$Genotype)
gen_name <- levels(BAdata_r$Genotype)

nrep <- nlevels(BAdata_r$Block)

# estimating the heritability

repetability <- matrix(0,notu,ntr)
rownames(repetability) <- otuname
colnames(repetability) <- treatmentname

for (i in 1:notu)
{
  for (j in 1:ntr)
  {
      subphdata <- BAdata_r[which(BAdata_r$Treatment==treatmentname[j]),]
      subphdata$Block <- factor(subphdata$Block)
      subphdata$Plot <- factor(subphdata$Plot)
      fitmm_r <- asreml(fixed = as.formula(paste(otuname[i],"~ 1")),
                        random = ~ Genotype + Block + at(Block):Plot,
                        data = subphdata)
      
      vc <- summary(fitmm_r)$varcomp
     
      var_G <- vc["Genotype","component"]
      var_e <- vc["units!R","component"]
      
      repetability[i,j] <- var_G/(var_G+var_e/nrep)
  }
}

write.table(repetability,"Result_heritability_bacteria_asv_root.txt",quote=FALSE)


# the BLUPs

for (j in 1:ntr)
{
  blue_mat <- matrix(0,ngen,notu)
  rownames(blue_mat) <- gen_name
  colnames(blue_mat) <- otuname
  
  subphdata <- BAdata_r[which(BAdata_r$Treatment==treatmentname[j]),]
  subphdata$Block <- factor(subphdata$Block)
  subphdata$Plot <- factor(subphdata$Plot)
  
  for (i in 1:notu)
  {
    fitmm_r <- asreml(fixed = as.formula(paste(otuname[i],"~ 1")),
                      random = ~ Genotype + Block + at(Block):Plot,
                      data = subphdata)
    
    blue <- predict(fitmm_r,classify="Genotype")$pvals
    rownames(blue) <- as.character(blue$Genotype)
    blue_mat[,i] <- blue[gen_name,"predicted.value"]
  }
  
  write.table(blue_mat,paste0("BLUE_bacteria_asv_root_",treatmentname[j],".txt"),quote=FALSE)
}


## The rhizosphere sample


BAdata_s <- read.table(paste0(datapwd,"Bacteria_data_rhizosphere.txt"),header=TRUE)

notu <- ncol(BAdata_s)-6
otuname <- colnames(BAdata_s)[-1:-6]

BAdata_s$Germplasm <- factor(BAdata_s$Germplasm)
BAdata_s$Treatment <- factor(BAdata_s$Treatment)
BAdata_s$Compartment <- factor(BAdata_s$Compartment)
BAdata_s$Block <- factor(BAdata_s$Block)
BAdata_s$Plot <- factor(BAdata_s$Plot)

treatmentname <- levels(BAdata_s$Treatment)
ntr <- nlevels(BAdata_s$Treatment)

ngen <- nlevels(BAdata_s$Genotype)
gen_name <- levels(BAdata_s$Genotype)

nrep <- nlevels(BAdata_s$Block)

# the heritability

repetability <- matrix(0,notu,ntr)
rownames(repetability) <- otuname
colnames(repetability) <- treatmentname

for (i in 1:notu)
{
  for (j in 1:ntr)
  {
      subphdata <- BAdata_s[which(BAdata_s$Treatment==treatmentname[j]),]
      subphdata$Block <- factor(subphdata$Block)
      subphdata$Plot <- factor(subphdata$Plot)
      subphdata$Genotype <- factor(subphdata$Genotype)
      fitmm_r <- asreml(fixed = as.formula(paste(otuname[i],"~ 1")),
                        random = ~ Genotype + Block + at(Block):Plot,
                        data = subphdata)
      
      vc <- summary(fitmm_r)$varcomp
     
      var_G <- vc["Genotype","component"]
      var_e <- vc["units!R","component"]
      
      repetability[i,j] <- var_G/(var_G+var_e/nrep)
  }
}

write.table(repetability,"Result_heritability_bacteria_asv_rhizosphere.txt",quote=FALSE)

# the BLUPs

for (j in 1:ntr)
{
  blue_mat <- matrix(0,ngen,notu)
  rownames(blue_mat) <- gen_name
  colnames(blue_mat) <- otuname
  
  subphdata <- BAdata_s[which(BAdata_s$Treatment==treatmentname[j]),]
  subphdata$Block <- factor(subphdata$Block)
  subphdata$Plot <- factor(subphdata$Plot)
  
  for (i in 1:notu)
  {
    fitmm_s <- asreml(fixed = as.formula(paste(otuname[i],"~ 1")),
                        random = ~ Genotype + Block + at(Block):Plot,
                        data = subphdata)
    
    blue <- predict(fitmm_s,classify="Genotype")$pvals
    rownames(blue) <- as.character(blue$Genotype)
    blue_mat[,i] <- blue[gen_name,"predicted.value"]
  }
    write.table(blue_mat,paste0("BLUE_bacteria_asv_rhizosphere_",treatmentname[j],".txt"),quote=FALSE)
}

