## to run the code please first cuctomize the following file directory

#setwd("~/GWAS/")
#datapwd <- "~/Input_data/"

library(BGLR)

source("GWAS_functions.R")

###------------------------------------------------------
#      GWAS for ionome data
###------------------------------------------------------

phdata <- read.table(paste0(datapwd,"PHdata_ionome.txt"),header=TRUE)
dim(phdata)
head(phdata)
rownames(phdata) <- phdata$Sample.ID

traitname <- colnames(phdata)[-1:-6]
nt <- length(traitname)

treatmentname <- levels(phdata$Treatment)
ntr <- nlevels(phdata$Treatment)

GNdata <- read.table(paste0(datapwd,"SNP_matrix_imputed_filtered.txt"),header=TRUE)
GNdata <- t(as.matrix(GNdata))
ngen <- nrow(GNdata)
nmar <- ncol(GNdata)
genoname <- rownames(GNdata)
marname <- colnames(GNdata)

MKinfo <- read.table(paste0(datapwd,"SNP_info_filtered.txt"),header=TRUE)

Kmat <- read.table(paste0(datapwd,"VanRaden_G_matrix_marker_filtered.txt"),header=TRUE)
Kmat <- as.matrix(Kmat)

Xmat <- read.table(paste0(datapwd,"Design_matrix_germplasm_effect.txt"),header=TRUE)
Xmat <- as.matrix(Xmat)

herit <- read.table(paste0(datapwd,"Result_heritability_Ionome.txt"),header=TRUE)

for (i in 1:nt)
{
  PHdata <- read.table(paste0(datapwd,"BLUE_",traitname[i],".txt"),header=TRUE)
  for (j in 1:ntr)
  {
    if (herit[traitname[i],treatmentname[j]] < 0.1)
    {
      next
    }else{
      gwas.res <- FastGWAS_P3D_basic_v2(y=PHdata[genoname,treatmentname[j]],
                                        M=GNdata[genoname,],
                                        G=Kmat[genoname,genoname],
                                        X=Xmat[genoname,])
      write.table(gwas.res,paste0("Result_GWAS_",traitname[i],"_",treatmentname[j],".txt"),quote=FALSE)
    }
  }
}
