## to run the code please first cuctomize the following file directory

#setwd("~/GWAS/")
#datapwd <- "~/Input_data/"

library(BGLR)

source("GWAS_functions.R")

###------------------------------------------------------
#      GWAS for fitness traits
###------------------------------------------------------

phinfo <- read.table(paste0(datapwd,"PHdata_fitness.txt"),header=TRUE)

traitname <- colnames(phinfo)[7:9]
nt <- length(traitname)

treatmentname <- levels(phinfo$Treatment)
ntr <- nlevels(phdata$Treatment)

GNdata <- read.table(paste0(datapwd,"SNP_matrix_imputed_filtered.txt"),header=TRUE)
dim(GNdata)
GNdata <- t(as.matrix(GNdata))
ngen <- nrow(GNdata)
nmar <- ncol(GNdata)
genoname <- rownames(GNdata)
marname <- colnames(GNdata)

MKinfo <- read.table(paste0(datapwd,"SNP_info_filtered.txt"),header=TRUE)
dim(MKinfo)
head(MKinfo)

# generating the VanRaden G-matrix

Kmat <- GetVanRadenKinMat(GNdata)
dim(Kmat)
Kmat[1:5,1:5]
mean(diag(Kmat))

Kmat_new <- MakeSemiPDmatrix(Kmat)
rownames(Kmat_new) <- genoname
colnames(Kmat_new) <- genoname
write.table(Kmat_new,"VanRaden_G_matrix_marker_filtered.txt",quote=FALSE)


# generating the coefficient matrix for the covariates (subpopulation effects)

info_germplasm <- rep("N",ngen)
for(i in 1:ngen)
{
  subdta <- phinfo[which(phinfo$Genotype==genoname[i]),]
  subdta$Germplasm <- factor(subdta$Germplasm)
  info_germplasm[i] <- levels(subdta$Germplasm)
}
res <- data.frame(Genotype=genoname,Germplasm=info_germplasm)
write.table(res,"Germplasm_info.txt",quote=FALSE)

ngerm <- nlevels(res$Germplasm)
germ_name <- levels(res$Germplasm)

X <- matrix(0,ngen,ngerm-1)
rownames(X) <- genoname
colnames(X) <- germ_name[-ngerm]
for (i in 1:(ngerm-1))
{
  X[which(res$Germplasm==germ_name[i]),i] <- 1
}
write.table(X,"Design_matrix_germplasm_effect.txt",quote=FALSE)


# performing GWAS

for (i in 1:nt)
{
  PHdata <- read.table(paste0(datapwd,"BLUE_",traitname[i],".txt"),header=TRUE)
  for (j in 1:ntr)
  {
    gwas.res <- FastGWAS_P3D_basic_v2(y=PHdata[genoname,treatmentname[j]],
                                   M=GNdata[genoname,],
                                   G=Kmat_new,
                                   X=X)
    write.table(gwas.res,paste0("Result_GWAS_",traitname[i],"_",treatmentname[j],".txt"),quote=FALSE)
  }
}

