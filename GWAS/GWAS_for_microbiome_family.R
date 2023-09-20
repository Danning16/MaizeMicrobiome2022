## to run the code please first cuctomize the following file directory

#setwd("~/GWAS/")
#datapwd <- "~/Input_data/"

library(BGLR)

source("GWAS_functions.R")

###------------------------------------------------------
#      GWAS for microbiome family data
###------------------------------------------------------


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

microbiometype <- c("fungi","bacteria")
nmicro <- length(microbiometype)

sampletype <- c("root","rhizosphere")
nsample <- length(sampletype)

treatmentname <- c("CK","D","LN","LP")
ntreatment <- length(treatmentname)

for (i in 1:nmicro)
{
  for (j in 1:nsample)
  {
    herit <- read.table(paste0(datapwd,"Result_heritability_",microbiometype[i],"_family_",sampletype[j],".txt"),header=TRUE)
    for (k in 1:ntreatment)
    {
      PHdata <- read.table(paste0(datapwd,"BLUE_",microbiometype[i],"_family_",sampletype[j],"_",treatmentname[k],".txt"),header=TRUE)
      vec_herit <- herit[,treatmentname[k]]
      pos <- which(vec_herit >= 0.1)
      asvname_filtered <- rownames(herit)[pos]
      nasv <- length(asvname_filtered)
      
      pval_gwas <- matrix(0,nmar,nasv)
      colnames(pval_gwas) <- asvname_filtered
      rownames(pval_gwas) <- marname
      
      for (l in 1:nasv)
      {
        gwas.res <- FastGWAS_P3D_basic_v2(y=PHdata[genoname,asvname_filtered[l]],
                                          M=GNdata[genoname,],
                                          G=Kmat,
                                          X=Xmat)
        pval_gwas[,asvname_filtered[l]] <- gwas.res[marname,"P_value"]
        cat("GWAS for",microbiometype[i],sampletype[j],treatmentname[k],"family",l,"completed\n")
      } 
      
      write.table(pval_gwas,paste0("Result_GWAS_",microbiometype[i],"_family_",sampletype[j],"_",treatmentname[k],".txt"),quote=FALSE)
      
    } # end loop k
  } # end loop j
} # end loop i




