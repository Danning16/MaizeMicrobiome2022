## to run the code please first cuctomize the following file directory

#setwd("~/GWAS/")
#datapwd <- "~/Input_data/"

library(GMMAT)

#---------------------------------------------------------------------------
#   GWAS with presense-absence mode for microbiome ASVs
#---------------------------------------------------------------------------

## transfering the phenotypic data into the PA mode

microbiometype <- c("fungi","bacteria")
nmicro <- length(microbiometype)

sampletype <- c("root","rhizosphere")
nsample <- length(sampletype)

treatmentname <- c("CK","D","LN","LP")
ntreatment <- length(treatmentname)

checkpa <- function(x){
  y <- (x > 0)
  z <- length(which(y==TRUE))
  return(ifelse(z>1,1,0))
}

for (i in 1:nmicro)
{
  for (j in 1:nsample)
  {
    rawdta <- read.table(paste0(datapwd,microbiometype[i],"_data_",sampletype[j],".txt"),header=TRUE)

    for (k in 1:ntreatment)
    {
      subdta <- rawdta[which(rawdta$Treatment==treatmentname[k],),]
      dim(subdta)
      asvname <- colnames(subdta)[-1:-6]
      nasv <- length(asvname)
      
      gypname <- levels(subdta$Genotype)
      ngen <- length(gypname)
      
      dta <- list()
      for (s in 1:3)
      {
        dta[[s]] <- subdta[which(subdta$Block==s),]
        rownames(dta[[s]]) <- as.character(dta[[s]]$Genotype)
        dta[[s]] <- dta[[s]][gypname,]
      }
      
      phdata_pa <- matrix(0,ngen,nasv)
      rownames(phdata_pa) <- gypname
      colnames(phdata_pa) <- asvname
      
      for (l in 1:nasv)
      {
        cdta <- dta[[1]][,asvname[l],drop=FALSE]
        for (s in 2:3)
        {
          cdta <- cbind(cdta,dta[[s]][,asvname[l],drop=FALSE])
        }
        dim(cdta)
        head(cdta)
        
        phdata_pa[,asvname[l]] <- apply(cdta,1,checkpa)
      } # end loop l
      
      write.table(phdata_pa,paste0("PAdata_",microbiometype[i],"_asv_",sampletype[j],"_",treatmentname[k],".txt"),quote=FALSE)
    } # end loop k
    
  } # end loop j
  
} # end loop i


## performing PA-GWAS using the GMMAT package

Info_germ <- read.table(paste0(datapwd,"Germplasm_info.txt"),header=TRUE)
head(Info_germ)
rownames(Info_germ) <- as.character(Info_germ$Genotype)

geno_name <- rownames(Info_germ)
ngen <- length(geno_name)
Kmat <- Kmat[geno_name,geno_name]

GNdata <- read.table(paste0(datapwd,"SNP_matrix_imputed_filtered.txt"),header=TRUE)
GNdata <- GNdata[,geno_name]
GNdata[1:5,1:5]
write.table(GNdata,"input_MKdata_for_PAGWAS.txt",quote=FALSE,sep="\t",
            row.names = FALSE,col.names = FALSE)
marname <- rownames(GNdata)
nmar <- length(marname)

MKinfo <- read.table(paste0(datapwd,"SNP_info_filtered.txt"),header=TRUE)

Kmat <- read.table(paste0(datapwd,"VanRaden_G_matrix_marker_filtered.txt"),header=TRUE)
Kmat <- as.matrix(Kmat)

Xmat <- read.table(paste0(datapwd,"Design_matrix_germplasm_effect.txt"),header=TRUE)
Xmat <- as.matrix(Xmat)

for (i in 1:nmicro)
{
  for (j in 1:nsample)
  {
    
    herit <- read.table(paste0(datapwd,"Result_heritability_",microbiometype[i],"_asv_",sampletype[j],".txt"),header=TRUE)
    colnames(herit)[5] <- "across_treatment"
    
    for (k in 1:ntreatment)
    {
      PHdata <- read.table(paste0("PAdata_",microbiometype[i],"_asv_",sampletype[j],"_",treatmentname[k],".txt"),header=TRUE)
      ckfreq <- apply(as.matrix(PHdata),2,sum)/ngen
      vec_herit <- herit[,treatmentname[k]]
      pos <- which(vec_herit >= 0.1 & (ckfreq < 0.95 & ckfreq > 0.05))
      asvname_filtered <- rownames(herit)[pos]
      nasv <- length(asvname_filtered)
      
      res <- matrix(0,nmar,nasv)
      rownames(res) <- marname
      colnames(res) <- asvname_filtered
      
      for (l in 1:nasv)
      {
        pheno <- cbind(Info_germ,PHdata[geno_name,asvname_filtered[l]])
        colnames(pheno)[3] <- "ASV"
        
        tryCatch({model0 <- glmmkin(ASV ~ Germplasm, data = pheno, kins = Kmat, id = "Genotype",
                          family = binomial(link = "logit"))
        
        infile <- "input_MKdata_for_PAGWAS.txt"
        outfile <- paste0("GMMAT_score_test",microbiometype[i],"_",sampletype[j],"_",
                          treatmentname[k],"_",asvname_filtered[l],".txt")
        glmm.score(model0, infile = infile, outfile = outfile, infile.nrow.skip = 0,
                   infile.ncol.skip = 0, 
                   infile.ncol.print = 0)
        
        tempres <- read.table(outfile,header=TRUE)
        head(tempres)
        
        res[,asvname_filtered[l]] <- tempres$PVAL},
        
        error=function(e){cat("Error\n")
          res[,asvname_filtered[l]] <- NA
        }
        )
        
        cat("GMMAT score test for",microbiometype[i],sampletype[j],
            treatmentname[k],asvname_filtered[l],"completed\n")
      }# end loop l 
      
      write.table(res,paste0("result_PAGWAS_",microbiometype[i],"_asv_",sampletype[j],"_",
                             treatmentname[k],"_CofasFac.txt"),quote=FALSE)
      
    } # end loop k
    
  } # end loop j
  
} # end loop i

