## to run the code please first cuctomize the following file directory

#setwd("~/Genomic_prediction/")
#datapwd <- "~/Input_data/"

library(BGLR)


###----------------------------------------------------------------------------------
#      Prediction for ionomes using genomic data and microbiome ASVs
###----------------------------------------------------------------------------------

microbiometype <- c("fungi","bacteria")
nmicro <- length(microbiometype)

sampletype <- c("root","rhizosphere")
nsample <- length(sampletype)

treatmentname <- c("CK","D","LN","LP")
ntreatment <- length(treatmentname)

Kmat <- read.table(paste0(datapwd,"VanRaden_G_matrix_marker_filtered_for_GP.txt"),header=TRUE)
Kmat <- as.matrix(Kmat)

name_geno <- rownames(Kmat)
ngen <- length(name_geno)

herit_ph <- read.table(paste0(datapwd,"Result_heritability_Ionome.txt"),header=TRUE)
traitname <- rownames(herit_ph)
ntrait <- length(traitname)

sdd <- function(x){
  y <- x-mean(x)
  return(y/sd(x))
}

scenario_name_1 <- c("BA_RO","BA_RH","FU_RO","FU_RH","BA","FU","RO","RH","ALL")
scenario_name_2 <- paste("G",scenario_name_1,sep="_")
scenario_name <- c(scenario_name_1,scenario_name_2)
nscenario <- length(scenario_name)

res_pa <- matrix(0,ntrait*ntreatment,nscenario)
colnames(res_pa) <- scenario_name
rownames(res_pa) <- paste(rep(traitname,each=ntreatment),rep(treatmentname,ntrait),sep="_")

res_sd <- res_pa

for (l in 1:ntrait)
{
  PHdata <- read.table(paste0("BLUE_",traitname[l],".txt"),header=TRUE)
  phdata <- PHdata[name_geno,]
  
  for (k in 1:ntreatment)
  {
    
    y=phdata[,k,drop=TRUE]
    ngen <- length(y)
    
    Kmat_micro <- list()
    for (i in 1:nmicro)
    {
      
      Kmat_micro[[i]] <- list()
      
      for (j in 1:nsample)
      {
        herit <- read.table(paste0("Result_heritability_",microbiometype[i],"_asv_",sampletype[j],".txt"),header=TRUE)
        Mdata <- read.table(paste0("BLUE_",microbiometype[i],"_asv_",sampletype[j],"_",treatmentname[k],".txt"),header=TRUE)
        
        pos <- which(herit[,k]>0.1)
        Mdata_r <- Mdata[,pos]
        Mdata_r <- Mdata_r[name_geno,]
        
        Mdata_sdd <- apply(Mdata_r,2,sdd)
        Gm <- Mdata_sdd %*% t(Mdata_sdd)
        Gm <- Gm/mean(diag(Gm))
        
        Kmat_micro[[i]][[j]] <- Gm
      } # end loop j
      
    } # end loop i
    
    names(Kmat_micro) <- microbiometype
    for (i in 1:nmicro)
    {
      names(Kmat_micro[[i]]) <- sampletype
    } # end loop i
    
    ETAlist <- list()
    ETAlist[[1]] <- list(list(K=Kmat_micro$bacteria$root,model="RKHS"))
    ETAlist[[2]] <- list(list(K=Kmat_micro$bacteria$rhizosphere,model="RKHS"))
    ETAlist[[3]] <- list(list(K=Kmat_micro$fungi$root,model="RKHS"))
    ETAlist[[4]] <- list(list(K=Kmat_micro$fungi$rhizosphere,model="RKHS"))
    ETAlist[[5]] <- list(list(K1=Kmat_micro$bacteria$root,model="RKHS"),
                         list(K2=Kmat_micro$bacteria$rhizosphere,model="RKHS"))
    ETAlist[[6]] <- list(list(K1=Kmat_micro$fungi$root,model="RKHS"),
                         list(K2=Kmat_micro$fungi$rhizosphere,model="RKHS"))
    ETAlist[[7]] <- list(list(K1=Kmat_micro$bacteria$root,model="RKHS"),
                         list(K2=Kmat_micro$fungi$root,model="RKHS"))
    ETAlist[[8]] <- list(list(K1=Kmat_micro$bacteria$rhizosphere,model="RKHS"),
                         list(K2=Kmat_micro$fungi$rhizosphere,model="RKHS"))
    ETAlist[[9]] <- list(list(K1=Kmat_micro$bacteria$rhizosphere,model="RKHS"),
                         list(K2=Kmat_micro$fungi$rhizosphere,model="RKHS"),
                         list(K3=Kmat_micro$bacteria$root,model="RKHS"),
                         list(K4=Kmat_micro$fungi$root,model="RKHS"))
    ETAlist[[10]] <- list(list(K1=Kmat_micro$bacteria$root,model="RKHS"),
                          list(K2=Kmat,model="RKHS"))
    ETAlist[[11]] <- list(list(K1=Kmat_micro$bacteria$rhizosphere,model="RKHS"),
                          list(K2=Kmat,model="RKHS"))
    ETAlist[[12]] <- list(list(K1=Kmat_micro$fungi$root,model="RKHS"),
                          list(K2=Kmat,model="RKHS"))
    ETAlist[[13]] <- list(list(K1=Kmat_micro$fungi$rhizosphere,model="RKHS"),
                          list(K2=Kmat,model="RKHS"))
    ETAlist[[14]] <- list(list(K1=Kmat_micro$bacteria$root,model="RKHS"),
                          list(K2=Kmat_micro$bacteria$rhizosphere,model="RKHS"),
                          list(K3=Kmat,model="RKHS"))
    ETAlist[[15]] <- list(list(K1=Kmat_micro$fungi$root,model="RKHS"),
                          list(K2=Kmat_micro$fungi$rhizosphere,model="RKHS"),
                          list(K3=Kmat,model="RKHS"))
    ETAlist[[16]] <- list(list(K1=Kmat_micro$bacteria$root,model="RKHS"),
                          list(K2=Kmat_micro$fungi$root,model="RKHS"),
                          list(K3=Kmat,model="RKHS"))
    ETAlist[[17]] <- list(list(K1=Kmat_micro$bacteria$rhizosphere,model="RKHS"),
                          list(K2=Kmat_micro$fungi$rhizosphere,model="RKHS"),
                          list(K3=Kmat,model="RKHS"))
    ETAlist[[18]] <- list(list(K1=Kmat_micro$bacteria$rhizosphere,model="RKHS"),
                          list(K2=Kmat_micro$fungi$rhizosphere,model="RKHS"),
                          list(K3=Kmat_micro$bacteria$root,model="RKHS"),
                          list(K4=Kmat_micro$fungi$root,model="RKHS"),
                          list(K5=Kmat,model="RKHS"))
    names(ETAlist) <- scenario_name
    
    
    for (i in 1:nscenario)
    {
      ypred <- array(0,ngen)
      for (m in 1:ngen)
      {
        yNA <- y
        yNA[m] <- NA
        fit.gblup <- BGLR(yNA,
                          ETA=ETAlist[[i]],
                          nIter=5000,
                          burnIn = 2000,
                          verbose = FALSE,
                          saveAt = paste0("GP_G_and_Micro_predict_ionome_",
                                          traitname[l],"_",
                                          treatmentname[k],".txt")
        )
        ypred[m] <- fit.gblup$yHat[m]
      } # end loop m
      
      res_pa[ntreatment*(l-1)+k,scenario_name[i]] <- cor(y,ypred)
      
      nboot <- 1000
      corvec <- array(0,nboot)
      for (s in 1:nboot)
      {
        index <- sample(1:ngen,ngen,replace=TRUE)
        corvec[s] <- cor(y[index],ypred[index]) 
      }
      res_sd[ntreatment*(l-1)+k,scenario_name[i]] <- sd(corvec)
      
      cat("trait",l,"treatment",k,"scenario",i,"completed\n")
    } # end loop i
    
  } # end loop k
  
} # end loop l

write.table(res_pa,"results_G_and_microbiome_predict_Ionome_pa.txt",quote=FALSE)
write.table(res_sd,"results_G_and_microbiome_predict_Ionome_sd.txt",quote=FALSE)

