## to run the code please first cuctomize the following file directory

#setwd("~/Genomic_prediction/")
#datapwd <- "~/Input_data/"

library(BGLR)


###----------------------------------------------------------------------------------
#      Prediction for microbiome families using genomic data and environmental characters
###----------------------------------------------------------------------------------


Kmat <- read.table(paste0(datapwd,"VanRaden_G_matrix_marker_filtered_for_GP.txt"),header=TRUE)
Kmat <- as.matrix(Kmat)

Kmat_2 <- read.table(paste0(datapwd,"Covariance_matrix_EC_for_GP.txt"),header=TRUE)
Kmat_2 <- as.matrix(Kmat_2)

name_geno <- rownames(Kmat)
ngen <- length(name_geno)

microbiometype <- c("fungi","bacteria")
nmicro <- length(microbiometype)

sampletype <- c("root","rhizosphere")
nsample <- length(sampletype)

treatmentname <- c("CK","D","LN","LP")
ntreatment <- length(treatmentname)

#----------------------------------------------------
#    GBLUP 
#----------------------------------------------------

ETA  <- list(list(K=Kmat,model="RKHS"))

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
      if (length(pos)==0)
      {
        break
      }
      asvname_filtered <- rownames(herit)[pos]
      nasv <- length(asvname_filtered)
      
      res <- matrix(0,nasv,3)
      rownames(res) <- asvname_filtered
      colnames(res) <- c("Heritability","Prediction_ability","Standard_error")
      
      res[,"Heritability"] <- herit[asvname_filtered,treatmentname[k]]
      
      for (l in 1:nasv)
      {
        y=PHdata[name_geno,asvname_filtered[l]]
        ngen <- length(y)
        ypred <- array(0,ngen)
        for (m in 1:ngen)
        {
          yNA <- y
          yNA[m] <- NA
          fit.gblup <- BGLR(yNA,
                            ETA=ETA,
                            nIter=5000,
                            burnIn = 2000,
                            verbose = FALSE,
                            saveAt = paste0("GP_GBLUP_",microbiometype[i],"_family_",sampletype[j],"_",treatmentname[k]))
          ypred[m] <- fit.gblup$yHat[m]
        }
        
        res[l,"Prediction_ability"] <- cor(y,ypred)
        
        nboot <- 1000
        corvec <- array(0,nboot)
        for (s in 1:nboot)
        {
          index <- sample(1:ngen,ngen,replace=TRUE)
          corvec[s] <- cor(y[index],ypred[index]) 
        }
        res[l,"Standard_error"] <- sd(corvec)
        
        cat("GP for ASV",l,"completed\n")
      }# end loop l
      
      write.table(res,paste0("result_GP_GBLUP_",microbiometype[i],"_family_",sampletype[j],"_",treatmentname[k],".txt"))
      
    }# end loop k
    
  }# end loop j
  
}# end loop i
  


#----------------------------------------------------
#    EC-GBLUP 
#----------------------------------------------------

ETA_2  <- list(list(K=Kmat_2,model="RKHS"))


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
      if (length(pos)==0)
      {
        break
      }
      asvname_filtered <- rownames(herit)[pos]
      nasv <- length(asvname_filtered)
      
      res_2 <- matrix(0,nasv,3)
      rownames(res_2) <- asvname_filtered
      colnames(res_2) <- c("Heritability","Prediction_ability","Standard_error")
      
      res_2[,"Heritability"] <- herit[asvname_filtered,treatmentname[k]]
      
      for (l in 1:nasv)
      {
        y=PHdata[name_geno,asvname_filtered[l]]
        ngen <- length(y)
        ypred <- array(0,ngen)
        for (m in 1:ngen)
        {
          yNA <- y
          yNA[m] <- NA
          fit.gblup <- BGLR(yNA,
                            ETA=ETA_2,
                            nIter=5000,
                            burnIn = 2000,
                            verbose = FALSE,
                            saveAt = paste0("GP_ECBLUP_",microbiometype[i],"_family_",sampletype[j],"_",treatmentname[k]))
          ypred[m] <- fit.gblup$yHat[m]
        }
        
        res_2[l,"Prediction_ability"] <- cor(y,ypred)
        
        nboot <- 1000
        corvec <- array(0,nboot)
        for (s in 1:nboot)
        {
          index <- sample(1:ngen,ngen,replace=TRUE)
          corvec[s] <- cor(y[index],ypred[index]) 
        }
        res_2[l,"Standard_error"] <- sd(corvec)
        
        cat("GP for ASV",l,"completed\n")
      }# end loop l
      
      write.table(res_2,paste0("result_GP_ECBLUP_",microbiometype[i],"_family_",sampletype[j],"_",treatmentname[k],".txt"))
      
    }# end loop k
    
  }# end loop j
  
}# end loop i


#----------------------------------------------------
#    EC-GBLUP 
#----------------------------------------------------

ETA_3  <- list(list(K1=Kmat,model="RKHS"),
               list(K2=Kmat_2,model="RKHS"))


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
      if (length(pos)==0)
      {
        break
      }
      asvname_filtered <- rownames(herit)[pos]
      nasv <- length(asvname_filtered)
      
      res_3 <- matrix(0,nasv,3)
      rownames(res_3) <- asvname_filtered
      colnames(res_3) <- c("Heritability","Prediction_ability","Standard_error")
      
      res_3[,"Heritability"] <- herit[asvname_filtered,treatmentname[k]]
      
      for (l in 1:nasv)
      {
        y=PHdata[name_geno,asvname_filtered[l]]
        ngen <- length(y)
        ypred <- array(0,ngen)
        for (m in 1:ngen)
        {
          yNA <- y
          yNA[m] <- NA
          fit.gblup <- BGLR(yNA,
                            ETA=ETA_3,
                            nIter=5000,
                            burnIn = 2000,
                            verbose = FALSE,
                            saveAt = paste0("GP_ECGBLUP_",microbiometype[i],"_family_",sampletype[j],"_",treatmentname[k]))
          ypred[m] <- fit.gblup$yHat[m]
        }
        
        res_3[l,"Prediction_ability"] <- cor(y,ypred)
        
        nboot <- 1000
        corvec <- array(0,nboot)
        for (s in 1:nboot)
        {
          index <- sample(1:ngen,ngen,replace=TRUE)
          corvec[s] <- cor(y[index],ypred[index]) 
        }
        res_3[l,"Standard_error"] <- sd(corvec)
        
        cat("GP for ASV",l,"completed\n")
      }# end loop l
      
      write.table(res_3,paste0("result_GP_ECGBLUP_",microbiometype[i],"_family_",sampletype[j],"_",treatmentname[k],".txt"))
      
    }# end loop k
    
  }# end loop j
  
}# end loop i
