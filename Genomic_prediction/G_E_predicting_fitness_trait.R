
## to run the code please first cuctomize the following file directory

#setwd("~/Genomic_prediction/")
#datapwd <- "~/Input_data/"

library(BGLR)


###----------------------------------------------------------------------------------
#      Prediction for fitness trait using genomic data and environmental characters
###----------------------------------------------------------------------------------

Kmat <- read.table(paste0(datapwd,"VanRaden_G_matrix_marker_filtered_for_GP.txt"),header=TRUE)
Kmat <- as.matrix(Kmat)

Kmat_2 <- read.table(paste0(datapwd,"Covariance_matrix_EC_for_GP.txt"),header=TRUE)
Kmat_2 <- as.matrix(Kmat_2)

name_geno <- rownames(Kmat)
ngen <- length(name_geno)

herit <- read.table(paste0(datapwd,"Result_heritability_fitness.txt"),header=TRUE)

traitname <- rownames(herit)
ntrait <- length(traitname)

treatmentname <- colnames(herit)
ntreatment <- length(treatmentname)

#----------------------------------------------------
#    GBLUP 
#----------------------------------------------------

ETA  <- list(list(K=Kmat,model="RKHS"))

for (i in 1:ntrait)
{
  PHdata <- read.table(paste0(datapwd,"BLUE_",traitname[i],".txt"),header=TRUE)
  
  res <- matrix(0,ntreatment,3)
  rownames(res) <- treatmentname
  colnames(res) <- c("Heritability","Prediction_ability","Standard_error")
  
  res[,"Heritability"] <- t(herit[traitname[i],])
  
  for (k in 1:ntreatment)
  {
    y=PHdata[name_geno,treatmentname[k]]
    ngen <- length(y)
    ypred <- array(0,ngen)
    for (n in 1:ngen)
    {
      yNA <- y
      yNA[n] <- NA
      fit.gblup <- BGLR(yNA,
                        ETA=ETA,
                        nIter=5000,
                        burnIn = 2000,
                        verbose = FALSE,
                        saveAt = paste0("GP_GBLUP_pheno_",treatmentname[k]))
      ypred[n] <- fit.gblup$yHat[n]
    }
    
    res[k,"Prediction_ability"] <- cor(y,ypred)
    
    nboot <- 1000
    corvec <- array(0,nboot)
    for (s in 1:nboot)
    {
      index <- sample(1:ngen,ngen,replace=TRUE)
      corvec[s] <- cor(y[index],ypred[index]) 
    }
    res[k,"Standard_error"] <- sd(corvec)
    
    cat("GP for",traitname[i],treatmentname[k],"completed\n")
  }# end loop k
  
  write.table(res,paste0("result_GP_GBLUP_",traitname[i],".txt"))
  
}# end loop i
      


#----------------------------------------------------
#    EC-GBLUP 
#----------------------------------------------------

ETA_2  <- list(list(K=Kmat_2,model="RKHS"))


for (i in 1:ntrait)
{
  PHdata <- read.table(paste0(datapwd,"BLUE_",traitname[i],".txt"),header=TRUE)
  
  res_2 <- matrix(0,ntreatment,3)
  rownames(res_2) <- treatmentname
  colnames(res_2) <- c("Heritability","Prediction_ability","Standard_error")
  
  res_2[,"Heritability"] <- t(herit[traitname[i],])
  
  for (k in 1:ntreatment)
  {
    y=PHdata[name_geno,treatmentname[k]]
    ngen <- length(y)
    ypred <- array(0,ngen)
    for (n in 1:ngen)
    {
      yNA <- y
      yNA[n] <- NA
      fit.gblup <- BGLR(yNA,
                        ETA=ETA_2,
                        nIter=5000,
                        burnIn = 2000,
                        verbose = FALSE,
                        saveAt = paste0("GP_GBLUP_pheno_",treatmentname[k]))
      ypred[n] <- fit.gblup$yHat[n]
    }
    
    res_2[k,"Prediction_ability"] <- cor(y,ypred)
    
    nboot <- 1000
    corvec <- array(0,nboot)
    for (s in 1:nboot)
    {
      index <- sample(1:ngen,ngen,replace=TRUE)
      corvec[s] <- cor(y[index],ypred[index]) 
    }
    res_2[k,"Standard_error"] <- sd(corvec)
    
    cat("GP for",traitname[i],treatmentname[k],"completed\n")
  }# end loop k
  
  write.table(res_2,paste0("result_GP_ECBLUP_",traitname[i],".txt"))
  
}# end loop i



#----------------------------------------------------
#    EC-GBLUP 
#----------------------------------------------------

ETA_3  <- list(list(K1=Kmat,model="RKHS"),
               list(K2=Kmat_2,model="RKHS"))


for (i in 1:ntrait)
{
  PHdata <- read.table(paste0(datapwd,"BLUE_",traitname[i],".txt"),header=TRUE)
  
  res_3 <- matrix(0,ntreatment,3)
  rownames(res_3) <- treatmentname
  colnames(res_3) <- c("Heritability","Prediction_ability","Standard_error")
  
  res_3[,"Heritability"] <- t(herit[traitname[i],])
  
  for (k in 1:ntreatment)
  {
    y=PHdata[name_geno,treatmentname[k]]
    ngen <- length(y)
    ypred <- array(0,ngen)
    for (n in 1:ngen)
    {
      yNA <- y
      yNA[n] <- NA
      fit.gblup <- BGLR(yNA,
                        ETA=ETA_3,
                        nIter=5000,
                        burnIn = 2000,
                        verbose = FALSE,
                        saveAt = paste0("GP_GBLUP_pheno_",treatmentname[k]))
      ypred[n] <- fit.gblup$yHat[n]
    }
    
    res_3[k,"Prediction_ability"] <- cor(y,ypred)
    
    nboot <- 1000
    corvec <- array(0,nboot)
    for (s in 1:nboot)
    {
      index <- sample(1:ngen,ngen,replace=TRUE)
      corvec[s] <- cor(y[index],ypred[index]) 
    }
    res_3[k,"Standard_error"] <- sd(corvec)
    
    cat("GP for",traitname[i],treatmentname[k],"completed\n")
  }# end loop k
  
  write.table(res_3,paste0("result_GP_ECGBLUP_",traitname[i],".txt"))
  
}# end loop i


