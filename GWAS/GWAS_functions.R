library(BGLR)

GetVanRadenKinMat <- function(M){alf <- colMeans(M)/2
  W <- t(M)-2*alf
  G <- t(W)%*%W/(2*sum(alf*(1-alf)))
  return(G)
}

MakeSemiPDmatrix <- function(X){
  eigenX <- eigen(X)
  V <- eigenX$vectors
  d <- eigenX$values
  pos <- which(abs(d) < 1e-5)
  minp <- min(d[-pos])
  d[pos] <- 10^(floor(log10(minp)))
  return(V%*%diag(d)%*%t(V))
}

Fit_nullModel_BGLR <- function(y,M,G=NULL,X=NULL){
  res <- list()
  if (is.null(G))
  {
    G <- GetVanRadenKinMat(M)
  }              
  if (is.null(X))
  {
    ETA <- list(list(K=G,model="RKHS"))
  }else{
    ETA <- list(list(X=X,model="FIXED"),
                list(K=G,model="RKHS"))
  }
  fit_BGLR_null <- BGLR(y=y,
                        ETA=ETA,
                        nIter=10000,
                        burnIn=1000,
                        saveAt="fit_nullmodel_",
                        verbose=FALSE) 
  if (is.null(X))
  {
    res$vg <- fit_BGLR_null$ETA[[1]]$varU
  }else{
    res$vg <- fit_BGLR_null$ETA[[2]]$varU    
  }
  res$ve <- fit_BGLR_null$varE
  res$Gmat <- G
  return(res)         
}

Check_MarColinearWithCovar <- function(M,X){
  m <- ncol(M)
  res <- array(0,m)
  for (i in 1:m)
  {
    res[i] <- ifelse(qr(X)$rank < qr(cbind(X,M[,i]))$rank,1,0)
  }
  return(res)
}

FastGWAS_P3D_basic <- function(y,M,G=NULL,varg=NULL,vare=NULL,wt=NULL){ 
  n <- nrow(M)
  m <- ncol(M)
  GWAS_res <- matrix(0,m,3)
  colnames(GWAS_res) <- c("Estimate_Eff","Chisq_Statistic","P_value")
  rownames(GWAS_res) <- colnames(M)
  if (is.null(G))
  {
    G <- GetVanRadenKinMat(M)
  } 
  if (is.null(varg)|is.null(vare))
  {
    fit_null <- Fit_nullModel_BGLR(y,M,G)
    varg <- fit_null$vg
    vare <- fit_null$ve
  } 
  lambda <- varg/vare
  I <- diag(1,n)
  if (is.null(wt))
  {  
    D <- I
  }else{
    D <- diag(wt)
  }
  K <- G*lambda+D   
  L <- chol(K)                      
  Ltinv <- t(backsolve(L,I)) 
  dta <- Ltinv%*%cbind(y,rep(1,n),M)
  z <- dta[,1]
  u <- dta[,2]
  W <- dta[,-1:-2]
  upu <- c(t(u)%*%u)
  upy <- c(t(u)%*%z)
  wpw <- apply(W^2,2,sum)
  upwsq <- apply(W*as.vector(u),2,sum)^2
  f <- upu*wpw-upwsq
  beta <- (t(W)%*%(z*upu-u*upy))/f
  varb <- upu*vare/f
  chisqstat <- beta^2/varb
  GWAS_res[,1] <- beta
  GWAS_res[,2] <- chisqstat
  GWAS_res[,3] <- pchisq(chisqstat,df=1,lower.tail = FALSE) 
  return(GWAS_res)
}


FastGWAS_P3D_basic_v2 <- function(y,M,G=NULL,varg=NULL,vare=NULL,wt=NULL,X=NULL){                     
  n <- nrow(M)
  m <- ncol(M)
  GWAS_res <- matrix(0,m,3)
  colnames(GWAS_res) <- c("Estimate_Eff","t_Statistic","P_value")
  rownames(GWAS_res) <- colnames(M)
  
  if (is.null(G))
  {
    G <- GetVanRadenKinMat(M)
  } 
  if (is.null(varg)|is.null(vare))
  {
    fit_null <- Fit_nullModel_BGLR(y,M,G,X)
    varg <- fit_null$vg
    vare <- fit_null$ve
  }
  
  if (is.null(X))
  {
    Xtilde <- matrix(rep(1,n),n,1)
  }else{
    Xtilde <- cbind(rep(1,n),X)
  }
  ckvec <- Check_MarColinearWithCovar(M,Xtilde)
  colrpos <- which(ckvec==0)     
  
  lambda <- varg/vare
  I <- diag(1,n)
  if (is.null(wt))
  {  
    D <- I
  }else{
    D <- diag(wt)
  }
  K <- G*lambda+D   
  L <- chol(K)                      
  Ltinv <- t(backsolve(L,I)) 
  sdta <- Ltinv%*%cbind(y,Xtilde)
  z <- sdta[,1]
  U <- sdta[,-1,drop=FALSE]
  UpU <- t(U)%*%U
  UpUinv <- solve(UpU)
  
  betavec <- array(0,m)
  tstatvec <- array(0,m)
  pvalvec <- array(0,m)
  
  if (length(colrpos) > 0)
  {
    indexset <- setdiff(1:m,colrpos)
    betavec[colrpos] <- NA
    tstatvec[corlpos] <- NA
    pvalvec[colrpos] <- NA
  }else{
    indexset <- 1:m
  }
  
  for (l in indexset)
  {
    w <- Ltinv%*%M[,l,drop=FALSE]
    wpU <- t(w)%*%U
    Upy <- t(U)%*%z
    wpy <- c(t(w)%*%z)
    wpw <- c(t(w)%*%w)
    f <- wpw-wpU%*%UpUinv%*%t(wpU)
    beta <- (wpy-wpU%*%UpUinv%*%Upy)/f
    alpha <- solve(UpU-t(wpU)%*%wpU/wpw)%*%(Upy-t(wpU)*wpy/wpw)
    resid <- z-U%*%alpha-w%*%beta
    s2 <- sum(resid^2)/(n-1-ncol(Xtilde))
    tstat <- beta/sqrt(s2/f)
    betavec[l] <- beta
    tstatvec[l] <- tstat
    pvalvec[l] <- 2*pt(abs(tstat),df=n-1-ncol(Xtilde),lower.tail = FALSE) 
  }
  
  GWAS_res[,1] <- betavec
  GWAS_res[,2] <- tstatvec
  GWAS_res[,3] <- pvalvec
  
  return(GWAS_res)
}