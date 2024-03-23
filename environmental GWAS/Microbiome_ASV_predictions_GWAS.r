
# This code is modified based on Gates et al., 2019 (10.1101/706739)
# Landrace genotype files were obtained from https://github.com/danjgates/AdaptationScripts2

library(data.table)
library(dplyr)
library(qqman)

# Read predicted root matrix: ASV7 and ASV37 results were extracted from RandomForest output: FOAMpredictions_ASVsSUBSET.csv 
# Genotypes that match the genotype file were kept
microMat <- read.csv('data/Microbiome_mexi_ASV_7_37_RFprediction.csv', header = TRUE)
colnames(microMat)[1] <- 'SampleID'
names(microMat)

#specify the layers to run
layers<- c("ASV37_bacteria_rhizosphere_CK_pred", "ASV37_bacteria_rhizosphere_LN_pred",
           "ASV37_bacteria_rhizosphere_LP_pred", "ASV37_bacteria_root_CK_pred", "ASV37_bacteria_root_LN_pred",       
           "ASV37_bacteria_root_LP_pred", "ASV7_bacteria_rhizosphere_CK_pred",  "ASV7_bacteria_rhizosphere_LN_pred", 
           "ASV7_bacteria_root_D_pred",   "ASV7_bacteria_root_LN_pred")

#make a loop to run through and make the manhattans:
sapply(layers,function(layer){
    precip<-microMat[,layer]
    names(precip)<-microMat$SampleID

    #pop structure:
    load('Genotypes/HighFilteredReduced.Rimage')
    inters<-intersect(names(precip),rownames(genoMDS))
    genoMDS<-genoMDS[inters,]
    precip<-precip[inters]
    pops<-cmdscale(dist(genoMDS),k=5)
    
    colms<-lapply(1:10,function(chr){
        mat<-data.frame(fread(paste('Genotypes/Ch',chr,'Merged.hmp.txt',sep="")))
        nms<-sapply(mat[,1],function(x) strsplit(x,'.M')[[1]][1])
        dups<-duplicated(nms)
        xx<-mat[-which(dups==TRUE),]
        rownames(xx)<-sapply(xx[,1],function(x) strsplit(x,'.M')[[1]][1])
        mat<-xx
        mat<-mat[,-1] 
        mat<-mat[names(precip),]
        mono<-apply(mat,MARGIN=2,function(x) length(table(x)))
        mat<-mat[,-which(mono==1)]
        
        #read in SNPs:
        colm<-sapply(1:ncol(mat),function(x){
            mergetab<-data.frame(cbind(geno=mat[,x],Precip=precip,V1=pops[,1],V2=pops[,2],V3=pops[,3],V4=pops[,4],V5=pops[,5]))
            #extract p-values for genotypes
            pv<-summary(lm(Precip~geno+V1+V2+V2+V4+V5,data=mergetab))$coefficient[2,4]
            #extract t-values for genotypes for z-score calculation
            tv<-summary(lm(Precip~geno+V1+V2+V2+V4+V5,data=mergetab))$coefficient[2,3]
            pred<-predict(lm(Precip~geno+V1+V2+V2+V4+V5,data=mergetab),data.frame(geno=1,V1=0,V2=0,V3=0,V4=0,V5=0))-predict(lm(Precip~geno+V1+V2+V3+V4+V5,data=mergetab),data.frame(geno=0,V1=0,V2=0,V3=0,V4=0,V5=0))
            return(c(CHR=chr,BP=as.numeric(strsplit(colnames(mat)[x],'_')[[1]][2]),P=pv,prediction=pred,z_score = tv))
        })
        return(colm)
    })
    colm<-data.frame(rbind(t(colms[[1]]),t(colms[[2]]),t(colms[[3]]),t(colms[[4]]),t(colms[[5]]),t(colms[[6]]),t(colms[[7]]),t(colms[[8]]),t(colms[[9]]),t(colms[[10]])))

    write.csv(colm, file = paste0(layer,'_mexi_results.csv'))
    png(filename = paste(layer,'mexi_result.png',sep = "_"))
    manhattan(colm, main = paste(layer,'mexi',sep = '_'), chr="CHR", bp="BP", snp="SNP", p="P",suggestiveline = F, genomewideline = F)
    dev.off()
})


#------prepare a genotype file for downstream cimmyt GWAS analysis--#
micro<-microMat[,'SampleID'] 
names(micro)<-microMat$SampleID

# read genotypic file;genotype for SRN (remove na)
colms<-lapply(1:10,function(chr){
      mat<-data.frame(fread(paste('Genotypes/Ch',chr,'Merged.hmp.txt',sep="")))
      nms<-sapply(mat[,1],function(x) strsplit(x,'.M')[[1]][1])
      dups<-duplicated(nms)
      xx<-mat[-which(dups==TRUE),]
      rownames(xx)<-sapply(xx[,1],function(x) strsplit(x,'.M')[[1]][1])
      mat<-xx
      mat<-mat[,-1] 
      mat<-mat[names(micro),]
      mono<-apply(mat,MARGIN=2,function(x) length(table(x)))
      mat<-mat[,-which(mono==1)]
      mat<-t(as.matrix(mat))
      mat<-mat*2
      return(mat)
})

#save(colms, file = 'cimmyt_genotype_list.Rimage')
#load('cimmyt_genotype_list.Rimage')

#unlist file

colm<-data.frame(rbind(colms[[1]],colms[[2]],colms[[3]],colms[[4]],colms[[5]],colms[[6]],colms[[7]],colms[[8]],colms[[9]],colms[[10]]))

save(colm, file = 'data/cimmyt_genotype_file.Rimage')