library(gtools)
library(plotrix)
library(config)
args = commandArgs(trailingOnly=TRUE)
if(length(args)<1){
  stop("Config file not provided")
  # config <- config::get(file = "/data/linkage/HCHS_SOL/Projects/2021_gencor_cog_sleep/Code/runGenCor_ALL_Base_v1.config")
}else{
  config <- config::get(file = args[1])
}

#model 2
workdirname <- config$workdirname
pheno_names <- paste0(config$phenotypes_GCTA_prefix,".pheno.names")
outfileprefix <- config$outfileprefix
ofilename = paste0(outfileprefix,"_combine.txt")
VCprefix = config$phenotypes_GCTA_prefix


outfilename = paste0(workdirname,ofilename)
phenotypes <- read.table(pheno_names,sep="\t",header = F, row.names = 1)
phenotypes <- phenotypes[[1]]

VCg1 <- rep(NA,length(phenotypes));names(VCg1) <- phenotypes
VCg2 <- rep(NA,length(phenotypes));names(VCg2) <- phenotypes
VCe <- rep(NA,length(phenotypes));names(VCe) <- phenotypes

res <- list()
for(i in 1:length(phenotypes)){ # i=6
  f=paste0(VCprefix,"_",i,".hsq")
  pheno.1 <- phenotypes[i]
  print(paste0(f," - ", pheno.1))
  if(!file.exists(f)){
    print("Missing file, skipping...")
    next
  }
  data = readLines(f)
  rsG1 <- data[which(grepl("V\\(G1)\t", data))]
  if(length(rsG1)!=1){
    print("Something is wrong")
    break
  }
  rsG2 <- data[which(grepl("V\\(G2)\t", data))]
  if(length(rsG2)!=1){
    print("Something is wrong")
    break
  }
  rsE <- data[which(grepl("V\\(e)\t", data))]
  if(length(rsE)!=1){
    print("Something is wrong")
    break
  }
  rsN <- data[which(grepl("n\t", data))]
  if(length(rsN)!=1){
    print("Something is wrong")
    break
  }
  mnG1 <- as.numeric(strsplit(rsG1,"\t",fixed = T)[[1]][2])
  seG1 <- as.numeric(strsplit(rsG1,"\t",fixed = T)[[1]][3])
  mnG2 <- as.numeric(strsplit(rsG2,"\t",fixed = T)[[1]][2])
  seG2 <- as.numeric(strsplit(rsG2,"\t",fixed = T)[[1]][3])
  mnE <- as.numeric(strsplit(rsE,"\t",fixed = T)[[1]][2])
  seE <- as.numeric(strsplit(rsE,"\t",fixed = T)[[1]][3])
  nn <- as.numeric(strsplit(rsN,"\t",fixed = T)[[1]][2])

  VCg1[[pheno.1]] <- mnG1 / (mnG1 + mnG2 + mnE)
  VCg2[[pheno.1]] <- mnG2 / (mnG1 + mnG2 + mnE)
  VCe[[pheno.1]] <- mnE / (mnG1 + mnG2 + mnE)
}

print("Recalculating NormK based on GCTA Sigmas...")
res$NormK <- as.data.frame(matrix(NA,nrow = length(phenotypes), ncol = length(phenotypes)))
colnames(res$NormK) <- phenotypes; rownames(res$NormK) <- phenotypes;
res$NormHH <- as.data.frame(matrix(NA,nrow = length(phenotypes), ncol = length(phenotypes)))
colnames(res$NormHH) <- phenotypes; rownames(res$NormHH) <- phenotypes;
res$cSP <- as.data.frame(matrix(NA,nrow = length(phenotypes), ncol = length(phenotypes)))
colnames(res$cSP) <- phenotypes; rownames(res$cSP) <- phenotypes;
res$cSPpval <- as.data.frame(matrix(NA,nrow = length(phenotypes), ncol = length(phenotypes)))
colnames(res$cSPpval) <- phenotypes; rownames(res$cSPpval) <- phenotypes;

fK <- gsub(".txt",paste0("_K.txt"),outfilename)
dataK<-read.table(fK,header = TRUE,sep = '\t')
cn<-colnames(dataK)[2:length(dataK)];rn <- as.character(dataK[,1]);dataK<-dataK[,2:length(dataK)]
colnames(dataK)<-cn;rownames(dataK)<-rn

fHH <- gsub(".txt",paste0("_HH.txt"),outfilename)
dataHH<-read.table(fHH,header = TRUE,sep = '\t')
cn<-colnames(dataHH)[2:length(dataHH)];rn <- as.character(dataHH[,1]);dataHH<-dataHH[,2:length(dataHH)]
colnames(dataHH)<-cn;rownames(dataHH)<-rn

fE <- gsub(".txt",paste0("_E.txt"),outfilename)
dataE<-read.table(fE,header = TRUE,sep = '\t')
cn<-colnames(dataE)[2:length(dataE)];rn <- as.character(dataE[,1]);dataE<-dataE[,2:length(dataE)]
colnames(dataE)<-cn;rownames(dataE)<-rn

fSp <- gsub(".txt",paste0("_Spearman.txt"),outfilename)
dataSp<-read.table(fSp,header = TRUE,sep = '\t')
cn<-colnames(dataSp)[2:length(dataSp)];rn <- as.character(dataSp[,1]);dataSp<-dataSp[,2:length(dataSp)]
colnames(dataSp)<-cn;rownames(dataSp)<-rn

fN <- gsub(".txt",paste0("_n.txt"),outfilename)
dataN<-read.table(fN,header = TRUE,sep = '\t')
cn<-colnames(dataN)[2:length(dataN)];rn <- as.character(dataN[,1]);dataN<-dataN[,2:length(dataN)]
colnames(dataN)<-cn;rownames(dataN)<-rn

for(i in 1:length(phenotypes)){ # i=1
  for(j in 1:length(phenotypes)){ # j=6
    pheno.1 <- phenotypes[i]
    pheno.2 <- phenotypes[j]

    tK <- dataK[pheno.1,pheno.2]
    if(is.na(tK)){
      tNormK <- NA
    }else if(tK==-100){
        tNormK <- -100
    }else{
      tNormK <- sqrt(VCg1[[pheno.1]]) * sqrt(VCg1[[pheno.2]]) * tK / sqrt((VCg1[[pheno.1]]+VCg2[[pheno.1]]+VCe[[pheno.1]])*(VCg1[[pheno.2]]+VCg2[[pheno.2]]+VCe[[pheno.2]]))
    }
    res$NormK[pheno.1,pheno.2] <- tNormK
    tE <- dataE[pheno.1,pheno.2]
    if(is.na(tE)){
      tNormE <- NA
    }else if(tE==-100){
      tNormE <- -100
    }else{
      tNormE <- sqrt(VCe[[pheno.1]]) * sqrt(VCe[[pheno.2]]) * tE / sqrt((VCg1[[pheno.1]]+VCg2[[pheno.1]]+VCe[[pheno.1]])*(VCg1[[pheno.2]]+VCg2[[pheno.2]]+VCe[[pheno.2]]))
    }
    tHH <- dataHH[pheno.1,pheno.2]
    if(is.na(tHH)){
      tNormHH <- NA
    }else if(tHH==-100){
      tNormHH <- -100
    }else{
      tNormHH <- sqrt(VCg2[[pheno.1]]) * sqrt(VCg2[[pheno.2]]) * tHH / sqrt((VCg1[[pheno.1]]+VCg2[[pheno.1]]+VCe[[pheno.1]])*(VCg1[[pheno.2]]+VCg2[[pheno.2]]+VCe[[pheno.2]]))
    }
    res$NormHH[i,j] <- tNormHH
    tcSP <- tNormE + tNormK + tNormHH;
    if(!is.na(tcSP)){
      if(tcSP<= -100){
        tcSP <- -100
      }
    }
      
    if(is.na(tcSP)){
      tcSPpval <- NA
    }else if (tcSP<= -100){
      tcSPpval <- -100
    }else{
      z <- 0.5 * log((1+tcSP)/(1-tcSP))
      zse <- 1/sqrt(dataN[pheno.1,pheno.2]-3)
      tcSPpval <- min(pnorm(z, sd=zse), pnorm(z, lower.tail=F, sd=zse))*2
    }
    res$cSP[pheno.1,pheno.2] <- tcSP
    res$cSPpval[pheno.1,pheno.2] <- tcSPpval
    
  }  
}

tfname <- gsub(".txt",paste0("_GCTA_VCs_NormK.txt"),outfilename)
res2<-data.frame(res$NormK)
rownames(res2)<-phenotypes
colnames(res2)<-phenotypes
write.table(res2, file = tfname, sep = "\t", col.names=NA)

tfname <- gsub(".txt",paste0("_GCTA_VCs_NormHH.txt"),outfilename)
res2<-data.frame(res$NormHH)
rownames(res2)<-phenotypes
colnames(res2)<-phenotypes
write.table(res2, file = tfname, sep = "\t", col.names=NA)

tfname <- gsub(".txt",paste0("_GCTA_VCs_calcSpearman.txt"),outfilename)
res2<-data.frame(res$cSP)
rownames(res2)<-phenotypes
colnames(res2)<-phenotypes
write.table(res2, file = tfname, sep = "\t", col.names=NA)

tfname <- gsub(".txt",paste0("_GCTA_VCs_calcSpearmanPvalue.txt"),outfilename)
res2<-data.frame(res$cSPpval)
rownames(res2)<-phenotypes
colnames(res2)<-phenotypes
write.table(res2, file = tfname, sep = "\t", col.names=NA)


############### Generate Heritablities files
ttitle <- c("K.txt","HH.txt")
VCprefix = config$phenotypes_GCTA_prefix
workdirname <- config$workdirname
outfileprefix <- config$outfileprefix
ofilename = paste0(outfileprefix,"_combine.txt")
outfilename = paste0(workdirname,ofilename)

fNormK <- gsub(".txt",paste0("_GCTA_VCs_NormK.txt"),outfilename)

for (i in 1:length(ttitle)){ # i=1
  
  #read GCTA heritabilities
  workdirname <- strsplit(config$phenotypes_GCTA_prefix,"/",fixed = T)[[1]]
  workdirname <- paste0(paste0(c(workdirname[1:(length(workdirname)-1)]),collapse = "/"),"/")
  ddata <- data.frame(
    name=factor(config$phenotypes_fieldnames,levels = rev(config$phenotypes_fieldnames)),
    value=rep(NA,length(config$phenotypes_fieldnames)),
    sd=rep(NA,length(config$phenotypes_fieldnames))
  )
  rownames(ddata) <- ddata$name
  phenotypes <- read.table(paste0(config$phenotypes_GCTA_prefix,".pheno.names"),sep="\t",header = F, row.names = 1)
  phenotypes <- phenotypes[[1]]
  
  for( j in 1:length(phenotypes)){ #j=1
    fher=paste0(VCprefix,"_",j,".hsq")
    if(!file.exists(fher)){
      print("Missing file, skipping...")
      next
    }
    data = readLines(fher)
    
    rsG1 <- data[which(grepl("V\\(G1)/Vp\t", data))]
    if(length(rsG1)!=1){
      print("Something is wrong")
      break
    }
    rsG2 <- data[which(grepl("V\\(G2)/Vp\t", data))]
    if(length(rsG2)!=1){
      print("Something is wrong")
      break
    }
    mnG1 <- as.numeric(strsplit(rsG1,"\t",fixed = T)[[1]][2])
    seG1 <- as.numeric(strsplit(rsG1,"\t",fixed = T)[[1]][3])
    mnG2 <- as.numeric(strsplit(rsG2,"\t",fixed = T)[[1]][2])
    seG2 <- as.numeric(strsplit(rsG2,"\t",fixed = T)[[1]][3])
    if(ttitle[i]=="K.txt"){
      ddata[phenotypes[j],"value"] <- mnG1
      ddata[phenotypes[j],"sd"] <- seG1
    }else if(ttitle[i]=="HH.txt"){
      ddata[phenotypes[j],"value"] <- mnG2
      ddata[phenotypes[j],"sd"] <- seG2
    }else{
      stop("Unknown parameter")
    }
  }
  
  ddata$value[ddata$value < 0] <- 0
  ddata$value[ddata$value > 1] <- 1
  ddata$value[is.na(ddata$value)] <- 0
  ddata$sd[is.na(ddata$sd)] <- 0
  
  if(ttitle[i]=="K.txt"){
    write.table(ddata,gsub("NormK.txt","K_Heritabilities.txt",fNormK),row.names = F, sep = "\t")
  }else{
    write.table(ddata,gsub("NormK.txt","HH_Heritabilities.txt",fNormK),row.names = F, sep = "\t")
  }

}

print("Done")

