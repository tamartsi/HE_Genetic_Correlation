require(tidyr)
library(config)
args = commandArgs(trailingOnly=TRUE)
if(length(args)<1){
  stop("Config file not provided")
}else{
  config <- config::get(file = args[1])
}

#read GCTA VCs
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
for(i in 1:length(phenotypes)){ # i=1
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

#read all chunks, recalculate NormK, NormHH with the GCTA VCs
covariates0 = config$factor_covariates
covariates1 <- config$numeric_covariates
covariates = c(covariates0,covariates1)

outfileprefix <- config$outfileprefix
workdirname <- config$workdirname
ONEBATCH <- config$ONEBATCH
fstart <- length(covariates)+2
IDFIELD <- config$idfield_name


outdirname = workdirname
ofilename = paste0(outfileprefix,"_combine.txt")

print("Loading Phenotypes...")
phenotypes <- readRDS(config$phenotypes_filename)
rownames(phenotypes) <- phenotypes[[config$idfield_name]]
phenotypes <- phenotypes[,c(config$idfield_name,covariates,config$phenotypes_fieldnames)]

print("Loading relatedness matrices...")
mat_phi <- readRDS(config$kinship_matrix_filename)
mat_hh <- readRDS(config$household_matrix_filename)

cn <- colnames(mat_phi);cn <- cn[!is.na(cn)]
cn <- cn[cn %in% colnames(mat_hh)]
cn <- cn[cn %in% phenotypes[[IDFIELD]]]

rownames(phenotypes) <- phenotypes[[IDFIELD]]
pdata <- phenotypes[cn,]


phenotypes<-config$phenotypes_fieldnames
tpname <- paste0(workdirname,outfileprefix,"_Phenotypes",length(phenotypes),"_PerBatch",ONEBATCH,"_accessory.RData")
if(file.exists(tpname)){
  print("Loading accessory data")
  load(tpname)
}else{
  stop(paste0("Generating and saving accessory data: ",tpname))
}

NN <- ceiling(length(accD$tpairs$a) / ONEBATCH)

outfilename = paste0(outdirname,ofilename)

print(paste0("Identified ",length(phenotypes)," phenotypes and ", NN, " files"))

#wide format for later figure generation
res2=list()
res2$K=matrix(NA, nrow = length(phenotypes), ncol = length(phenotypes))
res2$NormK=matrix(NA, nrow = length(phenotypes), ncol = length(phenotypes))
res2$HH=matrix(NA, nrow = length(phenotypes), ncol = length(phenotypes))
res2$NormHH=matrix(NA, nrow = length(phenotypes), ncol = length(phenotypes))
res2$cSp=matrix(NA, nrow = length(phenotypes), ncol = length(phenotypes))
res2$Sp=matrix(NA, nrow = length(phenotypes), ncol = length(phenotypes))
res2$cSpPval=matrix(NA, nrow = length(phenotypes), ncol = length(phenotypes))
res2$BestpvalK=matrix(NA, nrow = length(phenotypes), ncol = length(phenotypes))
res2$BestpvalHH=matrix(NA, nrow = length(phenotypes), ncol = length(phenotypes))
res2$n=matrix(NA, nrow = length(phenotypes), ncol = length(phenotypes))

res <- NA
tcnt <- 0
for (N in 1:NN){ # N <- 1
  print(paste0("Loading batch ",N," of size ",ONEBATCH," (out of ",NN,")"))
  
  tinfname=paste0(outfileprefix,"_",NN,"batches_",ONEBATCH,"perbatch_batch_",N,".txt")
  if (length(list.files(path=workdirname,pattern=tinfname))<1){
    print("Missing batch, skipping...")
    next
  }
  tcnt <- tcnt + 1
  tdata<-read.table(paste0(workdirname,tinfname),header = TRUE,sep = '\t')
  cn<-colnames(tdata)[2:length(tdata)];rn <- as.character(tdata[,1]);
  tdata<-as.data.frame(tdata[,2:length(tdata)]);colnames(tdata)<-cn;rownames(tdata)<-rn
  tnames <- names(tdata)

  i=1  
  pheno.1 <- strsplit(cn[i],"\\.")[[1]][1]
  pheno.2 <- strsplit(cn[i],"\\.")[[1]][2]
  print(paste0(pheno.1," - ",pheno.2))

  tK <- as.numeric(tdata["GenCorr",cn[i]])
  if(is.na(tK)){
    tNormK <- NA
  }else if(tK==-100){
    tNormK <- -100
  }else{
    tNormK <- sqrt(VCg1[[pheno.1]]) * sqrt(VCg1[[pheno.2]]) * tK / sqrt((VCg1[[pheno.1]]+VCg2[[pheno.1]]+VCe[[pheno.1]])*(VCg1[[pheno.2]]+VCg2[[pheno.2]]+VCe[[pheno.2]]))
  }
  print(paste0(tNormK, " - ", tdata["NormGenCorr",cn[i]]))
  tdata["NormGenCorr",cn[i]] <- tNormK

  tE <- as.numeric(tdata["ResidCorr",cn[i]])
  if(is.na(tE)){
    tNormE <- NA
  }else if(tE==-100){
    tNormE <- -100
  }else{
    tNormE <- sqrt(VCe[[pheno.1]]) * sqrt(VCe[[pheno.2]]) * tE / sqrt((VCg1[[pheno.1]]+VCg2[[pheno.1]]+VCe[[pheno.1]])*(VCg1[[pheno.2]]+VCg2[[pheno.2]]+VCe[[pheno.2]]))
  }

  tHH <- as.numeric(tdata["HHCorr",cn[i]])
  if(is.na(tHH)){
    tNormHH <- NA
  }else if(tHH==-100){
    tNormHH <- -100
  }else{
    tNormHH <- sqrt(VCg2[[pheno.1]]) * sqrt(VCg2[[pheno.2]]) * tHH / sqrt((VCg1[[pheno.1]]+VCg2[[pheno.1]]+VCe[[pheno.1]])*(VCg1[[pheno.2]]+VCg2[[pheno.2]]+VCe[[pheno.2]]))
  }
  print(paste0(tNormHH, " - ", tdata["NormHHCorr",cn[i]]))
  tdata["NormHHCorr",cn[i]] <- tNormHH
  
  
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
    print(tcSP)
    z <- 0.5 * log((1+tcSP)/(1-tcSP))
    zse <- 1/sqrt(as.numeric(tdata["n_people",cn[i]])-3)
    tcSPpval <- min(pnorm(z, sd=zse), pnorm(z, lower.tail=F, sd=zse))*2
  }
  
  tdata["modelPheno",cn[i]] <- tcSP
  tdata["modelPheno_pval",cn[i]] <- tcSPpval
  
  #combine
  if(all(is.na(res))){
    res <- tdata
  }else{
    res <- cbind(res,tdata)
  }

  i1 <- as.numeric(match(pheno.1,phenotypes))
  i2 <- as.numeric(match(pheno.2,phenotypes))
  
  res2$K[cbind(i1,i2)] <- tdata["GenCorr",cn[i]]; res2$K[cbind(i2,i1)] <- tdata["GenCorr",cn[i]]; 
  res2$HH[cbind(i1,i2)] <- tdata["HHCorr",cn[i]]; res2$HH[cbind(i2,i1)] <- tdata["HHCorr",cn[i]]; 
  res2$NormK[cbind(i1,i2)] <- tdata["NormGenCorr",cn[i]]; res2$NormK[cbind(i2,i1)] <- tdata["NormGenCorr",cn[i]]; 
  res2$NormHH[cbind(i1,i2)] <- tdata["NormHHCorr",cn[i]]; res2$NormHH[cbind(i2,i1)] <- tdata["NormHHCorr",cn[i]]; 
  res2$n[cbind(i1,i2)] <- tdata["n_people",cn[i]]; res2$n[cbind(i2,i1)] <- tdata["n_people",cn[i]]; 
  res2$BestpvalK[cbind(i1,i2)] <- tdata["GenCorrQ_sig",cn[i]]; res2$BestpvalK[cbind(i2,i1)] <- tdata["GenCorrQ_sig",cn[i]]; 
  res2$BestpvalHH[cbind(i1,i2)] <- tdata["HHCorrQ_sig",cn[i]]; res2$BestpvalHH[cbind(i2,i1)] <- tdata["HHCorrQ_sig",cn[i]]; 
  res2$cSp[cbind(i1,i2)] <- tdata["modelPheno",cn[i]]; res2$cSp[cbind(i2,i1)] <- tdata["modelPheno",cn[i]]; 
  res2$Sp[cbind(i1,i2)] <- tdata["Spearman",cn[i]]; res2$Sp[cbind(i2,i1)] <- tdata["Spearman",cn[i]]; 
  res2$cSpPval[cbind(i1,i2)] <- tdata["modelPheno_pval",cn[i]]; res2$cSpPval[cbind(i2,i1)] <- tdata["modelPheno_pval",cn[i]]; 
  
}

if (tcnt<NN){
  print(paste0("Processed only ",tcnt," files out of ",NN,". Can't find the rest..."))
}else{
  print("Saving the resulting files...")
}  

#flip the dataframe, add columns for phnotypes, calculate FDRs
res <- as.data.frame(t(res))
rn <- rownames(res); rn <- gsub("\\.","-",rn); rownames(res) <- rn
p12 <- strsplit(rn,"-"); p12 <- do.call(rbind,p12)
res$pheno1 <- c(p12[,1])
res$pheno2 <- c(p12[,2])

#FDR adjustment
res$GenCorr_Fisher_pval_FDR <- p.adjust(res$Fisher_k_pval)
res$HHCorr_Fisher_pval_FDR <- p.adjust(res$Fisher_HH_pval)
res$GenCorr_F_pval_FDR <- p.adjust(res$GenCorr_F_pval)
res$HHCorr_F_pval_FDR <- p.adjust(res$HHCorr_F_pval)


tfname <- gsub(".txt",paste0("_longformat.txt"),outfilename)

usefields <- c("pheno1","pheno2","modelPheno","modelPheno_pval","n_people","NormGenCorr","GenCorr",
               "GenCorr_F_pval", "GenCorrSD_CI_low", "GenCorrSD_CI_high","GenCorr_F_pval_FDR","GenCorrQ_CI_low","GenCorrQ_CI_high","GenCorrQ_sig",
               "NormHHCorr","HHCorr","HHCorr_F_pval",
               "HHCorrSD_CI_low", "HHCorrSD_CI_high","HHCorr_F_pval_FDR","HHCorrQ_CI_low", "HHCorrQ_CI_high", "HHCorrQ_sig",
               "Fisher_k_pval","Fisher_k_CI_Low","Fisher_k_CI_High","GenCorr_Fisher_pval_FDR",
               "Fisher_HH_pval","Fisher_HH_CI_Low","Fisher_HH_CI_High","HHCorr_Fisher_pval_FDR",
               "ResidCorr")
res3 <- res[,usefields]

newfields <- c("pheno1","pheno2","modelPheno","modelPheno_pval","n_people","FractGenCorr","GenCorr",
               "GenCorr_F_pval", "GenCorrSD_CI_low", "GenCorrSD_CI_high","GenCorr_F_pval_FDR","GenCorrQ_CI_low","GenCorrQ_CI_high","GenCorrQ_sig",
               "FractHHCorr","HHCorr","HHCorr_F_pval",
               "HHCorrSD_CI_low", "HHCorrSD_CI_high","HHCorr_F_pval_FDR","HHCorrQ_CI_low", "HHCorrQ_CI_high", "HHCorrQ_sig",
               "GenCorr_Fisher_pval","GenCorr_Fisher_CI_low", "GenCorr_Fisher_CI_high","GenCorr_Fisher_pval_FDR",
               "HHCorr_Fisher_pval","HHCorr_Fisher_CI_low", "HHCorr_Fisher_CI_high","HHCorr_Fisher_pval_FDR",
               "ResidCorr")
colnames(res3) <- newfields

res3 <- res3[res3$pheno1!=res3$pheno2,]

write.table(res3, file = tfname, sep = "\t", row.names = F)

######################### Write separate files ###############################

tfname <- gsub(".txt",paste0("_K.txt"),outfilename)
res3<-data.frame(res2$K)
rownames(res3)<-phenotypes
colnames(res3)<-phenotypes
write.table(res3, file = tfname, sep = "\t", col.names=NA)

tfname <- gsub(".txt",paste0("_HH.txt"),outfilename)
res3<-data.frame(res2$HH)
rownames(res3)<-phenotypes
colnames(res3)<-phenotypes
write.table(res3, file = tfname, sep = "\t", col.names=NA)

tfname <- gsub(".txt",paste0("_NormK.txt"),outfilename)
res3<-data.frame(res2$NormK)
rownames(res3)<-phenotypes
colnames(res3)<-phenotypes
write.table(res3, file = tfname, sep = "\t", col.names=NA)

tfname <- gsub(".txt",paste0("_NormHH.txt"),outfilename)
res3<-data.frame(res2$NormHH)
rownames(res3)<-phenotypes
colnames(res3)<-phenotypes
write.table(res3, file = tfname, sep = "\t", col.names=NA)

tfname <- gsub(".txt",paste0("_BootstrapQ_sig_K.txt"),outfilename)
res3<-data.frame(res2$BestpvalK)
rownames(res3)<-phenotypes
colnames(res3)<-phenotypes
write.table(res3, file = tfname, sep = "\t", col.names=NA)

tfname <- gsub(".txt",paste0("_BootstrapQ_sig_HH.txt"),outfilename)
res3<-data.frame(res2$BestpvalHH)
rownames(res3)<-phenotypes
colnames(res3)<-phenotypes
write.table(res3, file = tfname, sep = "\t", col.names=NA)

tfname <- gsub(".txt",paste0("_n.txt"),outfilename)
res3<-data.frame(res2$n)
rownames(res3)<-phenotypes
colnames(res3)<-phenotypes
write.table(res3, file = tfname, sep = "\t", col.names=NA)

tfname <- gsub(".txt",paste0("_calcSpearman.txt"),outfilename)
res3<-data.frame(res2$cSp)
rownames(res3)<-phenotypes
colnames(res3)<-phenotypes
write.table(res3, file = tfname, sep = "\t", col.names=NA)

tfname <- gsub(".txt",paste0("_Spearman.txt"),outfilename)
res3<-data.frame(res2$Sp)
rownames(res3)<-phenotypes
colnames(res3)<-phenotypes
write.table(res3, file = tfname, sep = "\t", col.names=NA)

tfname <- gsub(".txt",paste0("_calcSpearmanPvalue.txt"),outfilename)
res3<-data.frame(res2$cSpPval)
rownames(res3)<-phenotypes
colnames(res3)<-phenotypes
write.table(res3, file = tfname, sep = "\t", col.names=NA)

############### Generate Heritablities files
print("Generate Heritablities files")

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

print("Done.")
