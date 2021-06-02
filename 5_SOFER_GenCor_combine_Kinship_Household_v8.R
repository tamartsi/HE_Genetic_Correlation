#require(dplyr)
require(tidyr)
library(config)
args = commandArgs(trailingOnly=TRUE)
if(length(args)<1){
  stop("Config file not provided")
  # config <- config::get(file = "/data/linkage/HCHS_SOL/Projects/2021_gencor_cog_sleep/Code/runGenCor_ALL_Base_v1.config")
}else{
  config <- config::get(file = args[1])
}

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

tempN = 0
res=list()
res$K=matrix(-100, nrow = length(phenotypes), ncol = length(phenotypes))
res$NormK=matrix(-100, nrow = length(phenotypes), ncol = length(phenotypes))
res$HH=matrix(-100, nrow = length(phenotypes), ncol = length(phenotypes))
res$NormHH=matrix(-100, nrow = length(phenotypes), ncol = length(phenotypes))
res$Fisher_K_pval=matrix(-100, nrow = length(phenotypes), ncol = length(phenotypes))
res$Fisher_HH_pval=matrix(-100, nrow = length(phenotypes), ncol = length(phenotypes))
res$Bootstrap_K_pval=matrix(-100, nrow = length(phenotypes), ncol = length(phenotypes))
res$Bootstrap_HH_pval=matrix(-100, nrow = length(phenotypes), ncol = length(phenotypes))
res$KSig1=matrix(-100, nrow = length(phenotypes), ncol = length(phenotypes))
res$KSig2=matrix(-100, nrow = length(phenotypes), ncol = length(phenotypes))
res$ESig1=matrix(-100, nrow = length(phenotypes), ncol = length(phenotypes))
res$ESig2=matrix(-100, nrow = length(phenotypes), ncol = length(phenotypes))
res$HHSig1=matrix(-100, nrow = length(phenotypes), ncol = length(phenotypes))
res$HHSig2=matrix(-100, nrow = length(phenotypes), ncol = length(phenotypes))
res$E=matrix(-100, nrow = length(phenotypes), ncol = length(phenotypes))
res$cSp=matrix(-100, nrow = length(phenotypes), ncol = length(phenotypes))
res$Sp=matrix(-100, nrow = length(phenotypes), ncol = length(phenotypes))
res$SpPval=matrix(-100, nrow = length(phenotypes), ncol = length(phenotypes))
res$BestpvalK=matrix(-100, nrow = length(phenotypes), ncol = length(phenotypes))
res$BestpvalHH=matrix(-100, nrow = length(phenotypes), ncol = length(phenotypes))
res$n=matrix(-100, nrow = length(phenotypes), ncol = length(phenotypes))
res$Bootstrap_K_sd=matrix(-100, nrow = length(phenotypes), ncol = length(phenotypes))
res$Bootstrap_HH_sd=matrix(-100, nrow = length(phenotypes), ncol = length(phenotypes))
for (N in (tempN+1):NN){ # N <- (tempN+1) +1
  print(paste0("Loading batch ",N," of size ",ONEBATCH," (out of ",NN,")"))
  
  tinfname=paste0(outfileprefix,"_",NN,"batches_",ONEBATCH,"perbatch_batch_",N,".txt")
  if (length(list.files(path=workdirname,pattern=tinfname))<1){
    print("Missing batch, skipping...")
    next
  }
  tdata<-read.table(paste0(workdirname,tinfname),header = TRUE,sep = '\t')
  cn<-colnames(tdata)[2:length(tdata)];rn <- as.character(tdata[,1]);tdata<-tdata[,2:length(tdata)];colnames(tdata)<-cn;rownames(tdata)<-rn
  tnames <- names(tdata)
  # print("Processing")
  t2 <- as.data.frame(t(tdata))
  t2$K <- as.numeric(as.character(t2$K));
  t2$E <- as.numeric(as.character(t2$E));
  t2$HH <- as.numeric(as.character(t2$HH));
  t2$Fisher_k_pval <- as.numeric(as.character(t2$Fisher_k_pval));
  t2$Fisher_k_pval_adj <- t2$Fisher_k_pval
  t2$Fisher_k_pval_adj[t2$`Fisher_k_CI_Low:Fisher_k_CI_High`=="NA:NA"] <- NA
  t2$Fisher_HH_pval <- as.numeric(as.character(t2$Fisher_HH_pval));
  t2$Fisher_HH_pval_adj <- t2$Fisher_HH_pval
  t2$Fisher_HH_pval_adj[t2$`Fisher_HH_CI_Low:Fisher_HH_CI_High`=="NA:NA"] <- NA
  t2$Bootstrap_k_pval <- as.numeric(as.character(t2$Bootstrap_k_pval));
  t2$Bootstrap_MeanK <- as.numeric(as.character(t2$Bootstrap_MeanK));
  t2$Bootstrap_SdK <- as.numeric(as.character(t2$Bootstrap_SdK));
  t2$Bootstrap_HH_pval <- as.numeric(as.character(t2$Bootstrap_HH_pval));
  t2$Bootstrap_MeanHH <- as.numeric(as.character(t2$Bootstrap_MeanHH));
  t2$Bootstrap_SdHH <- as.numeric(as.character(t2$Bootstrap_SdHH));
  #Some sanity checks
  t2$Bootstrap_k_pval_adj <- t2$Bootstrap_k_pval; t2$Bootstrap_k_pval_adj[t2$Bootstrap_k_pval_adj == -100] <- NA
  t2$Bootstrap_HH_pval_adj <- t2$Bootstrap_HH_pval; t2$Bootstrap_HH_pval_adj[t2$Bootstrap_HH_pval_adj == -100] <- NA
  t2$Bootstrap_k_pval_adj[(t2$Bootstrap_k_pval<0.05) & (abs(t2$Bootstrap_MeanK)*0.8 < t2$Bootstrap_SdK)] <- NA       # std not too large
  t2$Bootstrap_k_pval_adj[(t2$Bootstrap_k_pval<0.05) & (abs(t2$K)*0.8 > abs(t2$Bootstrap_MeanK))] <- NA    # bootstrap mean not too different from K and same direction
  t2$BestpvalK <- ifelse(is.na(t2$Bootstrap_k_pval_adj), t2$Fisher_k_pval_adj, t2$Bootstrap_k_pval_adj)    # Fill blanks of Bootstrap with Fisher
  t2$Bootstrap_HH_pval_adj[(t2$Bootstrap_HH_pval<0.05) & (abs(t2$Bootstrap_MeanHH)*0.8 < t2$Bootstrap_SdHH)] <- NA       # std not too large
  t2$Bootstrap_HH_pval_adj[(t2$Bootstrap_HH_pval<0.05) & (abs(t2$HH)*0.8 > abs(t2$Bootstrap_MeanHH))] <- NA    # bootstrap mean not too different from K and same direction
  t2$BestpvalHH <- ifelse(is.na(t2$Bootstrap_HH_pval_adj), t2$Fisher_HH_pval_adj, t2$Bootstrap_HH_pval_adj)    # Fill blanks of Bootstrap with Fisher
  
  t2$NormK <- as.numeric(as.character(t2$NormK)); 
  t2$NormHH <- as.numeric(as.character(t2$NormHH)); 
  t2 <- t2 %>% separate(`ErrorSigma1:ErrorSigma2`, c("Ve1", "Ve2"), ":");
  t2$Ve1 <- as.numeric(t2$Ve1);t2$Ve2 <- as.numeric(t2$Ve2);
  t2 <- t2 %>% separate(`KinshipSigma1:KinshipSigma2`, c("Vk1", "Vk2"), ":")
  t2$Vk1 <- as.numeric(t2$Vk1);t2$Vk2 <- as.numeric(t2$Vk2);
  t2 <- t2 %>% separate(`HouseholdSigma1:HouseholdSigma2`, c("Vh1", "Vh2"), ":")
  t2$Vh1 <- as.numeric(t2$Vh1);t2$Vh2 <- as.numeric(t2$Vh2);
  t2$pheno <- rownames(t2)
  t2 <- t2 %>% separate(pheno, c("pheno1", "pheno2"), "\\.")
  t2$i1 <- as.numeric(match(t2$pheno1,phenotypes))
  t2$i2 <- as.numeric(match(t2$pheno2,phenotypes))
  t2$cSP <- t2$Ve1 * t2$Ve2 * t2$E + t2$Vk1 * t2$Vk2 * t2$K
  t2$Sp <- as.numeric(as.character(t2$cor))
  t2$SpPval <- as.numeric(as.character(t2$`cor-p.val`));
  
  
  res$K[cbind(t2$i1,t2$i2)] <- t2$K; res$K[cbind(t2$i2,t2$i1)] <- t2$K; 
  res$HH[cbind(t2$i1,t2$i2)] <- t2$HH; res$HH[cbind(t2$i2,t2$i1)] <- t2$HH; 
  res$Fisher_K_pval[cbind(t2$i1,t2$i2)] <- t2$Fisher_k_pval_adj; res$Fisher_K_pval[cbind(t2$i2,t2$i1)] <- t2$Fisher_k_pval_adj; 
  res$Fisher_HH_pval[cbind(t2$i1,t2$i2)] <- t2$Fisher_HH_pval_adj; res$Fisher_HH_pval[cbind(t2$i2,t2$i1)] <- t2$Fisher_HH_pval_adj; 
  res$NormK[cbind(t2$i1,t2$i2)] <- t2$NormK; res$NormK[cbind(t2$i2,t2$i1)] <- t2$NormK; 
  res$NormHH[cbind(t2$i1,t2$i2)] <- t2$NormHH; res$NormHH[cbind(t2$i2,t2$i1)] <- t2$NormHH; 
  res$Bootstrap_K_pval[cbind(t2$i1,t2$i2)] <- t2$Bootstrap_k_pval_adj; res$Bootstrap_K_pval[cbind(t2$i2,t2$i1)] <- t2$Bootstrap_k_pval_adj; 
  res$Bootstrap_HH_pval[cbind(t2$i1,t2$i2)] <- t2$Bootstrap_HH_pval_adj; res$Bootstrap_HH_pval[cbind(t2$i2,t2$i1)] <- t2$Bootstrap_HH_pval_adj; 
  res$n[cbind(t2$i1,t2$i2)] <- t2$n; res$n[cbind(t2$i2,t2$i1)] <- t2$n; 
  res$BestpvalK[cbind(t2$i1,t2$i2)] <- t2$BestpvalK; res$BestpvalK[cbind(t2$i2,t2$i1)] <- t2$BestpvalK; 
  res$BestpvalHH[cbind(t2$i1,t2$i2)] <- t2$BestpvalHH; res$BestpvalHH[cbind(t2$i2,t2$i1)] <- t2$BestpvalHH; 
  res$E[cbind(t2$i1,t2$i2)] <- t2$E; res$E[cbind(t2$i2,t2$i1)] <- t2$E; 
  res$KSig1[cbind(t2$i1,t2$i2)] <- t2$Vk1; res$KSig1[cbind(t2$i2,t2$i1)] <- t2$Vk1; 
  res$KSig2[cbind(t2$i1,t2$i2)] <- t2$Vk2; res$KSig2[cbind(t2$i2,t2$i1)] <- t2$Vk1; 
  res$ESig1[cbind(t2$i1,t2$i2)] <- t2$Ve1; res$ESig1[cbind(t2$i2,t2$i1)] <- t2$Ve2; 
  res$ESig2[cbind(t2$i1,t2$i2)] <- t2$Ve2; res$ESig2[cbind(t2$i2,t2$i1)] <- t2$Ve1; 
  res$HHSig1[cbind(t2$i1,t2$i2)] <- t2$Vh1; res$HHSig1[cbind(t2$i2,t2$i1)] <- t2$Vh2; 
  res$HHSig2[cbind(t2$i1,t2$i2)] <- t2$Vh2; res$HHSig2[cbind(t2$i2,t2$i1)] <- t2$Vh1; 
  res$cSp[cbind(t2$i1,t2$i2)] <- t2$cSP; res$cSp[cbind(t2$i2,t2$i1)] <- t2$cSP; 
  res$Sp[cbind(t2$i1,t2$i2)] <- t2$Sp; res$Sp[cbind(t2$i2,t2$i1)] <- t2$Sp; 
  res$SpPval[cbind(t2$i1,t2$i2)] <- t2$SpPval; res$SpPval[cbind(t2$i2,t2$i1)] <- t2$SpPval; 
  res$Bootstrap_K_sd[cbind(t2$i1,t2$i2)] <- t2$Bootstrap_SdK; res$Bootstrap_K_sd[cbind(t2$i2,t2$i1)] <- t2$Bootstrap_SdK; 
  res$Bootstrap_HH_sd[cbind(t2$i1,t2$i2)] <- t2$Bootstrap_SdHH; res$Bootstrap_HH_sd[cbind(t2$i2,t2$i1)] <- t2$Bootstrap_SdHH; 
}

if (N<NN){
  print(paste0("Processed only ",N,"files out of ",NN,". Can't find the rest..."))
}else{
  print("Saving the resulting files...")
}  

tfname <- gsub(".txt",paste0("_K.txt"),outfilename)
res2<-data.frame(res$K)
rownames(res2)<-phenotypes
colnames(res2)<-phenotypes
write.table(res2, file = tfname, sep = "\t", col.names=NA)

tfname <- gsub(".txt",paste0("_HH.txt"),outfilename)
res2<-data.frame(res$HH)
rownames(res2)<-phenotypes
colnames(res2)<-phenotypes
write.table(res2, file = tfname, sep = "\t", col.names=NA)

tfname <- gsub(".txt",paste0("_NormK.txt"),outfilename)
res2<-data.frame(res$NormK)
rownames(res2)<-phenotypes
colnames(res2)<-phenotypes
write.table(res2, file = tfname, sep = "\t", col.names=NA)

tfname <- gsub(".txt",paste0("_NormHH.txt"),outfilename)
res2<-data.frame(res$NormHH)
rownames(res2)<-phenotypes
colnames(res2)<-phenotypes
write.table(res2, file = tfname, sep = "\t", col.names=NA)

tfname <- gsub(".txt",paste0("_Fisher_K_pval.txt"),outfilename)
res2<-data.frame(res$Fisher_K_pval)
rownames(res2)<-phenotypes
colnames(res2)<-phenotypes
write.table(res2, file = tfname, sep = "\t", col.names=NA)

tfname <- gsub(".txt",paste0("_Fisher_HH_pval.txt"),outfilename)
res2<-data.frame(res$Fisher_HH_pval)
rownames(res2)<-phenotypes
colnames(res2)<-phenotypes
write.table(res2, file = tfname, sep = "\t", col.names=NA)

tfname <- gsub(".txt",paste0("_Bootstrap_K_pval.txt"),outfilename)
res2<-data.frame(res$Bootstrap_K_pval)
rownames(res2)<-phenotypes
colnames(res2)<-phenotypes
write.table(res2, file = tfname, sep = "\t", col.names=NA)

tfname <- gsub(".txt",paste0("_Bootstrap_HH_pval.txt"),outfilename)
res2<-data.frame(res$Bootstrap_HH_pval)
rownames(res2)<-phenotypes
colnames(res2)<-phenotypes
write.table(res2, file = tfname, sep = "\t", col.names=NA)

tfname <- gsub(".txt",paste0("_BestpvalK.txt"),outfilename)
res2<-data.frame(res$BestpvalK)
rownames(res2)<-phenotypes
colnames(res2)<-phenotypes
write.table(res2, file = tfname, sep = "\t", col.names=NA)

tfname <- gsub(".txt",paste0("_BestpvalHH.txt"),outfilename)
res2<-data.frame(res$BestpvalHH)
rownames(res2)<-phenotypes
colnames(res2)<-phenotypes
write.table(res2, file = tfname, sep = "\t", col.names=NA)

tfname <- gsub(".txt",paste0("_KSig1.txt"),outfilename)
res2<-data.frame(res$KSig1)
rownames(res2)<-phenotypes
colnames(res2)<-phenotypes
write.table(res2, file = tfname, sep = "\t", col.names=NA)

tfname <- gsub(".txt",paste0("_KSig2.txt"),outfilename)
res2<-data.frame(res$KSig2)
rownames(res2)<-phenotypes
colnames(res2)<-phenotypes
write.table(res2, file = tfname, sep = "\t", col.names=NA)

tfname <- gsub(".txt",paste0("_ESig1.txt"),outfilename)
res2<-data.frame(res$ESig1)
rownames(res2)<-phenotypes
colnames(res2)<-phenotypes
write.table(res2, file = tfname, sep = "\t", col.names=NA)

tfname <- gsub(".txt",paste0("_ESig2.txt"),outfilename)
res2<-data.frame(res$ESig2)
rownames(res2)<-phenotypes
colnames(res2)<-phenotypes
write.table(res2, file = tfname, sep = "\t", col.names=NA)

tfname <- gsub(".txt",paste0("_HHSig1.txt"),outfilename)
res2<-data.frame(res$HHSig1)
rownames(res2)<-phenotypes
colnames(res2)<-phenotypes
write.table(res2, file = tfname, sep = "\t", col.names=NA)

tfname <- gsub(".txt",paste0("_HHSig2.txt"),outfilename)
res2<-data.frame(res$HHSig2)
rownames(res2)<-phenotypes
colnames(res2)<-phenotypes
write.table(res2, file = tfname, sep = "\t", col.names=NA)

tfname <- gsub(".txt",paste0("_E.txt"),outfilename)
res2<-data.frame(res$E)
rownames(res2)<-phenotypes
colnames(res2)<-phenotypes
write.table(res2, file = tfname, sep = "\t", col.names=NA)

tfname <- gsub(".txt",paste0("_n.txt"),outfilename)
res2<-data.frame(res$n)
rownames(res2)<-phenotypes
colnames(res2)<-phenotypes
write.table(res2, file = tfname, sep = "\t", col.names=NA)

tfname <- gsub(".txt",paste0("_Bootstrap_K_sd.txt"),outfilename)
res2<-data.frame(res$Bootstrap_K_sd)
rownames(res2)<-phenotypes
colnames(res2)<-phenotypes
write.table(res2, file = tfname, sep = "\t", col.names=NA)

tfname <- gsub(".txt",paste0("_Bootstrap_HH_sd.txt"),outfilename)
res2<-data.frame(res$Bootstrap_HH_sd)
rownames(res2)<-phenotypes
colnames(res2)<-phenotypes
write.table(res2, file = tfname, sep = "\t", col.names=NA)

tfname <- gsub(".txt",paste0("_calcSpearman.txt"),outfilename)
res2<-data.frame(res$cSp)
rownames(res2)<-phenotypes
colnames(res2)<-phenotypes
write.table(res2, file = tfname, sep = "\t", col.names=NA)

tfname <- gsub(".txt",paste0("_Spearman.txt"),outfilename)
res2<-data.frame(res$Sp)
rownames(res2)<-phenotypes
colnames(res2)<-phenotypes
write.table(res2, file = tfname, sep = "\t", col.names=NA)

tfname <- gsub(".txt",paste0("_SpearmanPvalue.txt"),outfilename)
res2<-data.frame(res$SpPval)
rownames(res2)<-phenotypes
colnames(res2)<-phenotypes
write.table(res2, file = tfname, sep = "\t", col.names=NA)

print("Done.")
