require(RColorBrewer)
library(corrplot)
library(plotrix)
library(ggchicklet)
library(qgraph)
library(corpcor)
library(gdata)
library(reshape2)
library(ggpattern)
library(igraph)
library(ggplot2)
library(config)
args = commandArgs(trailingOnly=TRUE)
if(length(args)<1){
  stop("Config file not provided")
  # config <- config::get(file = "/data/linkage/HCHS_SOL/Projects/2021_gencor_cog_sleep/Code/runGenCor_ALL_Base_v1.config")
}else{
  config <- config::get(file = args[1])
}

CALC_ORDER = config$CALC_ORDER
CALC_ORDER_QGRAPH = config$CALC_ORDER_QGRAPH

tthres = config$pvalueThreshold
normthres = config$correlationThreshold
peoplethres = config$peopleNumThreshold

SOL_Order1 <- config$phenotypes_fieldnames

groupsNames <- list()
groupsNames[[config$phenotypes_group1_Name]] <- config$phenotypes_group1
groupsNames[[config$phenotypes_group2_Name]] <- config$phenotypes_group2

groups <- list()
groups[[config$phenotypes_group1_Name]] <- groupsNames[[config$phenotypes_group1_Name]]
groups[[config$phenotypes_group2_Name]] <- groupsNames[[config$phenotypes_group2_Name]]

for (j in names(groups)){ #j=names(groups)[1]
  groups[[j]] <- match(groups[[j]],SOL_Order1)
}

groups <- groups[c(config$phenotypes_group1_Name,config$phenotypes_group2_Name)]

ttitle <- c("K.txt","HH.txt")

workdirname <- config$workdirname
outfileprefix <- config$outfileprefix
ofilename = paste0(outfileprefix,"_combine.txt")
outfilename = paste0(workdirname,ofilename)

fls <- c(gsub(".txt",paste0("_GCTA_VCs_NormK.txt"),outfilename),
         gsub(".txt",paste0("_K.txt"),outfilename),
         gsub(".txt",paste0("_BestpvalK.txt"),outfilename),
         gsub(".txt",paste0("_GCTA_VCs_calcSpearman.txt"),outfilename),
         gsub(".txt",paste0("_Spearman.txt"),outfilename),
         gsub(".txt",paste0("_GCTA_VCs_calcSpearmanPvalue.txt"),outfilename),
         gsub(".txt",paste0("_n.txt"),outfilename)
        )

pheno_names <- paste0(config$phenotypes_GCTA_prefix,".pheno.names")
VCprefix = config$phenotypes_GCTA_prefix

layoutfile = config$QGRAPH_LAYOUT_FILE
ttls <- c("Phenotypic Correlation","Genetic Correlation","Normalized Genetic Correlation","Genetic Heritability")
tflds <- SOL_Order1

gri <- c()
gr <- ""

toutfname <- gsub(".txt",".pdf",outfilename)
if(length(gri)>0){
  toutfname <- gsub("_SOLOrder_",paste0("_SOLOrder_",gr,"_"),toutfname)
}
toutfname <- gsub(".pdf",paste0("_filtPval",tthres,"_filtThres",normthres,"_pplThres",peoplethres,".pdf"),toutfname)

print(toutfname)

pdf(toutfname, height=20, width=20)

###########################################################################################################
for (i in 1:length(ttitle)){ # i=1
  print(ttitle[i])

  flsI <- c()
  for (j in 1:length(fls)) {
    flsI <- c(flsI, gsub("K.txt",ttitle[i],fls[j]))
  }
  if(ttitle[i]=="HH.txt"){
    ttls <- gsub("GenCorr","HouseholdCorr",ttls)
  }
  
  fNormK=flsI[1]
  fK=flsI[2]
  fpval = flsI[3]
  fcSp = flsI[4]
  fSp = flsI[5]
  fSpPval = flsI[6]
  fN = flsI[7]

  if(ttitle[i]=="HH.txt"){
    ttls <- gsub("Genetic","Household",ttls)
  }
  
  ##### Read calculated Spearman ######
  datacSP<-read.table(fcSp,header = TRUE,sep = '\t')
  cn<-colnames(datacSP)[2:length(datacSP)];rn <- as.character(datacSP[,1]);datacSP<-datacSP[,2:length(datacSP)];
  colnames(datacSP)<-cn;rownames(datacSP)<-rn
  if(!CALC_ORDER){
    flds<-tflds
    datacSP <- datacSP[flds,flds]
  }
  
  mcSP<-as.matrix(datacSP)
  mcSP[mcSP < -10]<- NA
  mcSP[mcSP < -1]<- -1
  mcSP[mcSP > 1]<- 1
  
  #read calcSpearman pval
  datacSP_Pval<-read.table(fSpPval,header = TRUE,sep = '\t')
  cn<-colnames(datacSP_Pval)[2:length(datacSP_Pval)];rn <- as.character(datacSP_Pval[,1]);datacSP_Pval<-datacSP_Pval[,2:length(datacSP_Pval)];
  colnames(datacSP_Pval)<-cn;rownames(datacSP_Pval)<-rn
  if(!CALC_ORDER){
    flds<-tflds
    datacSP_Pval <- datacSP_Pval[flds,flds]
  }
  mcSP_Pval <- as.matrix(datacSP_Pval);
  
  ##### Read Spearman ######
  dataSP<-read.table(fSp,header = TRUE,sep = '\t')
  cn<-colnames(dataSP)[2:length(dataSP)];rn <- as.character(dataSP[,1]);dataSP<-dataSP[,2:length(dataSP)];
  colnames(dataSP)<-cn;rownames(dataSP)<-rn
  if(!CALC_ORDER){
    flds<-tflds
    dataSP <- dataSP[flds,flds]
  }
  
  mSP<-as.matrix(dataSP)
  mSP[mSP < -10]<- NA
  mSP[mSP < -1]<- -1
  mSP[mSP > 1]<- 1
  mSP[is.na(mSP)]<-0
  
  ##### Complete calculated Spearman with Spearman where NA ######
  mcSP <- ifelse(is.na(mcSP), mSP, mcSP)
  
  ##### Read number of people ######
  dataN<-read.table(fN,header = TRUE,sep = '\t')
  cn<-colnames(dataN)[2:length(dataN)];rn <- as.character(dataN[,1]);dataN<-dataN[,2:length(dataN)];
  colnames(dataN)<-cn;rownames(dataN)<-rn
  if(!CALC_ORDER){
    flds<-tflds
    dataN <- dataN[flds,flds]
  }
  
  mN<-as.matrix(dataN)
  mN[mN < 0]<- NA
  mN[is.na(mN)]<-0
  
  ##### Read NormK file ######
  dataNormK<-read.table(fNormK,header = TRUE,sep = '\t')
  cn<-colnames(dataNormK)[2:length(dataNormK)];rn <- as.character(dataNormK[,1]);dataNormK<-dataNormK[,2:length(dataNormK)];
  colnames(dataNormK)<-cn;rownames(dataNormK)<-rn
  if(!CALC_ORDER){
    flds<-tflds
    dataNormK <- dataNormK[flds,flds]
  }
  
  mNormK<-as.matrix(dataNormK)
  mNormK[mNormK < -10]<- NA
  mNormK[mNormK < -1]<- -1
  mNormK[mNormK > 1]<- 1
  mNormK[is.na(mNormK)]<-0
  
  n <- nrow(mNormK)
  
  dataK_Pval<-read.table(fpval,header = TRUE,sep = '\t')
  cn<-colnames(dataK_Pval)[2:length(dataK_Pval)];rn <- as.character(dataK_Pval[,1]);dataK_Pval<-dataK_Pval[,2:length(dataK_Pval)];
  colnames(dataK_Pval)<-cn;rownames(dataK_Pval)<-rn
  if(!CALC_ORDER){
    flds<-tflds
    dataK_Pval <- dataK_Pval[flds,flds]
  }
  mK_Pval <- as.matrix(dataK_Pval);
  
  ##### Read K file ######
  dataK<-read.table(fK,header = TRUE,sep = '\t')
  cn<-colnames(dataK)[2:length(dataK)];rn <- as.character(dataK[,1]);dataK<-dataK[,2:length(dataK)];
  colnames(dataK)<-cn;rownames(dataK)<-rn
  if(!CALC_ORDER){
    flds<-tflds
    dataK <- dataK[flds,flds]
  }
  
  mK<-as.matrix(dataK)
  mK[mK < -10]<- NA
  mK[mK < -1]<- -1
  mK[mK > 1]<- 1
  mK[is.na(mK)]<-0
  # diag(mK) <- 0
  
  #read GCTA heritabilities
  workdirname <- strsplit(config$phenotypes_GCTA_prefix,"/",fixed = T)[[1]]
  workdirname <- paste0(paste0(c(workdirname[1:(length(workdirname)-1)]),collapse = "/"),"/")
  ddata <- data.frame(
    name=factor(SOL_Order1,levels = rev(SOL_Order1)),
    value=rep(NA,length(SOL_Order1)),
    sd=rep(NA,length(SOL_Order1))
  )
  rownames(ddata) <- ddata$name
  phenotypes <- read.table(pheno_names,sep="\t",header = F, row.names = 1)
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
  
  rownames(ddata) <- ddata$name;ddata <- ddata[tflds,]
  ddata <- data.frame(
    name=factor(tflds,levels = rev(tflds)),
    value=ddata$value,
    sd=ddata$sd
  )
  
  # # print(gsub("NormK.txt","Heritabilities.txt",fNormK))
  # if(ttitle[i]=="K.txt"){
  #   write.table(ddata,gsub("NormK.txt","K_Heritabilities.txt",fNormK),row.names = F, sep = "\t")
  # }else{
  #   write.table(ddata,gsub("NormHH.txt","HH_Heritabilities.txt",fNormK),row.names = F, sep = "\t")
  # }
  ################################################## make filters 
  
  mNf <- abs(mN)
  mNf[mNf < peoplethres] <- 0
  mNf[mNf >= peoplethres] <- 1
  
  # make a filter out of pvalue
  mPf <- mK_Pval
  mPf[is.na(mPf)]<- 1
  mPf[mPf== -100]<- 1
  mPf[mPf >= tthres] <- 1
  mPf[mPf < tthres] <- 0
  mPf <- !mPf
  
  # make a filter out of NormK
  mNKf <- abs(mNormK)
  mNKf[mNKf >= normthres] <- 1
  mNKf[mNKf < normthres] <- 0
  
  # make a filter out of K
  mKf <- abs(mK)
  # mKf[mKf > 0] <- 1
  mKf[mKf >= normthres] <- 1
  mKf[mKf < normthres] <- 0
  
  #apply the filters
  final_NK <- mNormK * mPf * mNKf * mNf
  final_NK <- round(final_NK,2)
  
  final_K <- mK * mPf * mKf * mNf
  final_K <- round(final_K,2)
  
  
  ############################stats for groups
  if(length(names(groupsNames))>2){
    sdresPval <- matrix(1,nrow = length(groupsNames),ncol = length(groupsNames))
    sdres <- matrix(0,nrow = length(groupsNames),ncol = length(groupsNames))
    colnames(sdres) <- names(groupsNames); rownames(sdres) <- names(groupsNames);
    selMat <- final_K
    for(d in 1:length(groupsNames)){ # d=1
      for(s in 1:length(groupsNames)){ # s=3
        # if(s==d) next
        tsdr <- selMat[groupsNames[[d]],groupsNames[[s]]]
        sdres[d,s] <- length(which(tsdr!=0))
      }
    }
    ffnk <- selMat; ffnk[ffnk!=0] <- 1
    TraitProbs <- rowSums(ffnk)/sum(rowSums(ffnk))
    # #permutation test - sample pairs of groups (according to probability?)
    # #do that 1000 times; for each pair calculate probability
    #
    if(length(which(is.na(TraitProbs)))<1){
      NREPSSIM <- 1000
      sdresperm <- list()
      for(d in 1:NREPSSIM){ #d=1
        # no self allowed
        # g1 <- sample(x=1:dim(selMat)[1],size = sum(rowSums(ffnk))*2,replace = T,prob = TraitProbs)
        # g2 <- sample(x=1:dim(selMat)[1],size = sum(rowSums(ffnk))*2,replace = T,prob = TraitProbs)
        # gidx <- which(g1!=g2);g1 <- g1[gidx[1:sum(rowSums(ffnk))]];g2 <- g2[gidx[1:sum(rowSums(ffnk))]]
    
        # #self allowed
        g1 <- sample(x=1:dim(selMat)[1],size = sum(rowSums(ffnk)),replace = T,prob = TraitProbs)
        g2 <- sample(x=1:dim(selMat)[1],size = sum(rowSums(ffnk)),replace = T,prob = TraitProbs)
    
        tsdres <- matrix(0,nrow = dim(selMat)[1],ncol = dim(selMat)[1])
        colnames(tsdres) <- colnames(selMat); rownames(tsdres) <- rownames(selMat);
        for(s in 1:length(g1)){
          tsdres[g1[s],g2[s]] <- tsdres[g1[s],g2[s]]+1
        }
        tsdres2 <- matrix(0,nrow = length(groupsNames),ncol = length(groupsNames))
        colnames(tsdres2) <- names(groupsNames); rownames(tsdres2) <- names(groupsNames);
        for(w in 1:length(groupsNames)){ # w=1
          for(s in 1:length(groupsNames)){ # s=3
            # if(s==d) next
            tsdr <- tsdres[groupsNames[[w]],groupsNames[[s]]]
            tsdres2[w,s] <- length(which(tsdr!=0))
          }
        }
        sdresperm[[d]] <- tsdres2 #sum(rowSums(tsdres))
      }
      #
      colnames(sdresPval) <- names(groupsNames); rownames(sdresPval) <- names(groupsNames);
      for(d in 1:length(groupsNames)){ # d=1
        for(s in 1:length(groupsNames)){ # s=3
          dsi <- rep(0,NREPSSIM)
          for(w in 1:NREPSSIM){ #w=1
            dsi[w] <- sdresperm[[w]][d,s]
          }
          sdresPval[d,s] <- length(which(dsi>=sdres[d,s]))/NREPSSIM
        }
      }
    
      sdresres <- sdresPval; sdresres[lower.tri(sdresres)] <- 1;
      sdresres2 <- sdresres
      sdresres <- as.data.frame(as.table(sdresres))
      sdresres <- sdresres[sdresres$Freq <0.05,]
    }else{
      sdresPval <- matrix(1,nrow = length(groupsNames),ncol = length(groupsNames))
      colnames(sdresPval) <- names(groupsNames); rownames(sdresPval) <- names(groupsNames);
    }
  }else{
    sdresPval <- matrix(0,nrow = length(groupsNames),ncol = length(groupsNames))
  }
  # print(sdresres)
  

  ##########################################################################################  
  mcSPf <- mcSP_Pval
  mcSPf[is.na(mcSPf)]<- 1
  mcSPf[mcSPf== -100]<- 1
  mcSPf[mcSPf >= tthres] <- 1
  mcSPf[mcSPf < tthres] <- 0
  mcSPf <- !mcSPf
  final_mcSP <- mcSP * mcSPf * mNf
  
  
  # plot K 
  if(!CALC_ORDER){
    corrplot(final_K, method = "circle", p.mat = mK_Pval, type = "lower",
             tl.col = "black", tl.cex = 2, tl.srt = 45,  
             insig = "label_sig", sig.level = .05, pch.cex = 1, pch.col = "darkgrey",
             col = colorRampPalette(rev(brewer.pal(n = 10, name = "PuOr")))(200),bg = "azure1")
  }else{
    corrplot(final_K, method = "circle", p.mat = mK_Pval, type = "lower",
             tl.col = "black", tl.cex = 2, tl.srt = 45,  order = "hclust",
             insig = "label_sig", sig.level = .05, pch.cex = 1, pch.col = "darkgrey",
             col = colorRampPalette(rev(brewer.pal(n = 10, name = "PuOr")))(200),bg = "azure1")
  }
  mtext(text = ttls[2], side = 3, line = 1,cex = 3)
  
  
  # plot NormK
  if(!CALC_ORDER){
    corrplot(final_NK, method = "circle", p.mat = mK_Pval, type = "lower",
             tl.col = "black", tl.cex = 2, tl.srt = 45, 
             insig = "label_sig", sig.level = .05, pch.cex = 1, pch.col = "darkgrey",
             col = colorRampPalette(rev(brewer.pal(n = 10, name = "PuOr")))(200),bg = "azure1")
  }else{
    corrplot(final_NK, method = "circle", p.mat = mK_Pval, type = "lower",
             tl.col = "black", tl.cex = 2, tl.srt = 45, order = "hclust",
             insig = "label_sig", sig.level = .05, pch.cex = 1, pch.col = "darkgrey",
             col = colorRampPalette(rev(brewer.pal(n = 10, name = "PuOr")))(200),bg = "azure1")
  }
  mtext(text = ttls[3], side = 3, line = 1,cex = 3)
  
  # plot calc Spearman 
  if(!CALC_ORDER){
    corrplot(final_mcSP, method = "circle", type = "lower",
             tl.col = "black", tl.cex = 2, tl.srt = 45,  
             insig = "label_sig", sig.level = .05, pch.cex = 1, pch.col = "darkgrey",
             col = colorRampPalette(rev(brewer.pal(n = 10, name = "PuOr")))(200),bg = "azure1")
  }else{
    corrplot(final_mcSP, method = "circle", type = "lower",
             tl.col = "black", tl.cex = 2, tl.srt = 45,  order = "hclust",
             insig = "label_sig", sig.level = .05, pch.cex = 1, pch.col = "darkgrey",
             col = colorRampPalette(rev(brewer.pal(n = 10, name = "PuOr")))(200),bg = "azure1")
  }
  mtext(text = ttls[1], side = 3, line = 1,cex = 3)
  
  # plot the VCs
  dcols <- ddata$value; 
  oC <- brewer.pal(n = 9, name = "Oranges")
  oC <- colorRampPalette(oC)(20)
  obr <- 1/length(oC)
  for(ci in seq(length(oC),1,-1)){ #ci=8
    # print(paste0(ci*obr,"-",(ci-1)*obr))
    dcols[(ddata$value <= ci*obr) &(ddata$value > (ci-1)*obr)] <- oC[ci]
  }
  # dcols[ddata$value<=0.2] <- "darkorange1"; 
  # dcols[(ddata$value>0.2)&(ddata$value<=0.4)] <- "darkorange2";dcols[ddata$value>0.4] <- "darkorange3"; 
  mplt1 <- ggplot(ddata) +
    geom_chicklet(aes(x=name, y=value), radius = grid::unit(10, "pt"), fill=dcols,alpha=0.8) + 
    # geom_bar( aes(x=name, y=value), stat="identity",fill=dcols,alpha=0.8) +
    geom_errorbar( aes(x=name, ymin=value-sd, ymax=value+sd), width=0.4, colour="black", alpha=0.9, size=1.2) +
    scale_y_continuous(limits=c(-0.3,1.2))  + 
    labs(title=ttls[4])+
    coord_flip()
  mplt1 <- mplt1 + theme(text = element_text(size = 30))
  print(mplt1)
  
  
  ################################################## plot graphs
  #filter for one group of interest
  if(length(gri)>0){
    mxNK <- max(final_NK)
    fngri <- final_NK; fngri[fngri!=0] <- 0
    griKeep1 <- final_NK[gri,];griKeep2 <- final_NK[,gri]
    fngri[gri,] <- griKeep1;fngri[,gri] <- griKeep2
    final_NK <- fngri
    
    mxK <- max(final_K)
    fngri <- final_K; fngri[fngri!=0] <- 0
    griKeep1 <- final_K[gri,];griKeep2 <- final_K[,gri]
    fngri[gri,] <- griKeep1;fngri[,gri] <- griKeep2
    final_K <- fngri
    
    mxcSP <- max(final_mcSP)
    fngri <- final_mcSP; fngri[fngri!=0] <- 0
    griKeep1 <- final_mcSP[gri,];griKeep2 <- final_mcSP[,gri]
    fngri[gri,] <- griKeep1;fngri[,gri] <- griKeep2
    final_mcSP <- fngri
  }
  
  if(gr==""){
    gr="ALL"
  }
  #save the data
  write.table(final_mcSP,gsub("_calcSpearman.txt",
                              paste0("_calcSpearman_filtPval",tthres,"_filtThres",normthres,"_pplThres",peoplethres,"_",gr,".txt"),fcSp), sep = "\t")
  if(ttitle[i]=="K.txt"){
    write.table(final_NK,gsub("NormK.txt",paste0("NormK_filtPval",tthres,"_filtThres",normthres,"_pplThres",peoplethres,"_",gr,".txt"),fNormK), sep = "\t")
    write.table(final_NK,gsub("K.txt",paste0("K_filtPval",tthres,"_filtThres",normthres,"_pplThres",peoplethres,"_",gr,".txt"),fK), sep = "\t")
    write.table(sdresPval,gsub("K.txt",paste0("K_InteractionPvalueByGroup_filtPval",tthres,"_filtThres",normthres,"_pplThres",peoplethres,"_",gr,".txt"),fK), sep = "\t")
  }else{
    write.table(final_NK,gsub("NormHH.txt",paste0("NormHH_filtPval",tthres,"_filtThres",normthres,"_pplThres",peoplethres,"_",gr,".txt"),fNormK), sep = "\t")
    write.table(final_NK,gsub("HH.txt",paste0("HH_filtPval",tthres,"_filtThres",normthres,"_pplThres",peoplethres,"_",gr,".txt"),fK), sep = "\t")
    write.table(sdresPval,gsub("HH.txt",paste0("HH_InteractionPvalueByGroup_filtPval",tthres,"_filtThres",normthres,"_pplThres",peoplethres,"_",gr,".txt"),fK), sep = "\t")
  }
  
  diag(final_NK) <- 1
  diag(final_K) <- 1
  diag(final_mcSP) <- 1
  # 
  oC <- brewer.pal(n = 11, name = "Spectral")

  #plot network graphs
  if(ttitle[i]=="HH.txt"){
    if(CALC_ORDER_QGRAPH){
      CALC_ORDER_QGRAPH = F
      layoutfile = gsub("_HH.txt","_SOL_qgraph_res_v3.RData",fK)
    }
  }
  
  if(CALC_ORDER_QGRAPH){
    print("Generating new layout...")
    TOPMed_qgraph_res <- qgraph(final_K, label.cex = 1, label.scale = F,legend=F, vsize=2, color=oC, shape="ellipse", vsize2=2,
                                labels = colnames(final_NK), title = "Layout Generation, Disregard",DoNotPlot = TRUE,
                                posCol="darkorange3",negCol="purple",node.label.offset = c(0.5,-1.5),esize = 4,
                                groups=groups, curve = 0.5, curveAll = TRUE,edge.width = 3,label.font =2,
                                border.width = 2, border.color = "blue")
    if(ttitle[i]!="HH.txt"){
      save(TOPMed_qgraph_res, file = gsub("_K.txt","_SOL_qgraph_res_v3.RData",fK))
    }
  }
  if(!exists("TOPMed_qgraph_res")){
    if(file.exists(layoutfile)){
      print(paste0("Found previous layout, loading",layoutfile,"..."))
      load(layoutfile)
      if(!exists("TOPMed_qgraph_res")){
        TOPMed_qgraph_res <- TOPMed_qgraph_res_mutual
      }
    }
  }

  if(length(gri)>0){
    rratioNK <- mxNK / max(final_NK[final_NK<1]);rratioK <- mxK / max(final_K[final_K<1]);rratiocSP <- mxcSP / max(final_mcSP[final_mcSP<1]);
  }else{
    final_NKa <- abs(final_NK);final_Ka <- abs(final_K);final_mcSPa <- abs(final_mcSP);
    rratioNK <- 1 / max(final_NKa[final_NKa<1]);rratioK <- 1 / max(final_Ka[final_Ka<1]);rratiocSP <- 1 / max(final_mcSPa[final_mcSPa<1]);
  }
  
  scalefact <- ddata$value*8+1
  
  llab <- colnames(final_NK)
  qgK <- qgraph(final_K, layout = TOPMed_qgraph_res$layout, title = ttls[2],
         label.cex = 3, label.scale = F,legend=F, vsize=scalefact, color=oC, 
         labels = llab, shape="circle", 
         posCol="darkorange3",negCol="purple",node.label.offset = c(0.5,-1.5),esize = 1/rratioK,
         groups=groups, curve = 0.5, curveAll = TRUE,edge.width = 3,label.font =2,
         border.width = 1, border.color = "blue")

  qgNK <- qgraph(final_NK, layout = TOPMed_qgraph_res$layout, title = ttls[3],
         label.cex = 3, label.scale = F,legend=F, vsize=scalefact, color=oC, 
         labels = llab, shape="circle", 
         posCol="darkorange3",negCol="purple",node.label.offset = c(0.5,-1.5),esize = 1/rratioNK,
         groups=groups, curve = 0.5, curveAll = TRUE,edge.width = 3,label.font =2,
         border.width = 1, border.color = "blue")
  
  qgcSP <- qgraph(final_mcSP, layout = TOPMed_qgraph_res$layout, title = ttls[1],
         label.cex = 3, label.scale = F,legend=F, vsize=scalefact, color=oC, 
         labels = llab, shape="circle", 
         posCol="darkorange3",negCol="purple",node.label.offset = c(0.5,-1.5),esize = 1/rratiocSP,
         groups=groups, curve = 0.5, curveAll = TRUE,edge.width = 3,label.font =2,
         border.width = 1, border.color = "blue")
  
  final_GK <- sdresPval; final_GK[final_GK==0] <- 0.001;final_GK[final_GK > 0.05] <- 0; final_GK[final_GK!=0] <- 1 
  diag(final_GK) <- 0
  qgraph(final_GK, title = "Group Interaction",layout="circle", title="Group Connections",
        label.cex = 3, label.scale = F,legend=F, vsize=15, color=oC, 
        labels = colnames(final_GK), shape="circle", vTrans = 50,directed=F,
        posCol="black",negCol="purple",node.label.offset = c(0.5,-1.5),esize = 8,
        curve = 0.5, curveAll = TRUE,edge.width = 3,label.font =2,
              border.width = 1, border.color = "blue")

  # igK <- as.igraph(qgK)
  final_K0 <- abs(final_K); diag(final_K0)<-0;final_K0[final_K0>0]<-1
  igK <- graph_from_adjacency_matrix(final_K0,mode = "undirected")
  deg <- igraph::degree(igK, mode="all")
  final_K0 <- final_K; diag(final_K0)<-0
  qgK2 <- qgraph(final_K0, layout = TOPMed_qgraph_res$layout, title = "Node degree (number of outgoing connections)",
                 label.cex = 3, label.scale = F,legend=F, vsize=deg/2, color=oC, 
                 labels = llab, shape="circle",
                 posCol="grey",negCol="grey",node.label.offset = c(0.5,-1.5),esize = 1/rratioK,
                 groups=groups, curve = 0.5, curveAll = TRUE,edge.width = 3,label.font =2,
                 border.width = 1, border.color = "blue")
}

dev.off()

