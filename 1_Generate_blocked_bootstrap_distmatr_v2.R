library(hash)
library(config)
args = commandArgs(trailingOnly=TRUE)
if(length(args)<1){
  stop("Config file not provided")
  # config <- config::get(file = "/data/linkage/HCHS_SOL/Projects/2021_gencor_cog_sleep/Code/runGenCor_ALL_Base_v1.config")
}else{
  config <- config::get(file = args[1])
}

phenotypes <- readRDS(config$phenotypes_filename)
infilenames <- c(config$kinship_matrix_filename,config$household_matrix_filename)
outfilenames <- c(config$kinship_matrix_filename_blocked_lookup,config$household_matrix_filename_blocked_lookup)

for(fid in 1:length(infilenames)){ # fid = 1
  print(paste0("Processing ",infilenames[fid]," with relatedness distance of ",config$relatedness_distance))
  mmat <- readRDS(infilenames[fid])
  iids <- paste0(phenotypes[[config$idfield_name]])
  iids <- iids[iids %in% colnames(mmat)]
  mmat <- mmat[iids,iids]
  
  dimensions = colnames(mmat)
  m1 <- round(mmat,3)
  thres <- config$relatedness_distance
  m1[m1<thres] <- 0
  length(which(m1>0))
  # test <- as.dist(m1)
  # plot(hclust(test))
  
  #We define a block if any of the participants there are correlated above some threshold
  #First we identify for each person who is he correlated with
  grps <- hash()
  lkup <-hash()
  cn <- c()
  for(i in 1:length(dimensions)){
    if(is.null(lkup[[paste0(i)]])){
      lkup[[paste0(i)]] <- i
    }
    tl <- paste0(lkup[[paste0(i)]])
    connected <- which(m1[i,]>0)
    grps[[tl]] <- unique(c(grps[[tl]],i,connected))
    cn <- append(cn,length(connected))
    for (j in connected) lkup[[paste0(j)]] <- tl
  }
  
  gn <- c()
  for(k in keys(grps)) gn <- append(gn,length(grps[[k]]))
  print(paste0("Average person is connected to ", mean(cn)," people (max-",
               max(cn),",min-",min(cn),")"))
  
  # first every one person group is OK
  lookup <- hash()
  for (k in keys(grps)){
    if(length(grps[[k]])==1){
      lookup[[dimensions[grps[[k]]]]] <- dimensions[grps[[k]]]; del(k,grps)
    }
  }
  print(paste0(length(keys(lookup))," people with no relatedness"))
  
  #now merge - if two groups share even one person, they are related and need to be merged
  mgrps <-hash()
  while(length(keys(grps))>0){
    k<-keys(grps)[1]
    fnd=0;fndgrps <- c()
    for (k2 in keys(grps)){ # k2=keys(grps)[2]
      if (k==k2) next
      if ( length(intersect(grps[[k]],grps[[k2]]))>1 ){
        fnd=1; fndgrps <- append(fndgrps,k2)
        # print(paste0("Merging ",k," and ",k2," intersect ",length(intersect(grps[[k]],grps[[k2]]))))
        # print(grps[[k]]);print(grps[[k2]])
        grps[[k]] <- unique(c(grps[[k]],grps[[k2]]))
      }
    }
    if(!fnd){  #unique
      # print(paste0(length(grps)))
      mgrps[[k]] <- grps[[k]]
      del(k,grps)
    }else{
      for(k2 in fndgrps) del(k2,grps)
    }
  } 
  
  
  gn2 <- c()
  for(k in keys(mgrps)) gn2 <- append(gn2,length(mgrps[[k]]))
  print(paste0(" After merge there are ",length(mgrps)," groups (mean-",mean(gn2),",max-",max(gn2),
               ",min-",min(gn2),")"))
  #now build a lookup table for each sample that would tell which samples must go in a block with it
  for(k in keys(mgrps)){
    for(v in mgrps[[k]]){
      # if(v %in% keys(lookup)) print(paste0("Problem with ",v))
      lookup[[dimensions[v]]]<- dimensions[mgrps[[k]]]
    }
  }
  
  save(lookup, file=outfilenames[fid])
}
