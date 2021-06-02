library(R.utils)
library(GENESIS)
library(config)
args = commandArgs(trailingOnly=TRUE)
if(length(args)<1){
  stop("Config file not provided")
  # config <- config::get(file = "/data/linkage/HCHS_SOL/Projects/2021_gencor_cog_sleep/Code/runGenCor_ALL_Base_v1.config")
}else{
  config <- config::get(file = args[1])
}

infilenames <- c(config$kinship_matrix_filename,config$household_matrix_filename)
outfilenames <- c(config$kinship_matrix_filename_GCTA,config$household_matrix_filename_GCTA)

for(fid in 1:length(infilenames)){ # fid = 1

  # non-related individuals : if don't exist will be calculated
  Unrelated3rdDeg_fname <- config$Unrelated_individuals_filename
  #################################### Convert the GRM to GCTA format ############################
  if(config$CONVERT_MATRICES_toGCTA){
    infilename <- infilenames[fid]
    print(paste0("Converting the matrix ", infilename, " to GCTA format"))

    mat_grm <- readRDS(infilename)
    cn <- colnames(mat_grm);
    cn <- paste0(cn[!is.na(cn)])
    mat_grm <- mat_grm[cn,cn]
    if(!file.exists(config$Unrelated_individuals_filename)){
      print(paste0("Generating a list of unrelated individuals at the ", config$relatedness_distance, " kinship distance"))
      partres <- pcairPartition(mat_grm,kin.thresh = config$relatedness_distance)
      Unrelated3rdDeg <- partres$unrels
      save(Unrelated3rdDeg,file=config$Unrelated_individuals_filename)
    }else{
      print(paste0("Loading an existing list of unrelated individuals", config$Unrelated_individuals_filename))
      load(config$Unrelated_individuals_filename)
    }
    mat_grm <- mat_grm[Unrelated3rdDeg,Unrelated3rdDeg]
    
    print("Done. Converting to GCTA format ...")
  
    SNPnum <- 638486    # TOPMed freeze 8
    n <- dim(mat_grm)[1]
    ns <- (n+1)*n/2
    res <- matrix(0, nrow = ns, ncol = 4)
    rw <- 1
    for(i in 1:n){
      for(j in 1:i){
        res[rw,] <- c(i,j,SNPnum,mat_grm[i,j])
        rw <- rw+1
      }
    }
    
    print("Done. Saving ...")
  
    write.table(res,file = gzfile(paste0(outfilenames[fid],".grm.gz")),quote = FALSE, sep = "\t",
                row.names = FALSE,col.names = FALSE)
    write.table(cbind(colnames(mat_grm),colnames(mat_grm)),
                file = paste0(outfilenames[fid],".grm.id"),quote = FALSE, sep = "\t",
                row.names = FALSE,col.names = FALSE)
  }
}

if(config$CONVERT_MATRICES_toGCTA){
  #write two-matrix file
  fileConn<-file(config$two_matrix_GCTA_filename)
  writeLines(outfilenames, fileConn)
  close(fileConn)
}

#################################### Convert the phenotype data to GCTA format ############################
phenotypes <- readRDS(config$phenotypes_filename)

covariates0 = config$factor_covariates
covariates1 <- config$numeric_covariates
covariates = c(covariates0,covariates1)

IBD_order_TOPMed_NormK <- config$phenotypes_fieldnames

rownames(phenotypes) <- phenotypes[[config$idfield_name]]
phenotypes <- phenotypes[,c(config$idfield_name,covariates,config$phenotypes_fieldnames)]
fstart <- length(covariates)+2

print("Converting phenotypes ...")

#phenotypes
cn <- IBD_order_TOPMed_NormK
flds <- c(config$idfield_name,config$idfield_name,cn)

pdata <- phenotypes[,flds]
write.table(pdata,
            file = paste0(config$phenotypes_GCTA_prefix,".pheno"),quote = FALSE, sep = "\t",
            row.names = FALSE,col.names = FALSE)
write.table(IBD_order_TOPMed_NormK,file = paste0(config$phenotypes_GCTA_prefix,".pheno.names"),quote = FALSE, sep = "\t",row.names = TRUE,col.names = FALSE)


#factor covariates

#make sure its not single value
gc <- c()
for(c in covariates0){ # c = covariates0[3]
  if(length(unique(phenotypes[[c]]))>1){
    gc <- c(gc,c)
  }
}

flds <- c(config$idfield_name,config$idfield_name,gc)
pdata <- phenotypes[,flds]
iids <- complete.cases(pdata)
pdata <- pdata[iids,]
write.table(pdata,
            file = paste0(config$phenotypes_GCTA_prefix,".covar"),quote = FALSE, sep = "\t",
            row.names = FALSE,col.names = FALSE)

#quantitative covariates
flds <- c(config$idfield_name,config$idfield_name,covariates1)
pdata <- phenotypes[,flds]
write.table(pdata,
            file = paste0(config$phenotypes_GCTA_prefix,".qcovar"),quote = FALSE, sep = "\t",
            row.names = FALSE,col.names = FALSE)


#################################### ranknorm the phenotype  and readjust covariates + convert to GCTA format ############################
print("Done. Adjust + Ranknorm + convert phenotypes ...")

#adjust for covariates
print("Adjusting raw data for Covariates...")
pdata <- phenotypes
phenotypesN <- colnames(phenotypes)
for (i in fstart:length(phenotypesN)){ # i <- fstart
  pheno.1 = phenotypesN[i]
  print(pheno.1)
  model.formula <- as.formula(paste(paste(pheno.1,"~"), paste0(covariates,collapse="+")))
  if(length(which(!is.na(pdata[pheno.1])))<50){
    print("All NA, Skippig")
    pdata[,i]<-NA
    next
  }
  r2 <- pdata[,i]
  result <- tryCatch({
    mod <- lm(model.formula, data = pdata[,c(pheno.1,covariates)],na.action="na.exclude")
    r2 <- resid(mod,na.action="na.exclude")
  }, warning = function(war) {
    print("Adjusting failed")
  }, error = function(err) {
    print("Adjusting failed")
  }) # END tryCatch
  
  pdata[,i]<-r2
}

#ranknorm
print("Rank-normalizing...")
for (i in fstart:length(pdata)) {
  x <-pdata[,i]
  pdata[,i] <- scale(qnorm(rank(x, ties.method = "random",na.last="keep")/(length(x)+1)))
}

#phenotypes
print("Saving ranknorm...")
cn <- IBD_order_TOPMed_NormK
flds <- c(config$idfield_name,config$idfield_name,cn)
pdata <- pdata[,flds]
write.table(pdata,
            file = paste0(config$phenotypes_GCTA_prefix,".ranknorm.pheno"),quote = FALSE, sep = "\t",
            row.names = FALSE,col.names = FALSE)

#final run script
fileConn<-file(config$run_GCTA_script_filename)
iids <- 1:length(IBD_order_TOPMed_NormK)
rrows <- paste0("bsub -q ", config$LSF_QUEUE_NAME ," -J ",iids, ' -n 2 -R "rusage[mem=',config$LSF_MEMORY,']" -e ',config$GCTA_cluster_filename_prefix, "_",iids,".e",
                " -o ",config$GCTA_cluster_filename_prefix, "_",iids,".o",' "/PHShome/me704/gcta_1.93.2beta/gcta64 --reml --pheno ',
                paste0(config$phenotypes_GCTA_prefix,".ranknorm.pheno"),' --mpheno ',iids," --covar ",
                paste0(config$phenotypes_GCTA_prefix,".covar")," --qcovar ",paste0(config$phenotypes_GCTA_prefix,".qcovar"),
                " --mgrm-gz ",config$two_matrix_GCTA_filename," --thread-num 8 --out ",paste0(config$phenotypes_GCTA_prefix,"_",iids,'"'))
writeLines(rrows, fileConn)
close(fileConn)
Sys.chmod(config$run_GCTA_script_filename,"755")

print("Done. ")
