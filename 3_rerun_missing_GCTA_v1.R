library(config)
args = commandArgs(trailingOnly=TRUE)
if(length(args)<1){
  stop("Config file not provided")
  # config <- config::get(file = "/data/linkage/HCHS_SOL/Projects/2021_gencor_cog_sleep/Code/runGenCor_ALL_Base_v1.config")
}else{
  config <- config::get(file = args[1])
}

print("Looking for missing GCTA heritability files...")
IBD_order_TOPMed_NormK <- config$phenotypes_fieldnames
iids <- c()
for (i in 1:length(IBD_order_TOPMed_NormK)){
  if(!file.exists(paste0(config$phenotypes_GCTA_prefix,"_",i,".hsq"))){
    iids <- c(iids,i)
    print(paste0("Missing ", paste0(config$phenotypes_GCTA_prefix,"_",i,".hsq")))
  }
}

if(length(iids)>0){
  fileConn<-file(config$rerun_GCTA_script_filename)
  rrows <- paste0("bsub -q ", config$LSF_QUEUE_NAME ," -J ",iids, ' -n 2 -R "rusage[mem=',config$LSF_MEMORY,']" -e ',config$GCTA_cluster_filename_prefix, "_",iids,".e",
                  " -o ",config$GCTA_cluster_filename_prefix, "_",iids,".o",' "/PHShome/me704/gcta_1.93.2beta/gcta64 --reml --reml-no-constrain --pheno ',
                  paste0(config$phenotypes_GCTA_prefix,".ranknorm.pheno"),' --mpheno ',iids," --covar ",
                  paste0(config$phenotypes_GCTA_prefix,".covar")," --qcovar ",paste0(config$phenotypes_GCTA_prefix,".qcovar"),
                  " --mgrm-gz ",config$two_matrix_GCTA_filename," --thread-num 8 --out ",paste0(config$phenotypes_GCTA_prefix,"_",iids,'"'))
  writeLines(rrows, fileConn)
  close(fileConn)
  Sys.chmod(config$rerun_GCTA_script_filename,"755")
  
  print(paste0("Done. Please run ",config$rerun_GCTA_script_filename))
}else{
  print("Everything ran OK, nothing to redo")
}
