library(hash)
library(config)
library(psych)
library(magrittr) # Yana: added this library, got errors without (needs to be run every time you start R and want to use %>%)
library(dplyr) # Yana: added this library, got errors without (alternatively, this also loads %>%)
library(tidyr)

args = commandArgs(trailingOnly=TRUE)
if(length(args)<1){
  stop("Config file not provided")
}else{
  config <- config::get(file = args[1])
}

#calculate genetic correlation between two traits given One matrix
calc.gencorr <- function( bVC.1, bVC.2, covMatList, XsTXsinv.ClosedForm){
  
  # # # testing:
  # bVC.1 <- VC.1$residuals
  # bVC.2 <- VC.2$residuals
  # XsTXsinv.ClosedForm <- VC.1$XsTXsinv.ClosedForm
  
  #This intentionally allows for duplicates for Bootstrap
  iids<-names(bVC.1)[names(bVC.1) %in% intersect(names(bVC.1),names(bVC.2))]

  kinMat <- covMatList$Kinship
  kinMat <- kinMat[iids,iids]

  houseMat<-covMatList$Household
  houseMat <- houseMat[iids,iids]
  
  agencor = list()
  agencor$K <- NA
  agencor$E <- NA
  agencor$HH <- NA
  
  errorMat = diag(x = 1, nrow = nrow(kinMat), ncol = nrow(kinMat))

  wgts = XsTXsinv.ClosedForm["Kinship",] 
  if("Household" %in% names(wgts)){
    kinMat.adjusted = wgts["error"]*errorMat + wgts["Kinship"]*kinMat + wgts["Household"]*houseMat
  }else{
    kinMat.adjusted = wgts["error"]*errorMat + wgts["Kinship"]*kinMat
  }
  
  kVC.1 <- kinMat.adjusted %*% bVC.1[iids]; names(kVC.1) <- names(bVC.1)
  kVC.2 <- kinMat.adjusted %*% bVC.2[iids]; names(kVC.2) <- names(bVC.2)
  
  arg1 <- t(bVC.1[iids]) %*%  kVC.2
  arg2 <- t(bVC.1[iids]) %*%  kVC.1
  arg3 <- t(bVC.2[iids]) %*%  kVC.2
  agencor$K <- arg1/sqrt(arg2*arg3)
  
  if((arg1==arg2)&(arg2==arg3)){
    agencor$K=abs(agencor$K)
  }
  
  # This is a speedup for bootstrap
  agencor$VC.1 <- bVC.1[iids]
  agencor$VC.2 <- bVC.2[iids]
  agencor$kVC.1 <- kVC.1
  agencor$kVC.2 <- kVC.2
  
  
  wgts = XsTXsinv.ClosedForm["error",] 
  if("Household" %in% names(wgts)){
    kinMat.adjusted = wgts["error"]*errorMat + wgts["Kinship"]*kinMat + wgts["Household"]*houseMat
  }else{
    kinMat.adjusted = wgts["error"]*errorMat + wgts["Kinship"]*kinMat
  }
  
  arg1 <- t(bVC.1[iids]) %*%  kinMat.adjusted %*% bVC.2[iids]
  arg2 <- t(bVC.1[iids]) %*%  kinMat.adjusted %*% bVC.1[iids]
  arg3 <- t(bVC.2[iids]) %*%  kinMat.adjusted %*% bVC.2[iids]
  agencor$E <- arg1/sqrt(arg2*arg3)
  
  if((arg1==arg2)&&(arg2==arg3)){
    agencor$E=abs(agencor$E)
  }
  
  if("Household" %in% names(wgts)){
    wgts = XsTXsinv.ClosedForm["Household",] 
    kinMat.adjusted = wgts["error"]*errorMat + wgts["Kinship"]*kinMat + wgts["Household"]*houseMat
    
    hVC.1 <- kinMat.adjusted %*% bVC.1[iids]; names(hVC.1) <- names(bVC.1)
    hVC.2 <- kinMat.adjusted %*% bVC.2[iids]; names(hVC.2) <- names(bVC.2)
    
    arg1 <- t(bVC.1[iids]) %*%  hVC.2
    arg2 <- t(bVC.1[iids]) %*%  hVC.1
    arg3 <- t(bVC.2[iids]) %*%  hVC.2
    agencor$HH <- arg1/sqrt(arg2*arg3)
    
    if((arg1==arg2)&(arg2==arg3)){
      agencor$HH=abs(agencor$HH)
    }
    
    # This is a speedup for bootstrap
    agencor$hVC.1 <- hVC.1
    agencor$hVC.2 <- hVC.2
  }else{
    agencor$hVC.1 <- NA
    agencor$hVC.2 <- NA
  }
  
  return(agencor)
}

#Do a blocked bootstrap
calc.bootstrap <- function( sVC.1, sVC.2, aVC.1, aVC.2, lookup, Grps, Singles, nrep = 100){
  #debug
  # sVC.1 <- gencor$VC.1
  # sVC.2 <- gencor$VC.2
  # aVC.1 <- gencor$kVC.1
  # aVC.2 <- gencor$kVC.2
  # aVC.1 <- gencor$hVC.1
  # aVC.2 <- gencor$hVC.2
  # lookup <- lookupHH
  # Grps <- GrpsHH
  # Singles <- SinglesHH
  # nrep = 100
  
  #blocked
  gnames <- intersect(names(sVC.1),names(sVC.2))
  tresb <- c()
  if(nrep<1){
    return(tresb)
  }
  tGrps <- Grps[Grps %in% intersect(Grps,gnames)]; Numgrps <- length(tGrps)
  tSingles <- Singles[Singles %in% intersect(Singles,gnames)]
  for(k in 1:nrep){
    allsamp <- c()
    # add correct proportion of groups
    chooGrp <- sample(Grps,Numgrps,replace = TRUE)
    # for each group add all of its members if selected
    for(j in chooGrp){
      tt <- lookup[[j]]
      tt <- intersect(tt,gnames)
      if (length(tt)>0) allsamp <- append(allsamp,tt)
    }
    allsamp <- allsamp[allsamp %in% intersect(allsamp,gnames)]
    # add correct proportion of singles
    allsamp <- append(allsamp,sample(tSingles,length(names(sVC.1))-length(allsamp),replace = TRUE))
    #bring the total number to correct one
    if(length(allsamp)<length(names(sVC.1))){
      allsamp <- append(allsamp,sample(Singles,length(names(sVC.1))-length(allsamp),replace = TRUE))
    }else if(length(allsamp)>length(names(sVC.1))){ 
      allsamp <- sample(allsamp,length(names(sVC.1)),replace = FALSE)
    }
    
    iids <- allsamp
    arg1 <- aVC.1[iids] %*% sVC.2[iids]   
    arg2 <- aVC.1[iids] %*% sVC.1[iids]   
    arg3 <- aVC.2[iids] %*% sVC.2[iids]   
    gencor.K <- arg1/sqrt(arg2*arg3)
    tresb <- append(tresb,gencor.K)
  }
  # mean(tresb,na.rm = TRUE)
  # sd(tresb,na.rm = TRUE)
  # gencor$K
  return(tresb)
}

### calculate variance commponents from Y, W, and covMat.L
calc.vc.resid <- function(Y, W, covMat.L, verbose = FALSE, return.residuals = TRUE, max.iter = 20, eps = 1e-4, diag.eps = 0, vectorize = TRUE, cholesky = TRUE){
  # testing:
  # Y = Y.1; W = W.1
  # covMat.L = covmatrix.list
  # verbose = TRUE; return.residuals = TRUE
  # max.iter = 20; eps = 1e-4; diag.eps = 0.001
  # vectorize = TRUE; cholesky = TRUE
  
  n <- length(Y)
  
  if (verbose){message("NEW: Transforming covariance matrices in preparation for analysis...")}
  
  n.cov.mat <- length(covMat.L)
  names.cov.mat <- names(covMat.L)
  
  if(vectorize){
    # tamar's method: (may contribute to problem with large vectors e.g. 34M elements?)
    covMatVecList <- vector(mode = "list", length = n.cov.mat)
    names(covMatVecList) <- names.cov.mat
    for (i in 1:n.cov.mat){
      # covMatVecList[[i]] <- covMat.L[[i]][upper.tri(covMat.L[[i]], diag = FALSE)] 
      covMat = covMat.L[[i]]
      diag(covMat) = diag.eps # now including small diagonal epsilon element
      covMatVecList[[i]] <- covMat[upper.tri(covMat, diag = TRUE)]
    }
  } else {
    # matthew's method: (keep matrix rather than vectorize to avoid large vector allocation)
    covMatVecList <- covMat.L
    
  }
  # lapply(covMatVecList, sum)
  
  XtX <- matrix(n, nrow = 1 + n.cov.mat, ncol = 1 + n.cov.mat)
  for (i in 1:n.cov.mat){
    for (j in 1:n.cov.mat){
      if(vectorize){
        XtX[i + 1, j + 1] <- sum(covMatVecList[[i]]*covMatVecList[[j]])  + n
      } else {
        XtX[i + 1, j + 1] <- sum(covMatVecList[[i]]*covMatVecList[[j]])
      }
    }
  }
  rownames(XtX) = colnames(XtX) = c("error", names.cov.mat )
  XtXinv <- solve(XtX)
  
  ############# Michael #####################
  # This is the matrix for the closed form solution, will be used downstream
  Xs <- matrix(0,ncol = n.cov.mat+1, nrow = n*n)
  Xs[seq(1,n*n,n+1),1]<-1
  for (i in 1:n.cov.mat){
    Xs[,i+1] <- c(covMat.L[[i]]); 
  }
  
  XsTXs <- t(Xs) %*% Xs
  rownames(XsTXs) = colnames(XsTXs) = c("error", names(covMatList) )
  XsTXsinv <- solve(XsTXs)
  ###########################################
  
  if (verbose) message("Fitting the model...")
  
  ########  fitting the model:
  converged <- FALSE
  n.iter <- 0
  
  mod.lm <- lm(Y ~ -1 + W) # W is design matrix w/intercept
  beta.init <- coef(mod.lm)
  fits <- tcrossprod(W, t(beta.init))
  residM <- as.vector(Y - fits)
  sum.residM.sq <- sum(residM^2)
  # resids.mat <- residM %*% t(residM) # tamar
  resids.mat <- crossprod(t(residM)) # matthew
  
  if(vectorize){
    # tamar's method 
    # epsilon.vec <- resids.mat[upper.tri(resids.mat, diag = FALSE)]
    epsilon.vec <- resids.mat[upper.tri(resids.mat, diag = TRUE)] #to match cov mat with diag.eps
  } else {
    # matthew's method (keep matrix rather than vectorize to avoid large vector allocation)
    epsilon.vec <- resids.mat 
  }
  
  Xerr = NULL
  Xerr <- matrix(sum.residM.sq, nrow = 1 + n.cov.mat, ncol = 1)
  for (i in 1:n.cov.mat){
    Xerr[i + 1,1] <- Xerr[i + 1,1] + sum(covMatVecList[[i]]*epsilon.vec )
  }
  # SEE BELOW FOR CALCULATION CHECKING
  
  vc.est.init <- XtXinv  %*% Xerr
  vc.est.init <- pmax(vc.est.init, 0)       #Michael - this is required for Cholesky to work, but in fact causes error
  rownames(vc.est.init) <- c("error", names.cov.mat )
  
  if (verbose) {
    message("Initial variance component estimates ")
    print(vc.est.init)	
    message("Initial beta estimates ")
    print(beta.init)	
  }
  
  vc.est.cur = vc.est.init
  beta = beta.init
  
  while(!converged){
    n.iter <- n.iter + 1
    if(exists("WtVinvWInv")){
      W_prev<-W
      V_prev<-V
      WtVinvWInv_prev<-WtVinvWInv
    }  
    V <-   diag(rep(vc.est.cur[1], n)) 
    for (i in 1:n.cov.mat){
      V <- V + vc.est.cur[i + 1]*covMat.L[[i]]
    }
    
    possibleError <- tryCatch(
      expr = {
        if(cholesky){
          # Tamar's method (numeric problems with asymetry and error in cholesky?)
          # # cholesky decomposition inverse
          cholV <- chol(V)
          Vinv  <- chol2inv(cholV)
          VinvW <- crossprod(Vinv,W)
          cholWtVinvW <- chol(crossprod(W, VinvW))
          WtVinvWInv  <- chol2inv(cholWtVinvW)
          beta <- crossprod(WtVinvWInv, crossprod(VinvW,Y))
        } else {
          # Matthew's inefficient method
          # solve inverse
          Vinv <- solve(V)
          WtVinvW <- t(W) %*% Vinv %*% W
          WtVinvWInv <- solve( WtVinvW )
          beta <- crossprod(WtVinvWInv, t(W) %*% Vinv %*% Y)
        }    
      },
      error = function(e){ 
        print("Cholesky failed, using best previous values")
        if(exists("WtVinvWInv_prev")){
          W<-W_prev
          V<-V_prev
          WtVinvWInv<-WtVinvWInv_prev
        }  
        converged=TRUE
      }
    )    
    if(inherits(possibleError, "error")) next
    
    if (verbose) cat("iteration count ", n.iter, "current beta ", beta, "\n")
    fits <- tcrossprod(W, t(beta))
    residM <- as.vector(Y - fits)
    sum.residM.sq <- sum(residM^2)
    resids.mat <- residM %*% t(residM)
    
    if(vectorize){
      # tamar
      # epsilon.vec <- resids.mat[upper.tri(resids.mat, diag = FALSE)]
      epsilon.vec <- resids.mat[upper.tri(resids.mat, diag = TRUE)] #to match cov mat with diag.eps
      
    } else {
      # matthew
      epsilon.vec <- resids.mat 
    }
    
    Xerr <- matrix(sum.residM.sq, nrow = 1 + n.cov.mat, ncol = 1)
    for (i in 1:n.cov.mat){
      Xerr[i + 1,1] <- Xerr[i + 1,1] + sum( covMatVecList[[i]]*epsilon.vec )
    }
    
    vc.est.new <- XtXinv %*% Xerr
    vc.est.new <- pmax(vc.est.new, 0)     #Michael, this is not so good
    rownames(vc.est.new) <- c("error", names.cov.mat )
    if (verbose) {cat("iteration count ", n.iter, "current VC estimates:");print(vc.est.new)}
    
    if (max(abs(vc.est.new - vc.est.cur) < eps) | n.iter >= max.iter) converged <- TRUE else{
      vc.est.cur <- vc.est.new
    }
  }
  
  if (verbose) cat("Finished estimation, preparing return values...", "\n")
  ### convenient return values:
  
  #### now to test variance component: Q = vc.est.cur[i],  we need eigenvalues of matrix K, that has the matrix from 
  ### the quadratic form, times the variance of the residuals. 
  ## let My = (y -fits). Then var(My) = M %*% V %*% t(M), which comes to:
  if(exists("WtVinvWInv")){
    V.My <- V  -  W %*% tcrossprod(WtVinvWInv, W)
  }else{
    V.My <- NULL
  }
  
  vc.est.cur <- drop(vc.est.cur)
  prop.var <- vc.est.cur/sum(vc.est.cur)
  names(prop.var) <- names(vc.est.cur)
  
  if(return.residuals){
    residuals = residM
  } else {
    residuals = NULL
  }
  
  return(list(VC = vc.est.cur, propVar = prop.var, var.resids = V.My, XtXinv = XtXinv, XsTXsinv.ClosedForm = XsTXsinv, residuals = residuals, covarResids = mod.lm$residuals))
}

calc.heritability <- function(phenotype.1, phenotype.2, covariates.1, covariates.2, pdata.df, covmatrix.list, 
                              subj.subset = NULL, ID.name = NULL, GRM.name = "Kinship", GRM.only.TF = FALSE, vectorize.TF = TRUE){
  # # # testing:
  # phenotype.1 = pheno.1
  # phenotype.2 = pheno.2
  # pdata.df    = pdata
  # covmatrix.list = covMatList
  # covariates.1 = covariates1
  # covariates.2 = covariates2
  # subj.subset = NULL
  # ID.name = "scanID"
  # GRM.name = "Kinship"
  # GRM.only.TF = FALSE
  # vectorize.TF = TRUE
  # # 
  # # Get Subjects
  # #####
  
  # get pdata ID vector
  if(is.null(ID.name)){
    pIDs = rownames(pdata.df)
  } else {
    pIDs = as.matrix(pdata.df[,ID.name])
  }
  
  # subjects who are complete cases and are in all covmatrix.list entries
  avail.subjs = intersect( pIDs[complete.cases(pdata.df[,c(phenotype.1, phenotype.2, covariates.1, covariates.2)])], Reduce(intersect, lapply(covmatrix.list, rownames ) ) )
  # note: covmatrix.list subjects are not consistent. sapply(covMatList, dim)
  
  # get subject intersection
  if(!is.null(subj.subset)){
    # subjs = intersect( avail.subjs , subj.subset )
    subjs = subj.subset[subj.subset %in% intersect( avail.subjs , subj.subset )]  #keep duplicates
  } else {
    subjs = avail.subjs
  }
  
  # sort subjects either numerically or alphabetically
  ids = sort(subjs)  # sort alphabetically
  
  # n = length(unique(ids))
  n = length(ids)
  
  # Get cov matrices with correct sorted final subject ids
  #####
  for (i in 1:length(covmatrix.list)){
    covmatrix.list[[i]] <- covmatrix.list[[i]][match(ids,rownames(covmatrix.list[[i]])),match(ids,rownames(covmatrix.list[[i]]))]
  }
  sapply(covmatrix.list, dim)
  
  #################################################################################################################
  # Prep data for regression:
  
  if(!all(covariates.1 %in% colnames(pdata.df)) | !all(covariates.2 %in% colnames(pdata.df)) ){
    warning("not all covariates not found in data")
  }
  
  covariates.1.cols = covariates.1[covariates.1 %in% colnames(pdata.df)]
  covariates.2.cols = covariates.2[covariates.2 %in% colnames(pdata.df)]
  
  temp.dat.1 <- as.data.frame(pdata.df[,c(covariates.1.cols, phenotype.1)])
  temp.dat.1 <- temp.dat.1[match(ids, pIDs),]
  temp.dat.1 <- droplevels(temp.dat.1)
  
  # need to verify that the new data still need to be adjusted for all the covariates
  new.covariates.1.cols <- NULL
  for(cv in covariates.1.cols){
    if(length(unique(temp.dat.1[[cv]]))>1)
      new.covariates.1.cols <- c(new.covariates.1.cols,cv)
  }
  if(length(new.covariates.1.cols)<length(covariates.1.cols)){
    print("Dropping some covariates. Old ones : ");print(covariates.1.cols)
    print("New ones : ");print(new.covariates.1.cols)
    temp.dat.1 <- as.data.frame(pdata.df[,c(new.covariates.1.cols, phenotype.1)])
    temp.dat.1 <- temp.dat.1[match(ids, pIDs),]
    temp.dat.1 <- droplevels(temp.dat.1)
  }
  
  model.formula.1 <- as.formula(paste(paste(phenotype.1,"~"), paste(new.covariates.1.cols,collapse="+")))
  W.1             <- model.matrix(model.formula.1, data=temp.dat.1)
  Y.1             <- model.frame(model.formula.1, data=temp.dat.1)[,phenotype.1]
  mod.lm <- lm(Y.1 ~ -1 + W.1) # W is design matrix w/intercept
  if(any(is.na(mod.lm$coefficients))){ # weird case
    badCov <- names(mod.lm$coefficients[is.na(mod.lm$coefficients)])
    new.covariates.1a.cols <- NULL
    for (kk in 1:length(new.covariates.1.cols)) {
      if(length(grep(new.covariates.1.cols[kk],badCov))>0){
        print(paste0(" Identified linead dependency in: ",new.covariates.1.cols[kk]))
      }else{
        new.covariates.1a.cols <- c(new.covariates.1a.cols,new.covariates.1.cols[kk])
      }
    }
    print("New covariates: ");print(new.covariates.1a.cols)
    temp.dat.1 <- as.data.frame(pdata.df[,c(new.covariates.1a.cols, phenotype.1)])
    temp.dat.1 <- temp.dat.1[match(ids, pIDs),]
    temp.dat.1 <- droplevels(temp.dat.1)
    model.formula.1 <- as.formula(paste(paste(phenotype.1,"~"), paste(new.covariates.1a.cols,collapse="+")))
    W.1             <- model.matrix(model.formula.1, data=temp.dat.1)
    Y.1             <- model.frame(model.formula.1, data=temp.dat.1)[,phenotype.1]
  }
  
  temp.dat.2 <- as.data.frame(pdata.df[,c(covariates.2.cols, phenotype.2)])
  temp.dat.2 <- temp.dat.2[match(ids, pIDs),]
  temp.dat.2 <- droplevels(temp.dat.2)
  
  # need to verify that the new data still need to be adjusted for all the covariates
  new.covariates.2.cols <- NULL
  for(cv in covariates.2.cols){
    if(length(unique(temp.dat.2[[cv]]))>1)
      new.covariates.2.cols <- c(new.covariates.2.cols,cv)
  }
  if(length(new.covariates.2.cols)<length(covariates.2.cols)){
    print("Dropping some covariates. Old ones : ");print(covariates.2.cols)
    print("New ones : ");print(new.covariates.2.cols)
    temp.dat.2 <- as.data.frame(pdata.df[,c(new.covariates.2.cols, phenotype.2)])
    temp.dat.2 <- temp.dat.2[match(ids, pIDs),]
    temp.dat.2 <- droplevels(temp.dat.2)
  }
  
  model.formula.2 <- as.formula(paste(paste(phenotype.2,"~"), paste(new.covariates.2.cols,collapse="+")))
  W.2             <- model.matrix(model.formula.2, data=temp.dat.2)
  Y.2             <- model.frame(model.formula.2, data=temp.dat.2)[,phenotype.2]
  mod.lm <- lm(Y.2 ~ -1 + W.2) # W is design matrix w/intercept
  if(any(is.na(mod.lm$coefficients))){ # weird case
    badCov <- names(mod.lm$coefficients[is.na(mod.lm$coefficients)])
    new.covariates.2a.cols <- NULL
    for (kk in 1:length(new.covariates.2.cols)) {
      if(length(grep(new.covariates.2.cols[kk],badCov))>0){
        print(paste0(" Identified linead dependency in: ",new.covariates.2.cols[kk]))
      }else{
        new.covariates.2a.cols <- c(new.covariates.2a.cols,new.covariates.2.cols[kk])
      }
    }
    print("New covariates: ");print(new.covariates.2a.cols)
    temp.dat.2 <- as.data.frame(pdata.df[,c(new.covariates.2a.cols, phenotype.2)])
    temp.dat.2 <- temp.dat.2[match(ids, pIDs),]
    temp.dat.2 <- droplevels(temp.dat.2)
    model.formula.2 <- as.formula(paste(paste(phenotype.2,"~"), paste(new.covariates.2a.cols,collapse="+")))
    W.2             <- model.matrix(model.formula.2, data=temp.dat.2)
    Y.2             <- model.frame(model.formula.2, data=temp.dat.2)[,phenotype.2]
  }
  
  covmatrix.list.k <- list(K = covmatrix.list[[grep(GRM.name, names(covmatrix.list))]])
  
  
  if(GRM.only.TF){
    ##### calculate pooled VC estimates from GRM alone:
    VC.1.k <- calc.vc.resid(Y = Y.1, W = W.1, covMat.L = covmatrix.list.k, return.residuals = TRUE, max.iter = 20, eps = 1e-4, diag.eps = 0, vectorize = vectorize.TF, cholesky = TRUE)
    VC.2.k <- calc.vc.resid(Y = Y.2, W = W.2, covMat.L = covmatrix.list.k, return.residuals = TRUE, max.iter = 20, eps = 1e-4, diag.eps = 0, vectorize = vectorize.TF, cholesky = TRUE)
    
    VC.1 = VC.1.k; VC.2 = VC.2.k
    names(VC.1$residuals) <- rownames(covmatrix.list[[1]]); 
    names(VC.2$residuals) <- rownames(covmatrix.list[[1]]); 
    
  } else {
    tryCatch(
      expr = {
        ##### calculate pooled VC estimates from all covariance matrices:
        VC.1.all <- calc.vc.resid(Y = Y.1, W = W.1, covMat.L = covmatrix.list, return.residuals = TRUE, max.iter = 20, eps = 1e-4, diag.eps = 0, vectorize = vectorize.TF, cholesky = TRUE)
        VC.2.all <- calc.vc.resid(Y = Y.2, W = W.2, covMat.L = covmatrix.list, return.residuals = TRUE, max.iter = 20, eps = 1e-4, diag.eps = 0, vectorize = vectorize.TF, cholesky = TRUE)
        
        VC.1 = VC.1.all; VC.2 = VC.2.all
        names(VC.1$residuals) <- rownames(covmatrix.list[[1]]); 
        names(VC.2$residuals) <- rownames(covmatrix.list[[1]]); 
        
      },
      error = function(e){
        print("Cholesky failed, trying the slow way...")
        ##### calculate pooled VC estimates from all covariance matrices:
        VC.1.all <- calc.vc.resid(Y = Y.1, W = W.1, covMat.L = covmatrix.list, return.residuals = TRUE, max.iter = 20, eps = 1e-4, diag.eps = 0, vectorize = vectorize.TF, cholesky = FALSE)
        VC.2.all <- calc.vc.resid(Y = Y.2, W = W.2, covMat.L = covmatrix.list, return.residuals = TRUE, max.iter = 20, eps = 1e-4, diag.eps = 0, vectorize = vectorize.TF, cholesky = FALSE)
        
        VC.1 = VC.1.all; VC.2 = VC.2.all
        names(VC.1$residuals) <- rownames(covmatrix.list[[1]]); 
        names(VC.2$residuals) <- rownames(covmatrix.list[[1]]); 
      }
    )
  }
  
  return.list <- list(trait.1.her = VC.1$propVar["K"], trait.2.her = VC.2$propVar["K"], 
                      trait.1 = phenotype.1, trait.2 = phenotype.2, 
                      VC.1.list = VC.1, VC.2.list = VC.2,
                      covmatrix.list = covmatrix.list)
  return( return.list )
  
}

fisher_transform <- function(x){
  if (x <= -1) x <- -0.9999
  if (x >= 1) x <- 0.9999
  0.5*log((1+x)/(1-x))
}

fisher_inv_transform <- function(x){
  (exp(2*x) -1)/(exp(2*x) + 1)
}

calc.fisher <- function(gencor, n, conf = 0.05){
  if(is.na(gencor)){
    return(list(low = NA, high = NA))
  }
  if((gencor < -1)|(gencor > 1)){
    return(list(low = NA, high = NA))
  }
  qnorm_conf <- qnorm(conf/2, lower.tail = FALSE )
  z <- fisher_transform(gencor)
  low.z <- z - qnorm_conf/sqrt(n-3)
  high.z <- z + qnorm_conf/sqrt(n-3)
  
  low = fisher_inv_transform(low.z); high = fisher_inv_transform(high.z)
  return(list(low=low,high=high))
}

calc.fisher.CI <- function(zmean, zsd, conf = 0.05){
  qnorm_conf <- qnorm(conf/2, sd = zsd ,lower.tail = FALSE)
  low.z <- zmean - qnorm_conf
  high.z <- zmean + qnorm_conf
  
  low = fisher_inv_transform(low.z); high = fisher_inv_transform(high.z)
  return(list(low=low,high=high))
}

getCurrentFileLocation <-  function()
{
  this_file <- commandArgs() %>% 
    tibble::enframe(name = NULL) %>%
    tidyr::separate(col=value, into=c("key", "value"), sep="=", fill='right') %>%
    dplyr::filter(key == "--file") %>%
    dplyr::pull(value)
  if (length(this_file)==0)
  {
    this_file <- rstudioapi::getSourceEditorContext()$path
  }
  return(this_file)
}

################################### Begin ##############################################
PARTIAL <- config$PARTIAL

outfileprefix <- config$outfileprefix
workdirname <- config$workdirname

block_file1 <- config$kinship_matrix_filename_blocked_lookup
block_file2 <- config$household_matrix_filename_blocked_lookup

if (PARTIAL){
  pairfile <- config$pairfile
  grpfile <- config$grpfile
  load(grpfile)
}

print("Loading Phenotypes...")
phenotypes <- readRDS(config$phenotypes_filename)

covariates0 = config$factor_covariates
covariates1 <- config$numeric_covariates
covariates = c(covariates0,covariates1)

if((length(config$phenotypes_fieldnames)==1) & (length(grep('-',config$phenotypes_fieldnames))==1)){   #10-100
  sst <- as.numeric(strsplit(config$phenotypes_fieldnames,"-")[[1]][1])
  sse <- as.numeric(strsplit(config$phenotypes_fieldnames,"-")[[1]][2])
  IBD_order_TOPMed_NormK <- colnames(phenotypes)[sst:sse]
}else{
  IBD_order_TOPMed_NormK <- config$phenotypes_fieldnames
}

rownames(phenotypes) <- phenotypes[[config$idfield_name]]
phenotypes <- phenotypes[,c(config$idfield_name,covariates,IBD_order_TOPMed_NormK)]
fstart <- length(covariates)+2
IDFIELD <- config$idfield_name

print("Loading relatedness matrices...")
mat_phi <- readRDS(config$kinship_matrix_filename)
cn <- colnames(mat_phi);cn <- cn[!is.na(cn)]
cn <- cn[cn %in% phenotypes[[IDFIELD]]]

if(length(config$household_matrix_filename)>0){
  mat_hh <- readRDS(config$household_matrix_filename)
  cn <- cn[cn %in% colnames(mat_hh)]
}


rownames(phenotypes) <- phenotypes[[IDFIELD]]
pdata <- phenotypes[cn,]

mat_phi <- mat_phi[cn,cn]
if(length(config$household_matrix_filename)>0){
  mat_hh <- mat_hh[cn,cn]
}

if(length(config$household_matrix_filename)>0){
  covMatList<-list(as.matrix(mat_phi),as.matrix(mat_hh))
  names(covMatList)<-c("Kinship","Household")
}else{
  covMatList<-list(as.matrix(mat_phi))
  names(covMatList)<-c("Kinship")
}


ONEBATCH <- config$ONEBATCH
NREP <- config$BOOTSTRAP_REPEATS

#adjust for Weight_Norm and Age and Sex 
print("Adjusting raw data for Covariates...")
phenotypesN <- colnames(pdata)
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

phenotypes<-colnames(pdata)[fstart:length(colnames(pdata))]
tpname <- paste0(workdirname,outfileprefix,"_Phenotypes",length(phenotypes),"_PerBatch",ONEBATCH,"_accessory.RData")
if(file.exists(tpname)){
  print("Loading accessory data")
  load(tpname)
}else{
  print(paste0("Generating and saving accessory data: ",tpname))
  accD <- list()
  if (PARTIAL){   
    tpairs <- read.table(file = pairfile,sep="\t",header = TRUE)
    tpairs <- data.frame(a=tpairs$row,b=tpairs$col)
  }else{
    if((length(config$phenotypes_group1)==1) & (length(grep('-',config$phenotypes_group1))==1)){   #10-100
      sst <- as.numeric(strsplit(config$phenotypes_group1,"-")[[1]][1])
      sse <- as.numeric(strsplit(config$phenotypes_group1,"-")[[1]][2])
      flds1 <- colnames(pdata)[sst:sse]
    }else{
      flds1 <- config$phenotypes_group1
    }
    
    if((length(config$phenotypes_group2)==1) & (length(grep('-',config$phenotypes_group2))==1)){   #10-100
      sst <- as.numeric(strsplit(config$phenotypes_group2,"-")[[1]][1])
      sse <- as.numeric(strsplit(config$phenotypes_group2,"-")[[1]][2])
      flds2 <- colnames(pdata)[sst:sse]
    }else{
      flds2 <- config$phenotypes_group2
    }
    
    a <- match(flds1,phenotypes)
    b <- match(flds2,phenotypes)
    if((length(a)!=length(flds1))|(length(b)!=length(flds2))){
      stop("Unrecognized field names in config:phenotypes_group1 or config:phenotypes_group2")
    }
    tpairs <- expand.grid(a = a, b = b)
    tpairs.sort = t(apply(tpairs, 1, sort))
    tpairs <- tpairs[!duplicated(tpairs.sort),]
  }
  
  accD$tpairs <- tpairs
  save(accD,file=tpname)
}

print(paste0("Running ",Sys.getenv("LSB_JOBINDEX")))
args<-c(Sys.getenv("LSB_JOBINDEX"))
# print(paste0(" args ",args," - ", as.numeric(args[1])))
# quit(save = "no")
if ((args=="")|(args==0)){ #Generate runScript
  #final run script
  fileConn<-file(config$run_Gencor_script_filename)
  rrows <- c("#!/bin/sh",
             paste0("#BSUB -o ",config$workdirname,config$outfileprefix,"_%J_%I.out"),
             paste0("#BSUB -e ",config$workdirname,config$outfileprefix,"_%J_%I.err"),
             paste0("#BSUB -q ",config$LSF_QUEUE_NAME),
             paste0('#BSUB -R "rusage[mem=',config$LSF_MEMORY,']"'),
             paste0('#BSUB -J ',config$outfileprefix,'[1-',ceiling(length(accD$tpairs$a) / ONEBATCH),']'),
             "module load R/testversions/4.2.0",
             paste0("Rscript ",getCurrentFileLocation()," ",config$this_config_filename))
  writeLines(rrows, fileConn)
  close(fileConn)
  Sys.chmod(config$run_Gencor_script_filename,"755")
  
  fileConn<-file(gsub(".sh",".standalone.sh",config$run_Gencor_script_filename))
  iids <- 1:ceiling(length(accD$tpairs$a) / ONEBATCH)
  rrows <- c("#!/bin/sh",
             paste0("export LSB_JOBINDEX=",iids,"; Rscript ",getCurrentFileLocation()," ",config$this_config_filename))
  writeLines(rrows, fileConn)
  close(fileConn)
  Sys.chmod(gsub(".sh",".standalone.sh",config$run_Gencor_script_filename),"755")

  print(paste0('Generated runfile for cluster (',ceiling(length(accD$tpairs$a) / ONEBATCH),' jobs). Please run from command line via: bsub < ',config$run_Gencor_script_filename))
  
  quit(save = "no")
  # args<-c(5)
}

N <- as.numeric(args[1])
sstart<-(N-1)*ONEBATCH+1
ssend<-N*ONEBATCH

if (sstart>length(accD$tpairs$a)){
  stop(paste("Index too high - ", N, " exiting..."))
}

if(ssend>length(accD$tpairs$a)){
  ssend<-length(accD$tpairs$a)
}

print(paste0("Calculating batch ",N," of size ",ONEBATCH," (out of ",ceiling(length(accD$tpairs$a) / ONEBATCH),"): samples ",sstart,"-",ssend," out of ",length(accD$tpairs$a)))

toutfname=paste0(outfileprefix,"_",ceiling(length(accD$tpairs$a) / ONEBATCH),"batches_",ONEBATCH,"perbatch_batch_",N,".txt")
if (length(list.files(path=workdirname,pattern=toutfname))>0){
  stop("Already calculated this batch, skipping...")
}


##################################### Load and calc data for blocked bootstrap

if((length(block_file2)>0) & (NREP>0)){
  load(block_file1)
  lookupK <- lookup
  NumsinglesK=0;NumgrpsK=0
  SinglesK<-c();GrpsK<-c()
  g <- hash();
  kl <- intersect(keys(lookupK),colnames(mat_phi))
  for( k in kl){
    if(length(lookupK[[k]])<2){
      NumsinglesK <- NumsinglesK + 1
      SinglesK <- append(SinglesK,k)
    }else{
      if(!(paste0(lookupK[[k]],collapse = "-") %in% keys(g))){
        NumgrpsK <- NumgrpsK + 1
        GrpsK <- append(GrpsK,k)
        g[[paste(lookupK[[k]],collapse = "-")]] = k
      }
    }
  }
}

if((length(block_file2)>0) & (NREP>0)){
  load(block_file2)
  lookupHH <- lookup
  NumsinglesHH=0;NumgrpsHH=0
  SinglesHH<-c();GrpsHH<-c()
  g <- hash();
  kl <- intersect(keys(lookupHH),colnames(mat_hh))
  for( k in kl){
    if(length(lookupHH[[k]])<2){
      NumsinglesHH <- NumsinglesHH + 1
      SinglesHH <- append(SinglesHH,k)
    }else{
      if(!(paste0(lookupHH[[k]],collapse = "-") %in% keys(g))){
        NumgrpsHH <- NumgrpsHH + 1
        GrpsHH <- append(GrpsHH,k)
        g[[paste(lookupHH[[k]],collapse = "-")]] = k
      }
    }
  }
}
###################################################################################

amodes <- list("GenCorr", "ResidCorr","NormGenCorr", "HHCorr","NormHHCorr")
res=matrix(NA, nrow = 40, ncol = (ssend-sstart+1))
rownames(res)<-c(unlist(amodes),
                 "KinshipSigma1","KinshipSigma2","ErrorSigma1","ErrorSigma2","HouseholdSigma1","HouseholdSigma2","n_people","Spearman","Spearman_pval","modelPheno","modelPheno_pval",
                 "Fisher_k_CI_Low","Fisher_k_CI_High","Fisher_k_pval","Fisher_HH_CI_Low","Fisher_HH_CI_High","Fisher_HH_pval",
                 "Bootstrap_k_Nrep","GenCorrSD_CI_low","GenCorrSD_CI_high","GenCorr_F_pval","GenCorrQ_CI_low","GenCorrQ_CI_high","GenCorrQ_sig","Bootstrap_k_vals",
                 "Bootstrap_HH_Nrep","HHCorrSD_CI_low","HHCorrSD_CI_high","HHCorr_F_pval","HHCorrQ_CI_low","HHCorrQ_CI_high","HHCorrQ_sig","Bootstrap_HH_vals",
                 "K_neff","HH_neff")
resnames <- rep("", ssend-sstart+1)
pIDs = as.matrix(pdata[,IDFIELD])
for (i in sstart:ssend){  # i = sstart+1
  pheno.1 = phenotypes[accD$tpairs$a[i]]   # accD$tpairs$a[i] <- 1
  pheno.2 = phenotypes[accD$tpairs$b[i]]   # accD$tpairs$b[i] <- 2
  avail.subjs = intersect( pIDs[complete.cases(pdata[,c(pheno.1, pheno.2, covariates, covariates)])], Reduce(intersect, lapply(covMatList, rownames ) ) )
  print(paste0("Calculating VCs for ",pheno.1,"-",pheno.2," (",i-sstart+1," out of ",ssend-sstart+1,") based on ",length(avail.subjs), " samples"))
  
  if(length(avail.subjs)<100){
    print(paste0("Skipping based on too few samples - ",length(na.omit(pdata[[pheno.1]]))," and ",length(na.omit(pdata[[pheno.2]]))," samples: "))
    resnames[i-sstart+1]<-paste(pheno.1,pheno.2,sep=";")
    next
  }
  ######################### partial Genecor calc ####################################
  if(PARTIAL){
    covariates1 <- unlist(Gencor_grps[[pheno.1]]); names(covariates1) <- NULL; covariates1 <- covariates1[!covariates1 %in% c(pheno.1,pheno.2)]
    covariates2 <- unlist(Gencor_grps[[pheno.2]]); names(covariates2) <- NULL; covariates2 <- covariates2[!covariates2 %in% c(pheno.1,pheno.2)]
    covariates1 <- c(covariates,covariates1);covariates2 <- c(covariates,covariates2);
  }else{
    covariates1 <- covariates; covariates2 <- covariates;
  }
  VCs = calc.heritability(
    phenotype.1 = pheno.1,
    phenotype.2 = pheno.2,
    covariates.1 = covariates1,
    covariates.2 = covariates2,
    pdata.df = pdata,
    covmatrix.list = covMatList,
    subj.subset = NULL,
    ID.name =  IDFIELD,
    GRM.name = "Kinship",
    GRM.only.TF = FALSE )
  
  VCs$VC.1.list$VC[VCs$VC.1.list$VC==0]=0.0001
  VCs$VC.2.list$VC[VCs$VC.2.list$VC==0]=0.0001
  
  print(paste0("Sigmas based on ",length(VCs$VC.1.list$residuals)," samples: "))
  print(VCs$VC.1.list$VC)
  print(VCs$VC.2.list$VC)
  
  VC.1 = VCs$VC.1.list
  VC.2 = VCs$VC.2.list
  ######################################################################
  gencor=calc.gencorr(VC.1$residuals, VC.2$residuals, covMatList, VC.1$XsTXsinv.ClosedForm)
  res["KinshipSigma1",i-sstart+1] <- VC.1$VC["Kinship"]
  res["KinshipSigma2",i-sstart+1] <- VC.2$VC["Kinship"]
  if(length(config$household_matrix_filename)>0){
    res["HouseholdSigma1",i-sstart+1] <- VC.1$VC["Household"]
    res["HouseholdSigma2",i-sstart+1] <- VC.2$VC["Household"]
  }
  res["ErrorSigma1",i-sstart+1] <- VC.1$VC["error"]
  res["ErrorSigma2",i-sstart+1] <- VC.2$VC["error"]
  res["n_people",i-sstart+1] <- length(avail.subjs)
  res["GenCorr",i-sstart+1] <- gencor$K
  res["ResidCorr",i-sstart+1] <- gencor$E
  if(length(config$household_matrix_filename)>0){
      res["HHCorr",i-sstart+1] <- gencor$HH
  }
  res["Spearman",i-sstart+1] <- cor.test(VC.1$residuals, VC.2$residuals)$estimate
  res["Spearman_pval",i-sstart+1] <- cor.test(VC.1$residuals, VC.2$residuals)$p.value
  resnames[i-sstart+1]<-paste(pheno.1,pheno.2,sep=";")
  
  ####
  s1=VC.1$VC["Kinship"];s2=VC.2$VC["Kinship"]
  if(length(config$household_matrix_filename)>0){
    h1=VC.1$VC["Household"];h2=VC.2$VC["Household"]
  }else{
    h1=0;h2=0;
  }
  e1=VC.1$VC["error"];e2=VC.2$VC["error"]
  if(length(config$household_matrix_filename)>0){
    res["modelPheno",i-sstart+1] <- e1*e2*gencor$E + s1*s2*gencor$K + h1*h2*gencor$HH
  }else{
    res["modelPheno",i-sstart+1] <- e1*e2*gencor$E + s1*s2*gencor$K
  }

  n = length(avail.subjs); r= as.numeric(res["modelPheno",i-sstart+1])
  if(is.na(r)){
    p <- NA
  }else{
    t <- (r*sqrt(n-2))/sqrt(1-r^2)
    p <- 2*(1 - pt(abs(t),(n-2)))
  }
  res["modelPheno_pval",i-sstart+1] <- p
  
  
  normgencor=(sqrt(s1)*sqrt(s2)*gencor$K) / sqrt((s1+h1+e1)*(s2+h2+e2))
  res["NormGenCorr",i-sstart+1]=normgencor
  
  kinMat.zero.diag <- mat_phi[avail.subjs,avail.subjs]; diag(kinMat.zero.diag) <- 0
  kinMat.adj_n <- sqrt(sum(diag(tcrossprod(kinMat.zero.diag))))
  
  k_CI <- unlist(calc.fisher(gencor$K,kinMat.adj_n,conf = 0.05));
  fisher.vals <- r.test(n = kinMat.adj_n+3, r12 = gencor$K)
  
  
  res["Fisher_k_CI_Low",i-sstart+1] <- k_CI["low"]
  res["Fisher_k_CI_High",i-sstart+1] <- k_CI["high"]
  res["Fisher_k_pval",i-sstart+1] <- fisher.vals$p
  
  print(paste(pheno.1,"-",pheno.2,": K - ",gencor$K," (Fisher_CI: ",k_CI["low"],"-",k_CI["high"],";pval - ",
              res["Fisher_k_pval",i-sstart+1]))
  
  ####
  if(length(config$household_matrix_filename)>0){
    normgencor=(sqrt(h1)*sqrt(h2)*gencor$HH) / sqrt((s1+h1+e1)*(s2+h2+e2))
    res["NormHHCorr",i-sstart+1]=normgencor
    
    hhMat.zero.diag <- mat_hh[avail.subjs,avail.subjs]; diag(hhMat.zero.diag) <- 0
    hhMat.adj_n <- sqrt(sum(diag(tcrossprod(hhMat.zero.diag))))
    
    hh_CI <- unlist(calc.fisher(gencor$HH,hhMat.adj_n,conf = 0.05));
    fisher.vals <- r.test(n = hhMat.adj_n+3, r12 = gencor$HH)
  
    res["Fisher_HH_CI_Low",i-sstart+1] <- hh_CI["low"]
    res["Fisher_HH_CI_High",i-sstart+1] <- hh_CI["high"]
    res["Fisher_HH_pval",i-sstart+1] <- fisher.vals$p
    
    print(paste(pheno.1,"-",pheno.2,": HH - ",gencor$HH," (Fisher_CI: ",hh_CI["low"],"-",hh_CI["high"],";pval - ",
                res["Fisher_HH_pval",i-sstart+1]))
  }else{
    hhMat.adj_n <- NA
  }
  
  res["K_neff",i-sstart+1] <- kinMat.adj_n
  res["HH_neff",i-sstart+1] <- hhMat.adj_n
  ######################### Bootstrap VC calc ####################################
  if((!is.na(gencor$K)) & (NREP>0)){
    print(paste0("Bootstrapping K ",NREP," times..."))
    tresK <- calc.bootstrap(gencor$VC.1,gencor$VC.2,gencor$kVC.1,gencor$kVC.2,lookupK, GrpsK, SinglesK, nrep = NREP)
    
    s <- as.numeric(tresK)
    s[is.nan(s)] <- NA
    s[s == -100] <- NA
    s[s < -1] <- -1
    s[s > 1] <- 1
    
    s2 <- s[!is.na(s)]
    
    z <- sapply(s2, fisher_transform)
    mu <-fisher_transform(gencor$K)
    zmean <- mean(z);zsd <- sd(z); # mean(s2);sd(s2)
    lh <- calc.fisher.CI(zmean,zsd)
    p.value <- pchisq((zmean/zsd)^2, df =1, lower.tail = FALSE)
    
    res["Bootstrap_k_Nrep",i-sstart+1] <- NREP
    res["GenCorrSD_CI_low",i-sstart+1] <- lh$low
    res["GenCorrSD_CI_high",i-sstart+1] <- lh$high
    res["GenCorr_F_pval",i-sstart+1] <- p.value
    res["GenCorrQ_CI_low",i-sstart+1] <- quantile(s2,probs = c(0.025,0.975),na.rm = T)[1]
    res["GenCorrQ_CI_high",i-sstart+1] <- quantile(s2,probs = c(0.025,0.975),na.rm = T)[2]
    res["GenCorrQ_sig",i-sstart+1] <- as.logical(!between(0, res["GenCorrQ_CI_low",i-sstart+1],res["GenCorrQ_CI_high",i-sstart+1]))
    res["Bootstrap_k_vals",i-sstart+1] <- paste0(round(tresK*1000,0)/1000,collapse = ";")
    
    print(paste("K: ",gencor$K," [",res["GenCorrQ_CI_low",i-sstart+1],",",
                res["GenCorrQ_CI_high",i-sstart+1],"] 95% sig - ",
                res["GenCorrQ_sig",i-sstart+1]))
  
  }
  
  if((!is.na(gencor$HH)) & (NREP>0)){
    ####
    print(paste0("Bootstrapping HH ",NREP," times..."))
    tresHH <- calc.bootstrap(gencor$VC.1,gencor$VC.2,gencor$hVC.1,gencor$hVC.2,lookupHH, GrpsHH, SinglesHH, nrep = NREP)
    
    s <- as.numeric(tresHH)
    s[is.nan(s)] <- NA
    s[s == -100] <- NA
    s[s < -1] <- -1
    s[s > 1] <- 1
    
    s2 <- s[!is.na(s)]
    
    z <- sapply(s2, fisher_transform)
    mu <-fisher_transform(gencor$HH)
    zmean <- mean(z);zsd <- sd(z); # mean(s2);sd(s2)
    lh <- calc.fisher.CI(zmean,zsd)
    p.value <- pchisq((zmean/zsd)^2, df =1, lower.tail = FALSE)
    
    res["Bootstrap_HH_Nrep",i-sstart+1] <- NREP
    res["HHCorrSD_CI_low",i-sstart+1] <- lh$low
    res["HHCorrSD_CI_high",i-sstart+1] <- lh$high
    res["HHCorr_F_pval",i-sstart+1] <- p.value
    res["HHCorrQ_CI_low",i-sstart+1] <- quantile(s2,probs = c(0.025,0.975),na.rm = T)[1]
    res["HHCorrQ_CI_high",i-sstart+1] <- quantile(s2,probs = c(0.025,0.975),na.rm = T)[2]
    res["HHCorrQ_sig",i-sstart+1] <- as.logical(!between(0, res["HHCorrQ_CI_low",i-sstart+1],res["HHCorrQ_CI_high",i-sstart+1]))
    res["Bootstrap_HH_vals",i-sstart+1] <- paste0(round(tresHH*1000,0)/1000,collapse = ";")
    
    print(paste("HH: ",gencor$HH," [",res["HHCorrQ_CI_low",i-sstart+1],",",
                res["HHCorrQ_CI_high",i-sstart+1],"] 95% sig - ",
                res["HHCorrQ_sig",i-sstart+1]))
  }
}

colnames(res)<-resnames
print("Saving the resulting file...")
write.table(res, file = paste0(workdirname,toutfname), sep = "\t", col.names=NA)

print("Done.")
