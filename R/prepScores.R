prepScores <- function(Z, formula, family = gaussian(), SNPInfo=NULL, snpNames = "Name", aggregateBy = "gene", kins = NULL, sparse= TRUE, data=parent.frame(), verbose = FALSE){
	#fit Null model
	if(is.null(SNPInfo)){ 
		warning("No SNP Info file provided: loading the Illumina HumanExome BeadChip. See ?SNPInfo for more details")
		load(paste(find.package("skatMeta"), "data", "SNPInfo.rda",sep = "/"))
		aggregateBy = "SKATgene"
	}

	if(!is.null(kins)){
		n <- dim(kins)[1]
		stopifnot(nrow(kins) == nrow(Z))
		stopifnot(family$family == "gaussian")
		if(sparse){
			kins[kins < 2 * 2^{-6}] <- 0
			kins <- forceSymmetric(kins)
			#oo <- order(getcc(kins))
		} else {
			#oo <- 1:n
		}
		
		#ioo <- match(1:n,oo)
		
		if(is.null(colnames(kins))){
			data$id <- 1:ncol(kins)
		} else {
			data$id <- colnames(kins)
		}
		
		nullmodel <- lmekin(formula=update(formula, '~.+ (1|id)'), data=data, varlist = 2*kins,method="REML")	
		nullmodel$theta <- c(nullmodel$vcoef$id*nullmodel$sigma^2,nullmodel$sigma^2)
	
		SIGMA <- nullmodel$theta[1]*2*kins+nullmodel$theta[2]*Diagonal(n)
		X1 <- model.matrix(lm(formula,data=data))
	
		s2 <- sum(nullmodel$theta)
		Om_i <- solve(SIGMA/s2)
	
		#rotate data:
		res <- as.vector(nullmodel$res)* s2 / nullmodel$theta[2]	
		nullmodel$family$var <- function(x){1}
	} else {
		nullmodel <- glm(formula=formula, family = family, data=data)
		res <- residuals(nullmodel, type = "response")
		X1 <- model.matrix(nullmodel)
		n <- nrow(X1)
	}
 
	env <- environment()
	##check format:
	invisible(check_format_skat(Z, SNPInfo, nullmodel,aggregateBy, snpNames))
	
	##match snps in Z with master list in SNPInfo file 
	mysnps <- colnames(Z)
	
	SNPInfo[,aggregateBy] <- as.character(SNPInfo[,aggregateBy])
	
	SItoZ <- which(colnames(Z) %in% SNPInfo[,snpNames])

	which.snps.Z <- colnames(Z) %in% SNPInfo[,snpNames]
	ZtoSI <- match(SNPInfo[,snpNames], mysnps[which.snps.Z])

	nsnps <- sum(!is.na(ZtoSI)) #sum(which.snps.Z)
	if(nsnps == 0){ 
		stop("no column names in Z match SNP names in the SNP Info file!")
	}
	
	if(verbose){
    	cat("\n Scoring... Progress:\n")
    	pb <- txtProgressBar(min = 0, max = nsnps, style = 3)
    	pb.i <- 0
    }
	
	##fit individual betas/se's
	maf0 <- colMeans(Z,na.rm=TRUE)[which.snps.Z]/2
  maf0[is.nan(maf0)] <- -1
  
	maf <- maf0[ZtoSI]
	names(maf) <- SNPInfo[,snpNames]
	
	scores <- apply(Z[,which.snps.Z, drop = FALSE],2,function(z){
		if(any(is.na(z))){
      if(all(is.na(z))) z <- rep(0,length(z))
      mz <- mean(z, na.rm=TRUE)
			z[is.na(z)] <- mz
		}
        if (verbose){
				assign("pb.i", get("pb.i",env)+1,env)
				if(get("pb.i", env)%%ceiling(nsnps/100) == 0) setTxtProgressBar(get("pb",env),get("pb.i",env))
			}
		sum(res*z)
		})[ZtoSI]
	scores[is.na(scores)] <- 0
	names(scores) <- SNPInfo[,snpNames]
	
	if(verbose) close(pb)
	
	#deal with monomorphic SNPs
	scores[maf == 0] <- 0
	
	#differentiate missing from monomorphic:
	maf[!(SNPInfo[,snpNames] %in% colnames(Z))] <- -1

	#split into genes
	scores 	<- 	split(scores, SNPInfo[,aggregateBy])
	maf 	<- 	split(maf, SNPInfo[,aggregateBy])
	
	##get matrices for projection
	X1 <- sqrt(nullmodel$family$var(nullmodel$fitted))*X1
	if(is.null(kins)){
    AX1 <- with(svd(X1),  v[,d > 0,drop=FALSE]%*%( (1/d[d>0])*t(u[, d > 0,drop=FALSE]))) 
	} else {
	  AX1 <- with(svd(t(X1)%*%Om_i%*%X1),  v[,d > 0,drop=FALSE]%*%( (1/d[d>0])*t(v[, d > 0,drop=FALSE])))%*%t(X1)%*%Om_i
	}
	ngenes <- length(unique(SNPInfo[,aggregateBy]))
	if(verbose){
    	cat("\n Calculating covariance... Progress:\n")
    	pb <- txtProgressBar(min = 0, max = ngenes, style = 3)
    	pb.i <- 0
    }
	##get covariance matrices:
	re <- tapply(SNPInfo[,snpNames], SNPInfo[,aggregateBy],function(snp.names){
		inds <- match(snp.names,colnames(Z))
		mcov <- matrix(0,length(snp.names),length(snp.names))
		if(length(na.omit(inds)) > 0){
			Z0 <- sqrt(nullmodel$family$var(nullmodel$fitted))*as.matrix(Z[,na.omit(inds),drop=FALSE])
			if(any(is.na(Z0))) Z0 <- apply(Z0,2,function(z){
			  if(all(is.na(z))) z <- rep(0,length(z))
       			mz <- mean(z, na.rm=TRUE)
				z[is.na(z)] <- mz
				z
			})
			if(!is.null(kins)){
				mcov[!is.na(inds), !is.na(inds)] <- as.matrix(t(Z0)%*%Om_i%*%Z0 - (t(Z0)%*%Om_i%*%X1)%*%(AX1%*%Z0))
			} else {
				mcov[!is.na(inds), !is.na(inds)] <- crossprod(Z0) - (t(Z0)%*%X1)%*%(AX1%*%Z0)
			}
		}
		rownames(mcov) <- colnames(mcov) <- snp.names
		if(verbose){
				assign("pb.i", get("pb.i",env)+1,env)
				if(get("pb.i", env)%%ceiling(ngenes/100) == 0) setTxtProgressBar(get("pb",env),get("pb.i",env))		  
		}
		return(forceSymmetric(Matrix(mcov,sparse=TRUE)))
	},simplify = FALSE)
	sey = sqrt(var(res)*(nrow(X1)-1)/(nrow(X1)-ncol(X1)) )
	if(family$family == "binomial") sey = 1
	if(!is.null(kins)) 	sey = sqrt(s2)

	##aggregate
	for(k in 1:length(re)){
		re[[k]] <- list("scores" = scores[[k]], "cov" = re[[k]], "n" =n, "maf" = maf[[k]], "sey" = sey ) 
	}
	if(verbose) close(pb)
	
	attr(re,"family") <-  family$family
	class(re) <- "seqMeta"
	return(re)
}


### internal function to get connected components (pedigrees) of a sparse kinship matrix - reordering to get a block diagonal matrix can improve speed dramatically

getcc <- function(M){
	membership <- rep(0,ncol(M))
	clust <- 1
	ii <- 1
	visited <- 1
	tovisit <- NULL
	membership[ii] <- clust
	
	while(any(membership == 0)){
		ne <- which(M[ii,] != 0 & membership != clust)
		if(length(ne) > 0){
			membership[ne] <- clust
			tovisit <- union(tovisit,ne)
		}
		if(length(tovisit) >= 1){
			ii <- tovisit[1]
			tovisit <- tovisit[-1]
		} else {
			clust <- clust+1
			ii <- which(membership == 0)[1]
		}
	}	
	return(membership)
}


prepScoresX <- function(Z, formula, male, family = gaussian(), SNPInfo=NULL, snpNames = "Name", aggregateBy = "gene", kins = NULL, sparse= TRUE, data=parent.frame(), verbose = FALSE){
  #fit Null model
  if(is.null(SNPInfo)){ 
    warning("No SNP Info file provided: loading the Illumina HumanExome BeadChip. See ?SNPInfo for more details")
    load(paste(find.package("skatMeta"), "data", "SNPInfo.rda",sep = "/"))
    aggregateBy = "SKATgene"
  }
  cl <- match.call()
  
  
  if(!is.null(kins)){
    n <- dim(kins)[1]
    stopifnot(nrow(kins) == nrow(Z))
    stopifnot(family$family == "gaussian")
    if(sparse){
      kins[kins < 2 * 2^{-6}] <- 0
      kins <- forceSymmetric(kins)
      #oo <- order(getcc(kins))
    } else {
      #oo <- 1:n
    }
    
    #ioo <- match(1:n,oo)
    
    if(is.null(colnames(kins))){
      data$id <- 1:ncol(kins)
    } else {
      data$id <- colnames(kins)
    }
    
    nullmodel <- lmekin(formula=update(formula, '~.+ (1|id)'), data=data, varlist = 2*kins,method="REML")	
    nullmodel$theta <- c(nullmodel$vcoef$id*nullmodel$sigma^2,nullmodel$sigma^2)
    
    SIGMA <- nullmodel$theta[1]*2*kins+nullmodel$theta[2]*Diagonal(n)
    X1 <- model.matrix(lm(formula,data=data))
    
    s2 <- sum(nullmodel$theta)
    Om_i <- solve(SIGMA/s2)
    
    #rotate data:
    res <- as.vector(nullmodel$res)* s2 / nullmodel$theta[2]	
    nullmodel$family$var <- function(x){1}
  } else {
    nullmodel <- glm(formula=formula, family = family, data=data)
    res <- residuals(nullmodel, type = "response")
    X1 <- model.matrix(nullmodel)
    n <- nrow(X1)
  }

  env <- environment()
  ##check format:
  invisible(check_format_skat(Z, SNPInfo, nullmodel,aggregateBy, snpNames))
  
  male <- eval(cl$male,data)
  if(length(male) != length(res)) stop("`male' not the same length as phenotype")
  if(!all(male %in% c(0,1))) stop("`male' must be coded as 0/1 or T/F")
  male <- as.logical(male)
  
  ##match snps in Z with master list in SNPInfo file 
  mysnps <- colnames(Z)
  
  SNPInfo[,aggregateBy] <- as.character(SNPInfo[,aggregateBy])
  
  SItoZ <- which(colnames(Z) %in% SNPInfo[,snpNames])
  
  which.snps.Z <- colnames(Z) %in% SNPInfo[,snpNames]
  ZtoSI <- match(SNPInfo[,snpNames], mysnps[which.snps.Z])
  
  nsnps <- sum(!is.na(ZtoSI)) #sum(which.snps.Z)
  if(nsnps == 0){ 
    stop("no column names in Z match SNP names in the SNP Info file!")
  }
  
  if(verbose){
    cat("\n Scoring... Progress:\n")
    pb <- txtProgressBar(min = 0, max = nsnps, style = 3)
    pb.i <- 0
  }
  
  ##calculate af, separated by gender
  maf0 <- (colSums(Z[male,which.snps.Z],na.rm=TRUE)/2 + colSums(Z[!male,which.snps.Z],na.rm=TRUE))/
    (colSums(!is.na(Z[male,which.snps.Z])) + 2*colSums(!is.na(Z[!male,which.snps.Z])))

  maf0[is.nan(maf0)] <- -1
  
  maf <- maf0[ZtoSI]
  names(maf) <- SNPInfo[,snpNames]
  
  scores <- apply(Z[,which.snps.Z, drop = FALSE],2,function(z){
    naz <- is.na(z)
    if(any(naz)){
      if(all(naz)){ 
        z <- rep(0,length(z))
      } else {
        mz <- (sum(z[male], na.rm=TRUE)/2+ sum(z[!male],na.rm=TRUE))/
        sum( (2-as.numeric(male))[!naz])
        z[naz] <- 2*mz
      }
    }
    if (verbose){
      assign("pb.i", get("pb.i",env)+1,env)
      if(get("pb.i", env)%%ceiling(nsnps/100) == 0) setTxtProgressBar(get("pb",env),get("pb.i",env))
    }
    sum(res*z)
  })[ZtoSI]
  scores[is.na(scores)] <- 0
  names(scores) <- SNPInfo[,snpNames]
  
  if(verbose) close(pb)
  
  #deal with monomorphic SNPs
  scores[maf == 0] <- 0
  
  #differentiate missing from monomorphic:
  maf[!(SNPInfo[,snpNames] %in% colnames(Z))] <- -1
  
  #split into genes
  scores 	<- 	split(scores, SNPInfo[,aggregateBy])
  maf 	<- 	split(maf, SNPInfo[,aggregateBy])
  
  ##get matrices for projection
  X1 <- sqrt(nullmodel$family$var(nullmodel$fitted))*X1
  if(is.null(kins)){
    AX1 <- with(svd(X1),  v[,d > 0,drop=FALSE]%*%( (1/d[d>0])*t(u[, d > 0,drop=FALSE]))) 
  } else {
    AX1 <- with(svd(t(X1)%*%Om_i%*%X1),  v[,d > 0,drop=FALSE]%*%( (1/d[d>0])*t(v[, d > 0,drop=FALSE])))%*%t(X1)%*%Om_i
  }
  ngenes <- length(unique(SNPInfo[,aggregateBy]))
  if(verbose){
    cat("\n Calculating covariance... Progress:\n")
    pb <- txtProgressBar(min = 0, max = ngenes, style = 3)
    pb.i <- 0
  }
  ##get covariance matrices:
  re <- tapply(SNPInfo[,snpNames], SNPInfo[,aggregateBy],function(snp.names){
    inds <- match(snp.names,colnames(Z))
    mcov <- matrix(0,length(snp.names),length(snp.names))
    if(length(na.omit(inds)) > 0){
      Z0 <- sqrt(nullmodel$family$var(nullmodel$fitted))*as.matrix(Z[,na.omit(inds),drop=FALSE])
      if(any(is.na(Z0))) Z0 <- apply(Z0,2,function(z){
        naz <- is.na(z)
        if(any(naz)){
          if(all(naz)){ 
            z <- rep(0,length(z))
          } else {
            mz <- (sum(z[male], na.rm=TRUE)/2+ sum(z[!male],na.rm=TRUE))/
              sum( (2-as.numeric(male))[!naz])
            z[naz] <- 2*mz
          }
        }
        z
      })
      Z0 <- Z0
      if(!is.null(kins)){
        mcov[!is.na(inds), !is.na(inds)] <- as.matrix(t(Z0)%*%Om_i%*%Z0 - (t(Z0)%*%Om_i%*%X1)%*%(AX1%*%Z0))
      } else {
        mcov[!is.na(inds), !is.na(inds)] <- crossprod(Z0) - (t(Z0)%*%X1)%*%(AX1%*%Z0)
      }
    }
    rownames(mcov) <- colnames(mcov) <- snp.names
    if(verbose){
      assign("pb.i", get("pb.i",env)+1,env)
      if(get("pb.i", env)%%ceiling(ngenes/100) == 0) setTxtProgressBar(get("pb",env),get("pb.i",env))		  
    }
    return(forceSymmetric(Matrix(mcov,sparse=TRUE)))
  },simplify = FALSE)
  sey = sqrt(var(res)*(nrow(X1)-1)/(nrow(X1)-ncol(X1)) )
  if(family$family == "binomial") sey = 1
  if(!is.null(kins)) 	sey = sqrt(s2)
  
  ##aggregate
  for(k in 1:length(re)){
    re[[k]] <- list("scores" = scores[[k]], "cov" = re[[k]], "n" =n, "maf" = maf[[k]], "sey" = sey ) 
  }
  if(verbose) close(pb)
  
  attr(re,"family") <-  family$family
  class(re) <- "seqMeta"
  return(re)
}
