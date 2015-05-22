# prepSNPInfo
#' @export
prepSNPInfo <- function(snpinfo, snpNames, aggregateBy, wt1=NULL, wt2=NULL) {
  
  dups <- duplicated(snpinfo[ , c(aggregateBy, snpNames)])
  snp_na <- is.na(snpinfo[ , snpNames])
  gene_na <- is.na(snpinfo[ , aggregateBy])
  
  idx <- which(!(dups | snp_na | gene_na))
  
  cols <- c(aggregateBy, snpNames)
  if (is.character(wt1) && length(wt1) > 0L) {
    cols <- c(cols, wt1)
  }
  if (is.character(wt2) && length(wt2) > 0L) {
    cols <- c(cols, wt2)
  }
  
  if (length(cols) > 2L) {
    si_unique <- unique(snpinfo[!(snp_na | gene_na) , cols])
    for (i in 3:length(cols)) {
      if (any(is.na(si_unique[ , i]))) {
        stop("Cannot have missing weights.  Please impute and re-run")
      }
    }           
    if (nrow(si_unique) != length(idx)) {
      stop("Non-unique weight for snp-gene pair")
    }    
  } 
  
  nmiss <- sum(snp_na | gene_na)
  if (nmiss > 0L) {
    warning("Removed ", nmiss, " missing snps and/or genes from SNPinfo file") 
  } 
  
  ndups <- sum(dups)
  if (ndups > 0L) {
    warning("Removed ", ndups," duplicate snp-gene pairs from SNPinfo file") 
  }   
  
  snpinfo[idx, cols]
}

# prepGenotype
prepGenotype <- function(Z) {

  if (!is.matrix(Z)) {
    Z <- as.matrix(Z)
  }
  
  if (!is.numeric(Z)) {
    stop("Genotypes must be numeric! Alleles should be coded 0/1/2")
  }

  if(min(Z, na.rm=TRUE) < 0 | max(Z, na.rm=TRUE) > 2) {
    warning("Genotype possibly not dosage matrix! Alleles should be coded 0/1/2")
  }
  
  if(anyNA(Z)) {
    warning("Some missing genotypes - will be imputed to average dose")
  }
  
  Z 
}

# inpute to mean
impute_to_mean <- function(Z) {
  
  if (!is.matrix(Z)) {
    Z <- as.matrix(Z)
  }
  
  if (anyNA(Z)) {
    MZ <- colMeans(Z, na.rm=TRUE)
    allNA <- which(is.nan(MZ))
    Z[ , allNA] <- 0
    MZ[allNA] <- 0
    ISNAZ <- is.na(Z) 
    idx <- which(ISNAZ)
    Z[idx] <- MZ[((idx-1)%/%nrow(Z))+1]
    
    Z
  } else {
    Z
  }
}

# prepPhenotype
create_model <- function(formula, type="continuous", kins=NULL, sparse=TRUE, data=parent.frame()) {
  
  fam <- if (type == "continuous") {
    "gaussian"    
  } else if (type == "binary") {
    "binomial"
  }
    
  if(!is.null(kins)){
    if (type != "continuous") {
      stop("Family data is currently only supported for continuous outcomes.")
    } 
    if(sparse){
      kins[kins < 2^{-5}] <- 0
      kins <- forceSymmetric(kins)
    }
    data$.id <- if(is.null(colnames(kins))){
      1:ncol(kins)
    } else {
      colnames(kins)
    }    
    
    nullmodel <- lmekin(formula=update(formula, '~.+ (1|.id)'), data=data, varlist = 2*kins,method="REML")  
    
    nullmodel$theta <- c(nullmodel$vcoef$id*nullmodel$sigma^2,nullmodel$sigma^2)   
    SIGMA <- nullmodel$theta[1] * 2 * kins + nullmodel$theta[2] * Diagonal(n)   
    s2 <- sum(nullmodel$theta)
    
    #rotate data:
    sef <- sqrt(nullmodel$family$var(nullmodel$fitted))
    X1 <- sef*model.matrix(lm(formula,data=data)) 
    nullmodel$family$var <- function(x){1}
    res <- as.vector(nullmodel$res)* s2 / nullmodel$theta[2]  
    Om_i <- solve(SIGMA/s2)
    # optimize calculations
    tX1_Om_i <- crossprod(X1, Om_i)
    AX1 <- with(svd(tX1_Om_i%*%X1),  v[,d > 0,drop=FALSE]%*%( (1/d[d>0])*t(v[, d > 0,drop=FALSE])))%*%tX1_Om_i
    #    AX1 <- with(svd(t(X1)%*%Om_i%*%X1),  v[,d > 0,drop=FALSE]%*%( (1/d[d>0])*t(v[, d > 0,drop=FALSE])))%*%t(X1)%*%Om_i    
    
    check_dropped_subjects(res, formula)
    list(res=res, family="gaussian", n=nrow(X1), sey=sqrt(s2), sef=sef, X1=X1, AX1=AX1, Om_i=Om_i)
  } else {
    nullmodel <- glm(formula=formula, family=fam, data=data)
    res <- residuals(nullmodel, type = "response")  
    check_dropped_subjects(res, formula) 
    sef <- sqrt(nullmodel$family$var(nullmodel$fitted))
    X1 <- sef*model.matrix(nullmodel) 
    AX1 <- with(svd(X1),  v[,d > 0,drop=FALSE]%*%( (1/d[d>0])*t(u[, d > 0,drop=FALSE])))     
    sey <- if (type == "continuous") {
      sqrt(var(res)*(nrow(X1) - 1)/(nrow(X1) - ncol(X1)) )
    } else if (type == "binary") {
      1
    }    
    list(res=res, family=fam, n=nrow(X1), sey=sey, sef=sef, X1=X1, AX1=AX1)
  }  
}

check_dropped_subjects <- function(res, formula) {
  if (!is.null(na.action(res))) { 
    stop(paste0("Some observations in '", 
                capture.output(print(formula)), 
                "' are missing...\n Complete data in the null model is required. Please remove, and subset genotypes accordingly"))
  }
  invisible(NULL)  
}

# check_format_skat
check_format_skat <- function(Z, SNPInfo, data, snpNames, kins){
  if(nrow(data) != nrow(Z)) {
    stop("Number of genotypes is not equal to number of phenotypes!")
  }
  
  if(!is.null(kins)) {
    if(nrow(kins) != nrow(Z)) {
      stop("Number of genotype subjects is not equal to the number in the kinship matrix!")
    }
  }

  snps <- intersect(colnames(Z), SNPInfo[, snpNames])
  nsnps <- length(snps)
  if(nsnps == 0L) {
    stop("Column names of Z must correspond to 'snpNames' field in SNPInfo file!")
  } 
  if(nsnps < ncol(Z)) {
    warning(paste(ncol(Z) - nsps, "snps are not in SNPInfo file!"))
  } 
  
  invisible(NULL)
}

calculate_maf <- function(Z) {
  
  MeanZ <- colMeans(Z, na.rm=TRUE)
  maf <- MeanZ/2

  #differentiate all missing from monomorphic
  maf[which(is.nan(MeanZ))] <- -1
  
  maf
}

calculate_scores <- function(Z, m) {
  
  scores <- colSums(m$res*Z)
  scores[is.na(scores)] <- 0
  names(scores) <- colnames(Z)
  
  scores  
}

calculate_cov <- function(Z, m, SNPInfo, snpNames, aggregateBy, kins) {
  ##get matrices for projection
  X1 <- m$X1
  AX1 <- m$AX1
  
  ##get covariance matrices:
  re <- tapply(SNPInfo[, snpNames], SNPInfo[, aggregateBy],function(snp.names){
    inds <- intersect(snp.names, colnames(Z))
    if(length(inds) > 0L) {
      mcov <- matrix(0,length(snp.names),length(snp.names), dimnames=list(snp.names, snp.names))
      Z0 <- m$sef*Z[, inds, drop=FALSE]
      Z0 <- impute_to_mean(Z0)     # not sure how Z0 can have missing values
      if(!is.null(kins)){
        tZ0_Omi <- crossprod(Z0%*%m$Om_i)
        mcov[inds, inds] <- as.matrix(tZ0_Omi%*%Z0 - (tZ0_Omi%*%X1)%*%(AX1%*%Z0))
      } else {
        mcov[inds, inds] <- crossprod(Z0) - crossprod(Z0,X1)%*%(AX1%*%Z0)
      }
      forceSymmetric(Matrix(mcov,sparse=TRUE))
    } else{
      Matrix(0, nrow=length(snp.names), ncol=length(snp.names), dimnames=list(snp.names, snp.names), sparse=TRUE)
    }
  },simplify = FALSE)
  
  re
}



fill_values <- function(x, r) { 
  cmn <- intersect(names(x), names(r)) 
  idx <- which(names(x) %in% names(r[cmn]))
  x[idx]<- r[names(x)[idx]]
  x
}


create_seqMeta <- function(re, scores, maf , m, SNPInfo, snpNames, aggregateBy) {
  
  maf_si <- rep_len(-1, nrow(SNPInfo))
  names(maf_si) <- SNPInfo[ , snpNames]  
  maf_si <- fill_values(maf_si, maf)
  maf_split   <-   split(maf_si, SNPInfo[ , aggregateBy])

  scores_si <- rep_len(0, nrow(SNPInfo))
  names(scores_si) <- SNPInfo[ , snpNames]  
  scores_si <- fill_values(scores_si, scores)
  scores_split   <-   split(scores_si, SNPInfo[ , aggregateBy])
  
  ##aggregate
  for(k in 1:length(re)){
    re[[k]] <- list("scores"=scores_split[[k]], "cov"=re[[k]], "n"=m$n, "maf"=maf_split[[k]], "sey"=m$sey ) 
  }
  
  attr(re,"family") <-  m$family
  class(re) <- "seqMeta"
  return(re)  
}


prepScores2 <- function(Z, formula, type="continuous", SNPInfo=NULL, snpNames="Name", aggregateBy="gene", kins=NULL, sparse=TRUE, data=parent.frame(), verbose=FALSE) {
  
  if(is.null(SNPInfo)){ 
    warning("No SNP Info file provided: loading the Illumina HumanExome BeadChip. See ?SNPInfo for more details")
    load(paste(find.package("skatMeta"), "data", "SNPInfo.rda",sep = "/"))
    aggregateBy = "SKATgene"
  } else {
    SNPInfo <- prepSNPInfo(SNPInfo, snpNames, aggregateBy)
  }
  
  check_format_skat(Z, SNPInfo, data, snpNames, kins)
  
  m <- create_model(formula, type, kins=kins, sparse=sparse, data=data) 
  
  maf <- calculate_maf(Z)
  Z <- impute_to_mean(Z)  
  scores <- colSums(m$res*Z)  
  
  re <- calculate_cov(Z, m, SNPInfo, snpNames, aggregateBy, kins)
  
  create_seqMeta(re, scores, maf, m, SNPInfo, snpNames, aggregateBy) 
}