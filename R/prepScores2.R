#' @title Prepare scores for region based (meta) analysis
#'   
#' @description This function computes and organizes the neccesary output to 
#'   efficiently meta-analyze SKAT and other tests. Note that the tests are 
#'   *not* computed by these functions. The output must be passed to one of 
#'   \code{\link[seqMeta]{skatMeta}}, \code{\link[seqMeta]{burdenMeta}}, or 
#'   \code{\link[seqMeta]{singlesnpMeta}}.
#'   
#'   Unlike the SKAT package which operates on one gene at a time, these 
#'   functions are intended to operate on many genes, e.g. a whole exome, to 
#'   facilitate meta analysis of whole genomes or exomes.
#'   
#' @param Z A genotype matrix (dosage matrix) - rows correspond to individuals 
#'   and columns correspond to SNPs. Use 'NA' for missing values. The column 
#'   names of this matrix should correspond to SNP names in the SNP information 
#'   file.
#' @param formula Base formula, of the kind used in glm() - typically of the 
#'   form y~covariate1 + covariate2.
#' @param family either 'gaussian', for continuous data, or 'binomial' for 0/1 
#'   outcomes. Binary outcomes are not currently supported for family data.
#' @param SNPInfo SNP Info file - must contain fields given in 'snpName' and 
#'   'aggregateBy'.
#' @param snpNames The field of SNPInfo where the SNP identifiers are found. 
#'   Default is 'Name'
#' @param aggregateBy The field of SNPInfo on which the skat results were 
#'   aggregated. Default is 'gene'. For single snps which are intended only for 
#'   single variant analyses, it is reccomended that they have a unique 
#'   identifier in this field. for the SNP identifier is found.  See Details.
#' @param data  data frame in which to find variables in the formula
#' @param kins  the kinship matrix for related individuals. Only supported for 
#'   family=gaussian(). See lmekin and the kinship2 package for more details.
#' @param sparse  whether or not to use a sparse Matrix approximation for dense 
#'   kinship matrices (defaults to TRUE).
#' @param male For analyzing the X chromosome, with prepScoresX, 'male' is the 
#'   gender vector (0/1 or F/T) indicating female/male. See details.
#' @param verbose  logical. whether or not to print the progress bar.
#'   
#' @details This function computes the neccesary information to meta analyze 
#'   SKAT analyses: the individual SNP scores, their MAF, and a covariance 
#'   matrix for each unit of aggregation. Note that the SKAT test is *not* 
#'   calculated by this function. The output must be passed to one of 
#'   \code{\link[seqMeta]{skatMeta}}, \code{\link[seqMeta]{burdenMeta}}, or 
#'   \code{\link[seqMeta]{singlesnpMeta}}.
#'   
#'   A crucial component of SKAT and other region-based tests is a common unit 
#'   of aggregation accross studies. This is given in the SNP information file 
#'   (argument \code{SNPInfo}), which pairs SNPs to a unit of aggregation 
#'   (typically a gene). The additional arguments \code{snpNames} and 
#'   \code{aggregateBy} specify the columns of the SNP information file which 
#'   contain these pairings. Note that the column names of the genotype matrix 
#'   \code{Z} must match the names given in the \code{snpNames} field.
#'   
#'   Using \code{prepScores2}, users are strongly recommended to use all SNPs, 
#'   even if they are monomorphic in your study. This is for two reasons; 
#'   firstly, monomorphic SNPs provide information about MAF across all studies;
#'   without providing the information we are unable to tell if a missing SNP 
#'   data was monomorphic in a study, or simply failed to genotype adequately in
#'   that study. Second, even if some SNPs will be filtered out of a particular 
#'   meta-analysis (e.g., because they are intronic or common) constructing 
#'   seqMeta objects describing all SNPs will reduce the workload for subsequent
#'   follow-up analyses.
#'   
#'   Note: to view results for a single study, one can pass a single seqMeta 
#'   object to a function for meta-analysis.
#'   
#' @return an object of class 'seqMeta'. This is a list, not meant for human 
#'   consumption, but to be fed to \code{skatMeta()} or another function. The 
#'   names of the list correspond to gene names. Each element in the list 
#'   contains
#'   
#'   \describe{
#'   \item{scores}{The scores (y-yhat)^t g}
#'   \item{cov}{The variance of the scores. When no covariates are used, this is the LD matrix.}
#'   \item{n}{The number of subjects}
#'   \item{maf}{The minor allele frequency}
#'   \item{sey}{The residual standard error.}
#'   }
#'   
#' @note For \code{prepCox}, the signed likelihood ratio statistic is used 
#'   instead of the score, as the score test is anti-conservative for 
#'   proportional hazards regression. The code for this routine is based on the 
#'   \code{coxph.fit} function from the \code{survival} package.
#'   
#'   Please see the package vignette for more details.
#'   
#' @author Arie Voorman, Jennifer Brody
#' 
#' @references Wu, M.C., Lee, S., Cai, T., Li, Y., Boehnke, M., and Lin, X. (2011) Rare Variant Association Testing for Sequencing Data Using the Sequence Kernel Association Test (SKAT). American Journal of Human Genetics.
#' 
#' Chen H, Meigs JB, Dupuis J. Sequence Kernel Association Test for Quantitative Traits in Family Samples. Genetic Epidemiology. (To appear)
#' 
#' Lin, DY and Zeng, D. On the relative efficiency of using summary statistics versus individual-level data in meta-analysis. Biometrika. 2010.
#' 
#' @seealso \code{\link[seqMeta]{skatMeta}} \code{\link[seqMeta]{burdenMeta}} \code{\link[seqMeta]{singlesnpMeta}} \code{\link[seqMeta]{skatOMeta}} \code{\link[survival]{coxph}}
#' @export
prepScores2 <- function(Z, formula, family="gaussian", SNPInfo=NULL, snpNames="Name", aggregateBy="gene", kins=NULL, sparse=TRUE, data=parent.frame(), male=NULL, verbose=FALSE) {
  
  if(is.null(SNPInfo)){ 
    warning("No SNP Info file provided: loading the Illumina HumanExome BeadChip. See ?SNPInfo for more details")
    load(paste(find.package("skatMeta"), "data", "SNPInfo.rda",sep = "/"))
    aggregateBy = "SKATgene"
  } else {
    SNPInfo <- prepSNPInfo(SNPInfo, snpNames, aggregateBy)
  }
  
  if (!is.null(male)) {
    cl <- match.call()
    male <- eval(cl$male, data)
    if (anyNA(male)) {
      stop("Missing data not allowed in 'male'")
    }
    if(!all(male %in% c(0,1))) {
      stop("`male' must be coded as 0/1 or T/F")
    }
    if(length(male) != nrow(Z)) {
      stop("`male' not the same length as nrows genotype")
    }
    male <- as.logical(male)
  } 
  
  check_inputs(Z, SNPInfo, data, snpNames, kins)
  
  m <- create_model(formula, family, kins=kins, sparse=sparse, data=data) 
  
  maf <- calculate_maf(Z, male)
  Z <- impute_to_mean(Z, male)  
  scores <- colSums(m$res*Z)  
  
  re <- calculate_cov(Z, m, SNPInfo, snpNames, aggregateBy, kins)
  
  create_seqMeta(re, scores, maf, m, SNPInfo, snpNames, aggregateBy) 
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

# impute to mean
impute_to_mean <- function(Z, male=NULL) {
  
  if (!is.matrix(Z)) {
    Z <- as.matrix(Z)
  }
  
  if (anyNA(Z)) {
    if (is.null(male)) {
      MZ <- colMeans(Z, na.rm=TRUE)    
    } else {
      MZ <- (colSums(Z[male, ],na.rm=TRUE) + 2*colSums(Z[!male, ],na.rm=TRUE))/
        (colSums(!is.na(Z[male, ])) + 2*colSums(!is.na(Z[!male,])))
      
    }
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
create_model <- function(formula, family="gaussian", kins=NULL, sparse=TRUE, data=parent.frame()) {
  
  if (!is.character(family) || is.function(family) || is.family(family)) {
    fam <- family$family   
  } else {
    fam=family
  }
  if (is.null(fam)) {
    print(fam)
    stop("'family' not recognized")
  }
  
  if(!is.null(kins)){
    if (fam != "gaussian") {
      stop("Family data is currently only supported for continuous outcomes.")
    } 
    if(sparse){
      kins[kins < 2^{-5}] <- 0
      kins <- forceSymmetric(kins)
    }
    data$id <- if(is.null(colnames(kins))){
      1:ncol(kins)
    } else {
      colnames(kins)
    }    
    
    nullmodel <- lmekin(formula=update(formula, '~.+ (1|id)'), data=data, varlist = 2*kins,method="REML")  
    
    nullmodel$theta <- c(nullmodel$vcoef$id*nullmodel$sigma^2,nullmodel$sigma^2)   
    SIGMA <- nullmodel$theta[1] * 2 * kins + nullmodel$theta[2] * Diagonal(nrow(kins))   
    s2 <- sum(nullmodel$theta)
    
    #rotate data:
    nullmodel$family$var <- function(x){1}
    sef <- sqrt(nullmodel$family$var(nullmodel$fitted))
    X1 <- sef*model.matrix(lm(formula,data=data)) 
    res <- as.vector(nullmodel$res)* s2 / nullmodel$theta[2]  
    Om_i <- solve(SIGMA/s2)
    # optimize calculations
     tX1_Om_i <- crossprod(X1, Om_i)
     AX1 <- with(svd(tX1_Om_i%*%X1),  v[,d > 0,drop=FALSE]%*%( (1/d[d>0])*t(v[, d > 0,drop=FALSE])))%*%tX1_Om_i
    
    check_dropped_subjects(res, formula)
    list(res=res, family=fam, n=nrow(X1), sey=sqrt(s2), sef=sef, X1=X1, AX1=AX1, Om_i=Om_i)
  } else {
    nullmodel <- glm(formula=formula, family=fam, data=data)
    res <- residuals(nullmodel, type = "response")  
    check_dropped_subjects(res, formula) 
    sef <- sqrt(nullmodel$family$var(nullmodel$fitted))
    X1 <- sef*model.matrix(nullmodel) 
    AX1 <- with(svd(X1),  v[,d > 0,drop=FALSE]%*%( (1/d[d>0])*t(u[, d > 0,drop=FALSE])))     
    sey <- if (fam == "gaussian") {
      sqrt(var(res)*(nrow(X1) - 1)/(nrow(X1) - ncol(X1)) )
    } else if (fam == "binomial") {
      1
    } else {
      stop("Only family type 'binomial' and 'gaussian' are currently supported.")
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

is.family <- function(x) {"family" %in% class(x)}

# check_format_skat
check_inputs <- function(Z, SNPInfo, data, snpNames, kins){
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
    warning(paste(ncol(Z) - nsnps, "snps are not in SNPInfo file!"))
  } 
  
  invisible(NULL)
}

# calculate_maf <- function(Z) {
#   
#   MeanZ <- colMeans(Z, na.rm=TRUE)
#   maf <- MeanZ/2
# 
#   #differentiate all missing from monomorphic
#   maf[which(is.nan(MeanZ))] <- -1
#   
#   maf
# }
calculate_maf <- function(Z, male=NULL) {
  
  if (is.null(male)) {
    maf <- colMeans(Z, na.rm=TRUE)/2.0
  } else {
    male <- as.logical(male)
    maf <- (colSums(Z[male, ],na.rm=TRUE)/2 + colSums(Z[!male, ],na.rm=TRUE))/
      (colSums(!is.na(Z[male, ])) + 2*colSums(!is.na(Z[!male,])))    
  }
  
  #differentiate all missing from monomorphic
  maf[is.nan(maf)] <- -1
  
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
#      Z0 <- impute_to_mean(Z0)     # not sure how Z0 can have missing values
      if(!is.null(kins)){
        tZ0_Omi <- crossprod(Z0, m$Om_i)
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