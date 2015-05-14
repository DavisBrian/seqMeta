skatMeta <- function(..., SNPInfo=NULL, wts = function(maf){dbeta(maf,1,25)}, method = "saddlepoint", snpNames = "Name", aggregateBy = "gene", mafRange = c(0,0.5), verbose=FALSE){
	cl <- match.call(expand.dots = FALSE)
	if(is.null(SNPInfo)){ 
		warning("No SNP Info file provided: loading the Illumina HumanExome BeadChip. See ?SNPInfo for more details")
		load(paste(find.package("seqMeta"), "data", "SNPInfo.rda",sep = "/"))
		aggregateBy = "SKATgene"
	}
	genelist <- na.omit(unique(SNPInfo[,aggregateBy]))
	cohortNames <- lapply(cl[[2]],as.character)
	ncohort <- length(cohortNames)
	
	ev <- parent.frame()
	classes <- unlist(lapply(cohortNames,function(name){class(get(name,envir=ev))})) 
	if(!all(classes == "seqMeta" | classes == "skatCohort") ){
	 	stop("an argument to ... is not a seqMeta object!")
	}
	
	res.strings <- data.frame("gene"=genelist,stringsAsFactors=F)
	res.numeric <- matrix(NA, nrow= nrow(res.strings),ncol =  length(c("p","Qmeta","cmaf","nmiss", "nsnps")))
	colnames(res.numeric) <- c("p","Qmeta","cmaf","nmiss", "nsnps")		
	
    if(verbose){
    	cat("\n Meta Analyzing... Progress:\n")
    	pb <- txtProgressBar(min = 0, max = length(genelist), style = 3)
    	pb.i <- 0
    }
    ri <- 0
    snp.names.list <- split(SNPInfo[,snpNames],SNPInfo[,aggregateBy])
	for(gene in genelist){		
		ri <- ri+1
		nsnps.sub <- length(snp.names.list[[gene]])
		
		n.miss <- n.total <- mscores <- maf <- numeric(nsnps.sub)
		big.cov <- matrix(0, nsnps.sub,nsnps.sub)
		
		vary.ave <- 0
		for(cohort.k in 1:ncohort){
			cohort.gene <- get(cohortNames[[cohort.k]],envir=ev)[[gene]]
			
			if(!is.null(cohort.gene)){
				#cohort.gene <- lapply(cohort.gene,function(x){replace(x,is.nan(x),0)})
				sub <- match(snp.names.list[[gene]],colnames(cohort.gene$cov))
				if(any(is.na(sub)) | any(sub != 1:length(sub), na.rm=TRUE) | length(cohort.gene$maf) > nsnps.sub){
							#if(any(is.na(sub))) warning("Some SNPs were not in SNPInfo file for gene ", gene," and cohort ",names(cohorts)[cohort.k])
							cohort.gene$cov <- as.matrix(cohort.gene$cov)[sub,sub,drop=FALSE]
							cohort.gene$cov[is.na(sub),] <- cohort.gene$cov[,is.na(sub)] <- 0
							
							cohort.gene$maf <- cohort.gene$maf[sub]
							cohort.gene$maf[is.na(sub)] <- -1
							
							cohort.gene$scores <- cohort.gene$scores[sub]
							cohort.gene$scores[is.na(sub)] <- 0
					}				
					
					n.total <- n.total + (cohort.gene$maf >= 0)*cohort.gene$n
					n.miss[cohort.gene$maf < 0] <- n.miss[cohort.gene$maf < 0] + cohort.gene$n
					cohort.gene$maf[cohort.gene$maf < 0] <- 0
					
					mscores <- mscores + cohort.gene$scores/cohort.gene$sey^2
					maf <- maf + 2*cohort.gene$maf*(cohort.gene$n)
					big.cov <- big.cov + cohort.gene$cov/cohort.gene$sey^2
					vary.ave <- vary.ave + max(cohort.gene$n,na.rm=T)*cohort.gene$sey^2
			}else{
				n.miss <- n.miss + get(cohortNames[[cohort.k]],envir=parent.frame())[[1]]$n
			} 
		}
		if(any(maf >0)){ 
			maf <- maf/(2*n.total)
			maf[is.nan(maf)] <- 0

			flip <- maf > 0.5
			mscores[flip] <- -mscores[flip]
			big.cov[flip,!flip] <- -big.cov[flip,!flip]
			big.cov[!flip,flip] <- -big.cov[!flip,flip]
			maf <- pmin(maf,1-maf)
		}
		if(is.function(wts)){
			tmpwts <- ifelse(maf > 0, wts(maf),0)
		} else if(is.character(wts)){
			tmpwts <- as.numeric(SNPInfo[SNPInfo[,aggregateBy]==gene,wts])
		} else {
			tmpwts <- rep(1,length(maf))
		}
		tmpwts <- as.numeric(tmpwts)
		
		if( !all(mafRange == c(0,0.5))){
		  keep <- (maf >= min(mafRange)) & (maf <= max(mafRange))
      	  tmpwts[!keep] <- 0
		}
    
		if(length(maf) > 0){
		  Qmeta <- sum((tmpwts*mscores)^2, na.rm=TRUE)
		  wcov <- tmpwts*t(t(big.cov)*tmpwts)
			lambda<-eigen(zapsmall(wcov),symmetric=TRUE)$values
    } else {
    	Qmeta <- 0
      lambda <- 0
    }
    	
    if(any(lambda > 0)){
    		p<-pchisqsum2(Qmeta,lambda,method=method)$p
		} else {
			p<-1
		}
		res.numeric[ri,"p"] = p
		res.numeric[ri,"Qmeta"] = Qmeta
		res.numeric[ri,"cmaf"] = sum(maf[tmpwts > 0],na.rm=TRUE)
		res.numeric[ri,"nsnps"] = sum(maf[tmpwts > 0] != 0, na.rm =T)
		res.numeric[ri,"nmiss"] = sum(n.miss, na.rm =T)
		if(verbose){
			pb.i <- pb.i+1
			setTxtProgressBar(pb, pb.i)
		}
	}	
	if(verbose) close(pb)
	return(cbind(res.strings,res.numeric))
}
