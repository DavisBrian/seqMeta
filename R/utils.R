# prepSNPInfo
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
