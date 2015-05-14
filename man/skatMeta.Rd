\name{skatMeta}
\alias{skatMeta}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Combine SKAT analyses from one or more studies
}
\description{
Takes as input `seqMeta` objects (from e.g. \code{\link{prepScores}}), and meta analyzes them.
}
\usage{
skatMeta(..., SNPInfo=NULL, wts = function(maf){dbeta(maf,1,25)}, method = "saddlepoint", 
	snpNames = "Name", aggregateBy = "gene", mafRange = c(0,0.5), verbose=FALSE)
}
\arguments{
  \item{\dots}{
	seqMeta objects
}
  \item{SNPInfo}{
	the SNP Info file.  This should contain 'Name' and 'gene' fields, which match the 'Name' and 'gene' fields of the SNP Info file used in each study. Only SNPs and genes in this table will be meta analyzed, so this may be used to restrict the analysis.
		}
  \item{wts}{
	Either a function to calculate testing weights, or a character specifying a vector of weights in the SNPInfo file. For skatMeta the default are the `beta' weights.
	}
	  \item{method}{
	p-value calculation method. Default is 'saddlepoint', 'integration' is the Davies method used in the SKAT package. See pchisqsum() for more details.
}
  \item{snpNames}{
	The field of SNPInfo where the SNP identifiers are found.  Default is 'Name'
	}
  \item{aggregateBy}{
    The field of SNPInfo on which the skat results were aggregated.  Default is 'gene'. For single snps which are intended only for single variant analyses, it is reccomended that they have a unique identifier in this field.  
  }
  \item{mafRange}{
  	Range of MAF's to include in the analysis (endpoints included).  Default is all SNPs (0 <= MAF <= 0.5).
  }
  \item{verbose}{
  	logical. Whether or not to print progress bars.
  }
}
\details{
\code{skatMeta} implements an efficient SKAT meta analysis by meta-analyzing scores statistics and their variances.

Note: all studies must use coordinated SNP Info files - that is, the SNP names and gene definitions must be the same. 

Please see the package vignette for more details.
}
\value{
a data frame with columns:
\item{gene}{
	Name of the gene.
}
\item{p}{
	p-value of the SKAT test.
}
\item{Q}{
	The SKAT Q-statistic, defined as sum_j w_jS_j, where S_j is the squared score for SNP j, and w_j is a weight.
}
\item{cmaf}{
	The cumulative minor allele frequency.
}
\item{nmiss}{
	The number of `missing` SNPs. For a gene with a single SNP this is the number of individuals which do not contribute to the analysis, due to studies that did not report results for that SNP. For a gene with multiple SNPs, is totalled over the gene. 
}
\item{nsnps}{
	The number of SNPs in the gene.
}

}
\references{
Wu, M.C., Lee, S., Cai, T., Li, Y., Boehnke, M., and Lin, X. (2011) Rare Variant Association Testing for Sequencing Data Using the Sequence Kernel Association Test (SKAT). American Journal of Human Genetics.
}
\author{
Arie Voorman, Jennifer Brody
}
\seealso{
\code{\link{prepScores}}
\code{\link{burdenMeta}}
\code{\link{singlesnpMeta}}
\code{\link{skatOMeta}}
%\code{\link{pchisqsum}}
}
\examples{
###load example data for two studies:
### see ?seqMetaExample	
data(seqMetaExample)

####run on each study:
cohort1 <- prepScores(Z=Z1, y~sex+bmi, SNPInfo = SNPInfo, data =pheno1)
cohort2 <- prepScores(Z=Z2, y~sex+bmi, SNPInfo = SNPInfo, kins=kins, data=pheno2)

#### combine results:
##skat
out <- skatMeta(cohort1, cohort2, SNPInfo = SNPInfo)
head(out)

\dontrun{
##T1 test
out.t1 <- burdenMeta(cohort1,cohort2, SNPInfo = SNPInfo, mafRange = c(0,0.01))
head(out.t1)

##single snp tests:
out.ss <- singlesnpMeta(cohort1,cohort2, SNPInfo = SNPInfo)
head(out.ss)

########################
####binary data

cohort1 <- prepScores(Z=Z1, ybin~1, family=binomial(), SNPInfo = SNPInfo, data =pheno1)
out.bin <- skatMeta(cohort1, SNPInfo = SNPInfo)
head(out.bin)


####################
####survival data
cohort1 <- prepCox(Z=Z1, Surv(time,status)~strata(sex)+bmi, SNPInfo = SNPInfo, 
	data =pheno1)
out.surv <- skatMeta(cohort1, SNPInfo = SNPInfo)
head(out.surv)


##### Compare with SKAT on full data set
require(SKAT)
n <- nrow(pheno1)
bigZ <- matrix(NA,2*n,nrow(SNPInfo))
colnames(bigZ) <- SNPInfo$Name
for(gene in unique(SNPInfo$gene)){
	snp.names <- SNPInfo$Name[SNPInfo$gene == gene]
	bigZ[1:n,SNPInfo$gene == gene][,
		snp.names \%in\% colnames(Z1)] <- Z1[,
		na.omit(match(snp.names,colnames(Z1)))]
	bigZ[(n+1):(2*n),SNPInfo$gene == gene][,
		snp.names \%in\% colnames(Z2)] <- Z2[,
		na.omit(match(snp.names,colnames(Z2)))]
}
pheno <- rbind(pheno1[,c("y","bmi","sex")],pheno2[,c("y","bmi","sex")])
obj <- SKAT_Null_Model(y~sex+bmi+gl(2,nrow(pheno1)),data=pheno)
skat.pkg.p <- c(by(SNPInfo$Name, SNPInfo$gene,function(snp.names){
		inds <- match(snp.names,colnames(bigZ))
		if(sum(!is.na(inds)) ==0 ) return(1)
		SKAT(bigZ[,na.omit(inds)],obj, is_check = TRUE, missing = 1)$p.value
	}))

head(cbind(out$p,skat.pkg.p))

#Note: SKAT ignores family strucutre, resulting in p-values that are systematically too small: 
plot(y=out$p,x=skat.pkg.p, ylab = "SKAT meta p-values", xlab = "SKAT p-values")
abline(0,1)

#ignore family structure:
cohort1 <- prepScores(Z=Z1, y~sex+bmi, SNPInfo = SNPInfo, data =pheno1)
cohort2 <- prepScores(Z=Z2, y~sex+bmi, SNPInfo = SNPInfo, data =pheno2)

out.nofam <- skatMeta(cohort1,cohort2,SNPInfo=SNPInfo)
plot(y=out.nofam$p,x=skat.pkg.p, ylab = "SKAT meta p-values", xlab = "SKAT p-values")
abline(0,1)
}

}

