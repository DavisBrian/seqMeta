\name{prepCondScores}
\alias{prepCondScores}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Run SKAT on data from a single cohort, conditional on specified SNP effects
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This function works exactly as \code{\link{prepScores}}, but with the additional argument `adjustments' specifying genes for which conditional analyses are desired, and which SNPs to condition on.
}
\usage{
prepCondScores(Z, formula, family = gaussian(), SNPInfo=NULL, 
	adjustments= NULL, snpNames = "Name", aggregateBy = "gene", 
	kins=NULL, sparse = TRUE, data=parent.frame())
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Z}{
	A genotype matrix (dosage matrix) - rows correspond to individuals and columns correspond to SNPs. Use 'NA' for missing values. The column names of this matrix should correspond to SNP names in the SNP information file.
	}
  \item{formula}{
	Base formula, of the kind used in glm() - typically of the form y~covariate1 + covariate2. For Cox models, the formula follows that of the coxph() function.
}
  \item{family}{
	Either gaussian(), for continuous data, or binomial() for 0/1 outcomes. Binary outcomes are not currently supported for family data.
	}
  \item{SNPInfo}{
	SNP Info file - must contain fields given in 'snpName' and 'aggregateBy'. 	
	}
  \item{adjustments}{
	A data frame of the same format at SNPInfo, pairing genes to analyze with snp 
  }
  \item{snpNames}{
	The field of SNPInfo where the SNP identifiers are found.  Default is 'Name'
	}
  \item{aggregateBy}{
    The field of SNPInfo on which the skat results were aggregated.  Default is 'gene'. For single snps which are intended only for single variant analyses, it is reccomended that they have a unique identifier in this field.  
  }
  \item{data}{
	data frame in which to find variables in the formula
  }	
  \item{kins}{
  	The kinship matrix for related individuals. See lmekin and the kinship2 package for more details
  }
  \item{sparse}{
  	whether or not to use a sparse Matrix approximation for dense kinship matrices (defaults to TRUE)
  }
}
\details{
	This function has the same syntax as \code{\link{prepCondScores}}, but requires an extra argument `adjustments`. This is a data frame of the same format as the SNPInfo, i.e. with a `snpNames` and `aggregateBy` columns. The function works by looping through the genes in the adjustment file, adding the corresponding SNPs to the null model.  For instance, if one wants to adjuste `gene1` for SNPs a and b (which need not be in gene 1), and `gene2' for SNPs c, the adjustments would be something like
	\code{adjustments = data.frame(Name = c("a","b","c"), gene = c("gene1","gene1","gene2"))}

See the examples for an illustration.
}
\value{
	an object of class 'seqMeta'. Note that unlike output from the function \code{\link{prepScores}}, the null models in each element of the list may be different. When meta analyzing these, it may be good to subset the SNPInfo file to the genes of interest.
}
\author{
Arie Voorman, Jennifer Brody
}

\seealso{
\code{\link{prepScores}}
\code{\link{skatMeta}}
\code{\link{burdenMeta}}
\code{\link{singlesnpMeta}}
}
\examples{
	
###load example data for two studies:
### see ?seqMetaExample	
data(seqMetaExample)

#specify adjustment variables
adjustments <- SNPInfo[c(1:3, 20,100), ]
adjustments


####run on each study:
cohort1.adj <- prepCondScores(Z=Z1, y~sex+bmi, SNPInfo = SNPInfo, 
	adjustments=adjustments, data =pheno1)
cohort2.adj <- prepCondScores(Z=Z2, y~sex+bmi, SNPInfo = SNPInfo, 
	adjustments=adjustments, kins=kins, data=pheno2)

SNPInfo.sub <- subset(SNPInfo, (SNPInfo$gene \%in\% adjustments$gene) & 
	!(SNPInfo$Name \%in\% adjustments$Name) )

#skat
out.skat <- skatMeta(cohort1.adj,cohort2.adj, SNPInfo = SNPInfo.sub)
head(out.skat)

##T1 test
out.t1 <- burdenMeta(cohort1.adj,cohort2.adj, SNPInfo = SNPInfo.sub, mafRange = c(0,0.01))
head(out.t1)

##single snp tests:
out.ss <- singlesnpMeta(cohort1.adj,cohort2.adj, SNPInfo = SNPInfo.sub)
head(out.ss)

}
