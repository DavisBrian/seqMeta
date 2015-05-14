\name{prepScores}
\alias{prepScores}
\alias{prepScoresX}
\alias{prepCox}
%\alias{skatCoxCohort}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Prepare scores for region based (meta) analysis  
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This function computes and organizes the neccesary output to efficiently meta-analyze SKAT and other tests. Note that the tests are *not* computed by these functions. The output must be passed to one of \code{\link{skatMeta}}, \code{\link{burdenMeta}}, or \code{\link{singlesnpMeta}}.

Unlike the SKAT package which operates on one gene at a time, these functions are intended to operate on many genes, e.g. a whole exome, to facilitate meta analysis of whole genomes or exomes.
}
\usage{
prepScores(Z, formula, family = gaussian(), SNPInfo=NULL, 
	snpNames = "Name", aggregateBy = "gene", kins=NULL, sparse = TRUE,
	data=parent.frame(), verbose = FALSE)

prepScoresX(Z, formula, male, family = gaussian(), SNPInfo=NULL, 
  snpNames = "Name", aggregateBy = "gene", kins=NULL, sparse = TRUE,
	data=parent.frame(), verbose = FALSE)

prepCox(Z, formula, SNPInfo=NULL, snpNames = "Name", aggregateBy = "gene", 
	data=parent.frame(), verbose =FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Z}{
	A genotype matrix (dosage matrix) - rows correspond to individuals and columns correspond to SNPs. Use 'NA' for missing values. The column names of this matrix should correspond to SNP names in the SNP information file.
	}
  \item{formula}{
	Base formula, of the kind used in glm() - typically of the form y~covariate1 + covariate2. For Cox models, the formula follows that of the coxph() function.
}
\item{male}{
  For analyzing the X chromosome, with prepScoresX, `male' is the gender (0/1 or F/T) indicating female/male. See details.
}
  \item{family}{
 	either gaussian(), for continuous data, or binomial() for 0/1 outcomes. Binary outcomes are not currently supported for family data.
	}
  \item{SNPInfo}{
	SNP Info file - must contain fields given in 'snpName' and 'aggregateBy'. 	
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
  \item{verbose}{
  	logical. whether or not to print the progress bar.
  }
  \item{kins}{
  	the kinship matrix for related individuals. Only supported for family=gaussian(). See lmekin and the kinship2 package for more details.
  }
  \item{sparse}{
  	whether or not to use a sparse Matrix approximation for dense kinship matrices (defaults to TRUE).
  }
}
\details{
	This function computes the neccesary information to meta analyze SKAT analyses: the individual SNP scores, their MAF, and a covariance matrix for each unit of aggregation. Note that the SKAT test is *not* calculated by this function. The output must be passed to one of \code{\link{skatMeta}},\code{\link{burdenMeta}}, or \code{\link{singlesnpMeta}}.
	
A crucial component of SKAT and other region-based tests is a common unit of aggregation accross studies. This is given in the SNP information file (argument \code{SNPInfo}), which pairs SNPs to a unit of aggregation (typically a gene). The additional arguments \code{Name} and \code{aggregateBy} specify the columns of the SNP information file which contain these pairings. Note that the column names of the genotype matrix \code{Z} must match the names given in the \code{Name} field.
	
Using \code{prepScores}, users are strongly recommended to use all SNPs, even if they are monomorphic in your study. This is for two reasons; firstly, monomorphic SNPs provide information about MAF across all studies; without providing the information we are unable to tell if a missing SNP data was monomorphic in a study, or simply failed to genotype adequately in that study.  Second, even if some SNPs will be filtered out of a particular meta-analysis (e.g., because they are intronic or common) construct ing seqMeta objects describing all SNPs will reduce the workload for subsequent follow-up analyses.

In order to analyze the X chromosome, use \code{prepScoresX}, with the additional argument \code{male}. Here, it is assumed that \code{Z} is coded 0/2 for males, 0/1/2 for females. \code{male} is needed to calculate of allele frequencies, used in imputation, and for weighting. Output from prepScoresX is handled just as output from prepScores. 

Note: to view results for a single study, one can pass a single seqMeta object to a function for meta-analysis.		
}
\value{
	an object of class 'seqMeta'. This is a list, not meant for human consumption, but to be fed to \code{skatMeta()} or another function.
	The names of the list correspond to gene names. Each element in the list contains
	\item{scores}{
		The scores 	(y-yhat)^t g
	}
	\item{cov}{
		The variance of the scores. When no covariates are used, this is the LD matrix.
	}
	\item{n}{
		The number of subjects
	}
	\item{maf}{
		The minor allele frequency
	}
	\item{sey}{
		The residual standard error.
	}


}
\references{
Wu, M.C., Lee, S., Cai, T., Li, Y., Boehnke, M., and Lin, X. (2011) Rare Variant Association Testing for Sequencing Data Using the Sequence Kernel Association Test (SKAT). American Journal of Human Genetics.
		
Chen H, Meigs JB, Dupuis J. Sequence Kernel Association Test for Quantitative Traits in Family Samples.  Genetic Epidemiology. (To appear)

Lin, DY and Zeng, D. On the relative efficiency of using summary statistics versus individual-level data in meta-analysis. Biometrika. 2010.
}
\author{
Arie Voorman, Jennifer Brody
}
\note{
	For \code{prepCox}, the signed likelihood ratio statistic is used instead of the score, as the score test is anti-conservative for proportional hazards regression. The code for this routine is based on the \code{coxph.fit} function from the \code{survival} package.

Please see the package vignette for more details.
}

\seealso{
\code{\link{skatMeta}}
\code{\link{burdenMeta}}
\code{\link{singlesnpMeta}}
\code{\link{skatOMeta}}
\code{\link{coxph}}
}
\examples{
	
###load example data for two studies:
### see ?seqMetaExample	
data(seqMetaExample)

####run on each cohort:
cohort1 <- prepScores(Z=Z1, y~sex+bmi, SNPInfo = SNPInfo, data =pheno1)
cohort2 <- prepScores(Z=Z2, y~sex+bmi, SNPInfo = SNPInfo, 
	kins=kins, data=pheno2)

#### combine results:
##skat
out <- skatMeta(cohort1, cohort2, SNPInfo = SNPInfo)
head(out)

##T1 test
out.t1 <- burdenMeta(cohort1,cohort2, SNPInfo = SNPInfo, mafRange = c(0,0.01))
head(out.t1)

##single snp tests:
out.ss <- singlesnpMeta(cohort1,cohort2, SNPInfo = SNPInfo)
head(out.ss)

\dontrun{
########################
####binary data

cohort1 <- prepScores(Z=Z1, ybin~1, family=binomial(), SNPInfo = SNPInfo, data =pheno1)
out <- skatMeta(cohort1, SNPInfo = SNPInfo)
head(out)


####################
####survival data
cohort1 <- prepCox(Z=Z1, Surv(time,status)~strata(sex)+bmi, SNPInfo = SNPInfo, 
	data =pheno1)
out <- skatMeta(cohort1, SNPInfo = SNPInfo)
head(out)
}
}
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
