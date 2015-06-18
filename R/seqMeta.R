#' seqMeta: Meta-Analysis of Region-Based Tests of Rare DNA Variants
#'
#' Computes necessary information to meta analyze region-based tests for rare genetic variants (e.g. SKAT, T1) in individual studies, and performs meta analysis.
#'
#'
#' To learn more about seqMeta, start with the vignettes:
#' \code{browseVignettes(package = "seqMeta")}
#'
#' @name seqMeta
#' @docType package
#' @useDynLib seqMeta
#' @import Matrix CompQuadForm survival
#' @importFrom coxme lmekin
#' @export prepCondScores skatMeta burdenMeta singlesnpMeta skatOMeta
#' @exportMethod c.seqMeta c.seqMeta
NULL

# Illumina HumanExome BeadChip SNP Information file
# 
# Contains standard Names and associated genes for the Illumina HumanExome BeadChip
# 
# "SNPInfo"