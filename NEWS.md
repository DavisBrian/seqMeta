seqMeta 1.6.0
==================

-   Minimum R version moved to 3.1.0
-   Duplicated SNP in snpinfo gene no longer pulls from the genotype matrix twice.
-   Monomorphic snps with caf != 0 were handled incorrectly.
-   Binomial models when genotypes imputed outside of seqMeta did not match when models were imputed by seqMeta.  Very slight differences in the covariance structure.
-   Replaced `any(is.na(Z))` with `anyNA(Z)`
-   Range test now checks that genotypes are [0, 2].
-   SNPInfo in seqMetaExamples had incorrect type of snpNames and aggregateBy.
-   Automatically convert (with warning) aggregateBy and snpName columns to type character if they are not already.
-   Added new function prepScores2.  See README.