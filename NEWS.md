<!-- NEWS.md is generated from NEWS.Rmd. Please edit that file -->
seqMeta 1.5.0.9000
==================

-   Migrated to git / github
-   Minimum R version moved to 3.1.0
-   Added automated tests to test prepSNPInfo and prepScores2

New Function
------------

### prepScores2

prepScores2 is a drop in replacement for prepScores. The only difference is the family argument should be text. `gaussian()` becomes `"gaussian"`. Although the former will **temporarily** still work. Eventually I'd like to merge prepScoresX and prepCox into this function as 90+% of the code is the same.

-   removes duplicate snp-gene pairs in SNPInfo(fixes issue \#0001)
-   Reorganized code
-   Removed second imputation of genotypes in the covariance calculation
-   Removed orphaned code
-   Moved data check before imputation
-   Moved checks closer to where the data is calculated / used
-   Replaced `any(is.na(Z))` with `anyNA(Z)` which is faster and reduces memory usage
-   Replaced `range` with `min` and `max` which is faster and reduces memory usage
-   Replaced genotype imputation `apply` code with a faster vectorized version.
-   Replaced matrix operations of form `t(x)%*%y` with `crossprod(x, y)` is faster and reduces memory usage
-   Removed duplicate calculations
-   Added Error Messages for:
-   Missing weight for snp-gene pair
-   Duplicate weights for snp-gene pair
-   Non-numeric genotypes
-   family != gaussian or binomial
-   Added Warning Messages for:
-   Removing missing snps and genes from SNPInfo
-   Removing duplicated snp-gene pairs

Todo
====

-   change genotype documentation to "matrix or an object which can be coerced to one with as.matrix"
-   add back verbose option to `prepScores2`
-   prepScores2 support for male
-   prepScores2 support for Survival
-   check kins and Z are in the same order?
-   need to verify `impute_to_mean` code does NOT make a copy of the genotype matrix if it is a matrix
