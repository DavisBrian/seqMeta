<!-- README.md is generated from README.Rmd. Please edit that file -->
seqMeta
=======

[![Build Status](https://travis-ci.org/DavisBrian/seqMeta.svg?branch=master)](https://travis-ci.org/DavisBrian) [![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/DavisBrian/seqMeta?branch=master)](https://ci.appveyor.com/project/DavisBrian/seqMeta) [![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/seqMeta)](https://CRAN.R-project.org/package=seqMeta)

Meta-Analysis of Region-Based Tests of Rare DNA Variants

Computes necessary information to meta analyze region-based tests for rare genetic variants (e.g. SKAT, T1) in individual studies, and performs meta analysis.

You can install:
----------------

-   the latest released version from CRAN with

    ``` r
    install.packages("seqMeta")
    ```

-   the latest development version from github with

    ``` r
    if (packageVersion("devtools") < 1.6) {
      install.packages("devtools")
    }
    devtools::install_github("DavisBrian/seqMeta")
    ```

If you encounter a clear bug, please file a minimal reproducible example on [github](https://github.com/DavisBrian/seqMeta/issues).

seqMet 1.6.6
------------

-   Export c.seqMeta
-   Fixed version 1.4.2 of CompQuadForm changing the *name* of the result
-   Fixed canonical form CRAN url

seqMeta 1.6.5
-------------

-   Explicitly passed formula to check\_format\_skat. Previously relied on calling environment to define formula.
-   Return `NA` if standard error is 0 in burdenMeta and singlesnpMeta.
-   prepScores and prepScores2 scaled theta incorrectly when using a kinship matrix.
-   Minimum coxme version moved to 2.2-4.
-   prepCox function calculates the projection matrix incorrectly for collinear variants in the same gene.
-   Fixed impute\_to\_mean from producing an error in odd circumstances.
-   Added verbose functionality to prepScores2.
-   Added examples to prepScores2 documentation.
-   Explicitly import all functions in packages other than base to comply with new CRAN policy.

seqMeta 1.6.0
-------------

-   Migrated to git / github
    -   Bug Reports and Feature Requests should be submitted [github](https://github.com/DavisBrian/seqMeta/issues).
-   Minimum R version moved to 3.1.0
-   Duplicated SNP in snpinfo gene no longer pulls from the genotype matrix twice.
-   Monomorphic snps with caf != 0 were handled incorrectly.
-   Binomial models when genotypes imputed outside of seqMeta did not match when models were imputed by seqMeta. Very slight differences in the covariance structure.
-   Replaced `any(is.na(Z))` with `anyNA(Z)`
-   Range test now checks that genotypes are \[0, 2\].
-   SNPInfo in seqMetaExamples had incorrect type of snpNames and aggregateBy.
-   Automatically convert (with warning) aggregateBy and snpName columns to type character if they are not already.
-   Added new function prepScores2

### prepScores2

prepScores2 is a drop in replacement for prepScores, prepScoresX and prepCox. The only difference is the family argument should be text. `gaussian()` becomes `"gaussian"`, `binomial()` becomes `"binomial"` and `"cox"` is used for survival models. prepScores2 is much faster in cases where genotype imputation occurs.

-   enforces assumption that when a gene is being anlyzed the same snp can only be in that gene once (same snp can still be in multiple genes or in the same gene in SNPInfo)
-   Reorganized code
-   Moved data checks before imputation of genotypes
-   Moved other checks closer to where the data is calculated / used
-   Improved speed and reduced memory usage by
    -   Replacing genotype imputation `apply` code with a faster vectorized version.
    -   Replacing `any(is.na(Z))` with `anyNA(Z)`
    -   Replacing `range` with `min` and `max`
    -   Replacing matrix operations of the form `t(x)%*%y` with `crossprod(x, y)`
    -   Removing duplicate calculations
    -   Removing orphaned code
    -   Removing second imputation of genotypes in the covariance calculation
    -   Converting genotype matrix of type data.frame to a matrix.
-   Added Error Messages for
    -   Missing weight for snp-gene pair
    -   Duplicate weights for snp-gene pair
    -   Non-numeric genotypes
    -   family != gaussian or binomial
-   Added Warning Messages for:
    -   Removing missing snps and genes from SNPInfo
    -   Removing duplicated snp-gene pairs
