seqMeta 1.6.0.9xxx
==================

-   Updated to comply with CRAN policy to explicitly import all functions in 
packages other than base (ex. stats, utils).   (#10)
-   check_format_skat relied on calling environment to define formula.  (#17)
-   burdenMeta and singlesnpMeta produced warnings if the standard error was 0.  Now returns \code{NA}. (#12, #13)
-   prepScores and prepScores2 scaled theta incorrectly when using a kinship matrix.  (#15, #16)
-   prep family of functions incorrectly stated that the returned maf of a seqMeta object was the "minor allele frequency" when in reality it is the "alternate allele frequency". (#14)