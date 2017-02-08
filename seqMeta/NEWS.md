seqMeta 1.6.6
==================

-   c.seqMeta was not exported (#26)
-   version 1.42 of CompQuadForm changed the *name* of the results  (#28)


seqMeta 1.6.5
==================

-   Updated to comply with CRAN policy to explicitly import all functions in 
packages other than base (ex. stats, utils).   (#10)
-   check_format_skat relied on calling environment to define formula.  (#17)
-   burdenMeta and singlesnpMeta produced warnings if the standard error was 0.
Now returns \code{NA}. (#12, #13)
-   prepScores and prepScores2 scaled theta incorrectly when using a kinship 
matrix.  (#15, #16)
-   prep family of functions incorrectly stated that the returned maf of a 
seqMeta object was the "minor allele frequency" when in reality it is the 
"alternate allele frequency". (#14)
-   make sure coxme version is at least 2.2-4.  The definition of variance 
component estimates in lmekin changed. (#20)
-   prepCox function calculates the projection matrix incorrectly for 
collinear variants in the same gene. (#21)
-   internal function impute_to_mean could generate an error. (#23)
-   added examples to prepScores2 documentation. (#7)
-   added verbose functionality to prepScores2. (#11)