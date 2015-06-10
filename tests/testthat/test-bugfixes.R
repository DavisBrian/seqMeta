context("bug fixes")

data(seqMetaExample)

###########
# ISSUE N #
###########
## Template



###########
# ISSUE 3 #
###########
## Binomial imputation bug



###########
# ISSUE 2 #
###########
## Monomophic SNP w/ maf !=0  (CAF == 1  || all hets) handled incorrectly by singlesnpMeta
test_that("Monomophic snps (all 3 cases) handled correctly)", {
  ps <- prepScores(Z=Z1, y~sex+bmi, SNPInfo = SNPInfo, data =pheno1)
  ps2a <- prepScores2(Z=Z1, y~sex+bmi, family="gaussian", SNPInfo = SNPInfo, data =pheno1)
  ps2b <- prepScores2(Z=Z1, y~sex+bmi, family=gaussian(), SNPInfo = SNPInfo, data =pheno1)
  
  expect_equal(ps2a, ps)
  expect_equal(ps2b, ps)
})




###########
# ISSUE 1 #
###########
## Duplicated SNP in snpinfo gene pulls from the genotype matrix twice

test_that("Duplicated SNPS in snpinfo gene only get counted once)", {
  ps <- prepScores(Z=Z1, y~sex+bmi, SNPInfo = SNPInfo, data =pheno1)
  ps2a <- prepScores2(Z=Z1, y~sex+bmi, family="gaussian", SNPInfo = SNPInfo, data =pheno1)
  ps2b <- prepScores2(Z=Z1, y~sex+bmi, family=gaussian(), SNPInfo = SNPInfo, data =pheno1)
  
  expect_equal(ps2a, ps)
  expect_equal(ps2b, ps)
})




