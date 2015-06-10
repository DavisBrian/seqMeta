context("prepScores2")

data(seqMetaExample)


test_that("prepScores2 equals prepScores (gaussian w/o family)", {
  ps <- prepScores(Z=Z1, y~sex+bmi, SNPInfo = SNPInfo, data =pheno1)
  ps2a <- prepScores2(Z=Z1, y~sex+bmi, family="gaussian", SNPInfo = SNPInfo, data =pheno1)
  ps2b <- prepScores2(Z=Z1, y~sex+bmi, family=gaussian(), SNPInfo = SNPInfo, data =pheno1)
  
  expect_equal(ps2a, ps)
  expect_equal(ps2b, ps)
})


test_that("prepScores2 equals prepScores (gaussian w/ family)", {
  ps <- prepScores(Z=Z2, y~sex+bmi, SNPInfo = SNPInfo, kins=kins, data=pheno2)
  ps2a <- prepScores2(Z=Z2, y~sex+bmi, family="gaussian", SNPInfo = SNPInfo, kins=kins, data=pheno2)
  ps2b <- prepScores2(Z=Z2, y~sex+bmi, family=gaussian(), SNPInfo = SNPInfo, kins=kins, data=pheno2)
  
  expect_equal(ps2a, ps)
  expect_equal(ps2b, ps)
})

test_that("prepScores2 equals prepScores (binomial)", {
  ps <- prepScores(Z=Z1, ybin~1, family=binomial(), SNPInfo=SNPInfo, data=pheno1)
  ps2a <- prepScores2(Z=Z1, ybin~1, family="binomial", SNPInfo = SNPInfo, data=pheno1)
  ps2b <- prepScores2(Z=Z1, ybin~1, family=binomial(), SNPInfo = SNPInfo, data=pheno1)
  
  expect_equal(ps2a, ps)
  expect_equal(ps2b, ps)
})

###############
## SAME TESTS but for X
###############
test_that("prepScores2 equals prepScoresX (gaussian w/o family)", {
  ps <- prepScoresX(Z=Z1, y~sex+bmi, male=pheno1$sex-1, SNPInfo = SNPInfo, data =pheno1)
  ps2a <- prepScores2(Z=Z1, y~sex+bmi, male=pheno1$sex-1, family="gaussian", SNPInfo = SNPInfo, data =pheno1)
  ps2b <- prepScores2(Z=Z1, y~sex+bmi, male=pheno1$sex-1, family=gaussian(), SNPInfo = SNPInfo, data =pheno1)
  
  expect_equal(ps2a, ps)
  expect_equal(ps2b, ps)
})

test_that("prepScores2 equals prepScoresX (gaussian w/ family)", {
  ps <- prepScoresX(Z=Z2, y~sex+bmi, male=pheno2$sex-1, SNPInfo = SNPInfo, kins=kins, data=pheno2)
  ps2a <- prepScores2(Z=Z2, y~sex+bmi, male=pheno2$sex-1, family="gaussian", SNPInfo = SNPInfo, kins=kins, data=pheno2)
  ps2b <- prepScores2(Z=Z2, y~sex+bmi, male=pheno2$sex-1, family=gaussian(), SNPInfo = SNPInfo, kins=kins, data=pheno2)
  
  expect_equal(ps2a, ps)
  expect_equal(ps2b, ps)
})

test_that("prepScores2 equals prepScoresX (binomial)", {
  ps <- prepScoresX(Z=Z1, ybin~1, male=pheno1$sex-1, family=binomial(), SNPInfo=SNPInfo, data=pheno1)
  ps2a <- prepScores2(Z=Z1, ybin~1, male=pheno1$sex-1, family="binomial", SNPInfo = SNPInfo, data=pheno1)
  ps2b <- prepScores2(Z=Z1, ybin~1, male=pheno1$sex-1, family=binomial(), SNPInfo = SNPInfo, data=pheno1)
  
  expect_equal(ps2a, ps)
  expect_equal(ps2b, ps)
})


###############
## SAME TESTS but for Survival
###############



###############
## Missing values
###############
# test missing all missing genotypes in a snp gives the same as if the snp was removed from the genotype matrix


