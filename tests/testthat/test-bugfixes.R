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


library(seqMeta)
data(seqMetaExample)

### Issue 2
# make a small snpinfo to show the problem
# si <- SNPInfo[SNPInfo$gene %in% "gene1", ]
# 
# snps_gene1 <- as.character(intersect(colnames(Z1), si$Name))
# Zgene1 <- Z1[ , snps_gene1]
# 
# # make a 3 monomophic test case
# Zgene1[ , "1000011"] <- 0
# Zgene1[ , "1000009"] <- 1
# Zgene1[ , "1000001"] <- 2
# 
# 
# cohort1 <- prepScores2(Z=Zgene1, y~sex+bmi, SNPInfo = si, data=pheno1)
# # maf of Zgene1[ , "1000011"] == 0
# # maf of Zgene1[ , "1000009"] == 0.5
# # maf of Zgene1[ , "1000001"] == 1
# 
# # scores of Zgene1[ , "1000011"] == 0
# # scores of Zgene1[ , "1000009"] == 0
# # scores of Zgene1[ , "1000001"] == 0
# 
# # cov of Zgene1["1000011", ] == 0
# # cov of Zgene1[, "1000011"] == 0
# # cov of Zgene1["1000009", ] == 0
# # cov of Zgene1[, "1000009"] == 0
# # cov of Zgene1["1000001", ] == 0
# # cov of Zgene1[, "1000001"] == 0
# 
# singlesnpMeta(cohort1, SNPInfo=si)
# # maf same as above
# # caf min(maf, 1-maf)
# # p = NA
# # beta = NA
# # se = NA
# # beta.cohor1 = NA
# # se.cohort1 = Inf
# 
# 
# cohort1 <- prepScores(Z=Zgene1, ybin~1, family=binomial(), SNPInfo=si, data=pheno1)
# cohort1
# singlesnpMeta(cohort1, SNPInfo=si)
# # same as guassian case
# 
# 
# 
# 
# cohort1 <- prepScoresX(Z=Zgene1, male=pheno1$sex-1, y~sex+bmi, SNPInfo = si, data =pheno1)
# cohort1
# singlesnpMeta(cohort1, SNPInfo=si)
# 
# 
# cohort1 <- prepScoresX(Z=Zgene1, ybin~1, male=pheno1$sex-1, family=binomial(), SNPInfo=si, data=pheno1)
# cohort1
# 
# 
# 
# 
# cohort1 <- prepCox(Z=Zgene1, Surv(time,status)~strata(sex)+bmi, SNPInfo = si, data =pheno1)
# cohort1
# singlesnpMeta(cohort1, SNPInfo=si)



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


# si <- SNPInfo[SNPInfo$gene %in% "gene1", ]
# si_dups <- rbind(si, si[1,])
# 
# snps_gene1 <- as.character(intersect(colnames(Z1), si$Name))
# 
# 
# 
# prepScores(Z=Z1[ , snps_gene1], y~sex+bmi, SNPInfo = si, data =pheno1)
# prepScores(Z=Z1[ , snps_gene1], y~sex+bmi, SNPInfo = si_dups, data =pheno1)
# 
# 
# prepScores(Z=Z1[ , snps_gene1], ybin~1, family=binomial(), SNPInfo=si, data=pheno1)
# prepScores(Z=Z1[ , snps_gene1], ybin~1, family=binomial(), SNPInfo=si_dups, data=pheno1)
# 
# 
# prepScoresX(Z=Z1[ , snps_gene1], ybin~1, male=pheno1$sex-1, family=binomial(), SNPInfo=si, data=pheno1)
# prepScoresX(Z=Z1[ , snps_gene1], ybin~1, male=pheno1$sex-1, family=binomial(), SNPInfo=si_dups, data=pheno1)
# 
# 
# prepCox(Z=Z1[ , snps_gene1], Surv(time,status)~strata(sex)+bmi, SNPInfo = si, data =pheno1)
# prepCox(Z=Z1[ , snps_gene1], Surv(time,status)~strata(sex)+bmi, SNPInfo = si_dups, data =pheno1)
# 
