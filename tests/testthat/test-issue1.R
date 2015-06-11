context("Issue 1")
# Duplicated SNP in a snpinfo gene pulls from the genotype matrix twice.
#
# This code affects all functions that take SNPInfo as an argument
#
# TBD
# - prepScores2 prepCox equivalent

data(seqMetaExample)

test_that("Original prep functions duplicated SNPS in snpinfo gene only get counted once)", {
  si <- SNPInfo[SNPInfo$gene %in% "gene1", ]
  si_dups <- rbind(si, si[1,])
  
  snps_gene1 <- as.character(intersect(colnames(Z1), si$Name))
  Zgene1 <- Z1[ , snps_gene1]
  
  #### prepScores ####
  cohort1 <- prepScores(Z=Zgene1, y~sex+bmi, SNPInfo = si, data =pheno1)
  cohort2 <- prepScores(Z=Zgene1, y~sex+bmi, SNPInfo = si_dups, data =pheno1)
  expect_equal(cohort1, cohort2)
  expect_equal(length(cohort1$gene1$scores), 15)
  expect_equal(length(cohort2$gene1$scores), 15)
  expect_equal(length(cohort1$gene1$maf), 15)
  expect_equal(length(cohort2$gene1$maf), 15)
  expect_equal(nrow(cohort1$gene1$cov), 15)
  expect_equal(nrow(cohort2$gene1$cov), 15)
  
  # singlesnpMeta 
  single1 <- singlesnpMeta(cohort1, SNPInfo=si, studyBetas=FALSE)
  single2 <- singlesnpMeta(cohort2, SNPInfo=si_dups, studyBetas=FALSE)
  expect_equal(single1, single2)
  expect_equal(nrow(single1), 15)
  expect_equal(nrow(single2), 15)
  
  # burdenMeta
  b1 <- burdenMeta(cohort1, SNPInfo=si)
  b2 <- burdenMeta(cohort2, SNPInfo=si_dups)
  expect_equal(b1, b2)
  
  # skatMeta
  s1 <- skatMeta(cohort1, SNPInfo=si)
  s2 <- skatMeta(cohort2, SNPInfo=si_dups)
  expect_equal(s1, s2)
  
  # skatOMeta
  sO1 <- skatMeta(cohort1, SNPInfo=si)
  sO2 <- skatMeta(cohort2, SNPInfo=si_dups)
  expect_equal(sO1, sO2)
  
  cohort1b <- prepScores(Z=Zgene1, ybin~1, family=binomial(), SNPInfo = si, data =pheno1)
  cohort2b <- prepScores(Z=Zgene1, ybin~1, family=binomial(), SNPInfo = si_dups, data =pheno1)
  expect_equal(cohort1b, cohort2b)
  expect_equal(length(cohort1b$gene1$scores), 15)
  expect_equal(length(cohort2b$gene1$scores), 15)
  expect_equal(length(cohort1b$gene1$maf), 15)
  expect_equal(length(cohort2b$gene1$maf), 15)
  expect_equal(nrow(cohort1b$gene1$cov), 15)
  expect_equal(nrow(cohort2b$gene1$cov), 15)
  
  #### prepScoresX ####
  cohort1 <- prepScoresX(Z=Zgene1, y~sex+bmi, male=pheno1$sex-1, SNPInfo = si, data =pheno1)
  cohort2 <- prepScoresX(Z=Zgene1, y~sex+bmi, male=pheno1$sex-1, SNPInfo = si_dups, data =pheno1)
  expect_equal(cohort1, cohort2)
  expect_equal(length(cohort1$gene1$scores), 15)
  expect_equal(length(cohort2$gene1$scores), 15)
  expect_equal(length(cohort1$gene1$maf), 15)
  expect_equal(length(cohort2$gene1$maf), 15)
  expect_equal(nrow(cohort1$gene1$cov), 15)
  expect_equal(nrow(cohort2$gene1$cov), 15)
  
  # singlesnpMeta 
  single1 <- singlesnpMeta(cohort1, SNPInfo=si, studyBetas=FALSE)
  single2 <- singlesnpMeta(cohort2, SNPInfo=si_dups, studyBetas=FALSE)
  expect_equal(single1, single2)
  expect_equal(nrow(single1), 15)
  expect_equal(nrow(single2), 15)
  
  # burdenMeta
  b1 <- burdenMeta(cohort1, SNPInfo=si)
  b2 <- burdenMeta(cohort2, SNPInfo=si_dups)
  expect_equal(b1, b2)
  
  # skatMeta
  s1 <- skatMeta(cohort1, SNPInfo=si)
  s2 <- skatMeta(cohort2, SNPInfo=si_dups)
  expect_equal(s1, s2)
  
  # skatOMeta
  sO1 <- skatMeta(cohort1, SNPInfo=si)
  sO2 <- skatMeta(cohort2, SNPInfo=si_dups)
  expect_equal(sO1, sO2)
  
  cohort1b <- prepScoresX(Z=Zgene1, ybin~1, male=pheno1$sex-1, family=binomial(), SNPInfo = si, data =pheno1)
  cohort2b <- prepScoresX(Z=Zgene1, ybin~1, male=pheno1$sex-1, family=binomial(), SNPInfo = si_dups, data =pheno1)
  expect_equal(cohort1b, cohort2b)
  expect_equal(length(cohort1b$gene1$scores), 15)
  expect_equal(length(cohort2b$gene1$scores), 15)
  expect_equal(length(cohort1b$gene1$maf), 15)
  expect_equal(length(cohort2b$gene1$maf), 15)
  expect_equal(nrow(cohort1b$gene1$cov), 15)
  expect_equal(nrow(cohort2b$gene1$cov), 15)
  
  
  #### prepCox ####
  cohort1 <- prepCox(Z=Zgene1, Surv(time,status)~strata(sex)+bmi, SNPInfo = si, data =pheno1)
  cohort2 <- prepCox(Z=Zgene1, Surv(time,status)~strata(sex)+bmi, SNPInfo = si_dups, data =pheno1)
  expect_equal(cohort1, cohort2)
  expect_equal(length(cohort1$gene1$scores), 15)
  expect_equal(length(cohort2$gene1$scores), 15)
  expect_equal(length(cohort1$gene1$maf), 15)
  expect_equal(length(cohort2$gene1$maf), 15)
  expect_equal(nrow(cohort1$gene1$cov), 15)
  expect_equal(nrow(cohort2$gene1$cov), 15)
  
  # singlesnpMeta 
  single1 <- singlesnpMeta(cohort1, SNPInfo=si, studyBetas=FALSE)
  single2 <- singlesnpMeta(cohort2, SNPInfo=si_dups, studyBetas=FALSE)
  expect_equal(single1, single2)
  expect_equal(nrow(single1), 15)
  expect_equal(nrow(single2), 15)
  
  # burdenMeta
  b1 <- burdenMeta(cohort1, SNPInfo=si)
  b2 <- burdenMeta(cohort2, SNPInfo=si_dups)
  expect_equal(b1, b2)
  
  # skatMeta
  s1 <- skatMeta(cohort1, SNPInfo=si)
  s2 <- skatMeta(cohort2, SNPInfo=si_dups)
  expect_equal(s1, s2)
  
  # skatOMeta
  sO1 <- skatMeta(cohort1, SNPInfo=si)
  sO2 <- skatMeta(cohort2, SNPInfo=si_dups)
  expect_equal(sO1, sO2)
  
  #### prepScores2 / prepScores equivalent ####
  cohort1 <- prepScores2(Z=Zgene1, y~sex+bmi, SNPInfo = si, data =pheno1)
  cohort2 <- prepScores2(Z=Zgene1, y~sex+bmi, SNPInfo = si_dups, data =pheno1)
  expect_equal(cohort1, cohort2)
  expect_equal(length(cohort1$gene1$scores), 15)
  expect_equal(length(cohort2$gene1$scores), 15)
  expect_equal(length(cohort1$gene1$maf), 15)
  expect_equal(length(cohort2$gene1$maf), 15)
  expect_equal(nrow(cohort1$gene1$cov), 15)
  expect_equal(nrow(cohort2$gene1$cov), 15)
  
  # singlesnpMeta 
  single1 <- singlesnpMeta(cohort1, SNPInfo=si, studyBetas=FALSE)
  single2 <- singlesnpMeta(cohort2, SNPInfo=si_dups, studyBetas=FALSE)
  expect_equal(single1, single2)
  expect_equal(nrow(single1), 15)
  expect_equal(nrow(single2), 15)
  
  # burdenMeta
  b1 <- burdenMeta(cohort1, SNPInfo=si)
  b2 <- burdenMeta(cohort2, SNPInfo=si_dups)
  expect_equal(b1, b2)
  
  # skatMeta
  s1 <- skatMeta(cohort1, SNPInfo=si)
  s2 <- skatMeta(cohort2, SNPInfo=si_dups)
  expect_equal(s1, s2)
  
  # skatOMeta
  sO1 <- skatMeta(cohort1, SNPInfo=si)
  sO2 <- skatMeta(cohort2, SNPInfo=si_dups)
  expect_equal(sO1, sO2)
  
  cohort1b <- prepScores2(Z=Zgene1, ybin~1, family=binomial(), SNPInfo = si, data =pheno1)
  cohort2b <- prepScores2(Z=Zgene1, ybin~1, family=binomial(), SNPInfo = si_dups, data =pheno1)
  expect_equal(cohort1b, cohort2b)
  expect_equal(length(cohort1b$gene1$scores), 15)
  expect_equal(length(cohort2b$gene1$scores), 15)
  expect_equal(length(cohort1b$gene1$maf), 15)
  expect_equal(length(cohort2b$gene1$maf), 15)
  expect_equal(nrow(cohort1b$gene1$cov), 15)
  expect_equal(nrow(cohort2b$gene1$cov), 15)
  })


test_that("prepScores2 duplicated SNPS in snpinfo gene only get counted once)", {
  si <- SNPInfo[SNPInfo$gene %in% "gene1", ]
  si_dups <- rbind(si, si[1,])
  
  snps_gene1 <- as.character(intersect(colnames(Z1), si$Name))
  Zgene1 <- Z1[ , snps_gene1]
  
  #### prepScores2 / prepScores equivalent ####
  cohort1 <- prepScores2(Z=Zgene1, y~sex+bmi, SNPInfo = si, data =pheno1)
  cohort2 <- prepScores2(Z=Zgene1, y~sex+bmi, SNPInfo = si_dups, data =pheno1)
  expect_equal(cohort1, cohort2)
  expect_equal(length(cohort1$gene1$scores), 15)
  expect_equal(length(cohort2$gene1$scores), 15)
  expect_equal(length(cohort1$gene1$maf), 15)
  expect_equal(length(cohort2$gene1$maf), 15)
  expect_equal(nrow(cohort1$gene1$cov), 15)
  expect_equal(nrow(cohort2$gene1$cov), 15)
  
  # singlesnpMeta 
  single1 <- singlesnpMeta(cohort1, SNPInfo=si, studyBetas=FALSE)
  single2 <- singlesnpMeta(cohort2, SNPInfo=si_dups, studyBetas=FALSE)
  expect_equal(single1, single2)
  expect_equal(nrow(single1), 15)
  expect_equal(nrow(single2), 15)
  
  # burdenMeta
  b1 <- burdenMeta(cohort1, SNPInfo=si)
  b2 <- burdenMeta(cohort2, SNPInfo=si_dups)
  expect_equal(b1, b2)
  
  # skatMeta
  s1 <- skatMeta(cohort1, SNPInfo=si)
  s2 <- skatMeta(cohort2, SNPInfo=si_dups)
  expect_equal(s1, s2)
  
  # skatOMeta
  sO1 <- skatMeta(cohort1, SNPInfo=si)
  sO2 <- skatMeta(cohort2, SNPInfo=si_dups)
  expect_equal(sO1, sO2)
  
  cohort1b <- prepScores2(Z=Zgene1, ybin~1, family=binomial(), SNPInfo = si, data =pheno1)
  cohort2b <- prepScores2(Z=Zgene1, ybin~1, family=binomial(), SNPInfo = si_dups, data =pheno1)
  expect_equal(cohort1b, cohort2b)
  expect_equal(length(cohort1b$gene1$scores), 15)
  expect_equal(length(cohort2b$gene1$scores), 15)
  expect_equal(length(cohort1b$gene1$maf), 15)
  expect_equal(length(cohort2b$gene1$maf), 15)
  expect_equal(nrow(cohort1b$gene1$cov), 15)
  expect_equal(nrow(cohort2b$gene1$cov), 15)
  
  #### prepScores2 / prepScoresX equivalent ####
  cohort1 <- prepScores2(Z=Zgene1, y~sex+bmi, male=pheno1$sex-1, SNPInfo = si, data =pheno1)
  cohort2 <- prepScores2(Z=Zgene1, y~sex+bmi, male=pheno1$sex-1, SNPInfo = si_dups, data =pheno1)
  expect_equal(cohort1, cohort2)
  expect_equal(length(cohort1$gene1$scores), 15)
  expect_equal(length(cohort2$gene1$scores), 15)
  expect_equal(length(cohort1$gene1$maf), 15)
  expect_equal(length(cohort2$gene1$maf), 15)
  expect_equal(nrow(cohort1$gene1$cov), 15)
  expect_equal(nrow(cohort2$gene1$cov), 15)
  
  # singlesnpMeta 
  single1 <- singlesnpMeta(cohort1, SNPInfo=si, studyBetas=FALSE)
  single2 <- singlesnpMeta(cohort2, SNPInfo=si_dups, studyBetas=FALSE)
  expect_equal(single1, single2)
  expect_equal(nrow(single1), 15)
  expect_equal(nrow(single2), 15)
  
  # burdenMeta
  b1 <- burdenMeta(cohort1, SNPInfo=si)
  b2 <- burdenMeta(cohort2, SNPInfo=si_dups)
  expect_equal(b1, b2)
  
  # skatMeta
  s1 <- skatMeta(cohort1, SNPInfo=si)
  s2 <- skatMeta(cohort2, SNPInfo=si_dups)
  expect_equal(s1, s2)
  
  # skatOMeta
  sO1 <- skatMeta(cohort1, SNPInfo=si)
  sO2 <- skatMeta(cohort2, SNPInfo=si_dups)
  expect_equal(sO1, sO2)
  
  cohort1b <- prepScores2(Z=Zgene1, ybin~1, male=pheno1$sex-1, family=binomial(), SNPInfo = si, data =pheno1)
  cohort2b <- prepScores2(Z=Zgene1, ybin~1, male=pheno1$sex-1, family=binomial(), SNPInfo = si_dups, data =pheno1)
  expect_equal(cohort1b, cohort2b)
  expect_equal(length(cohort1b$gene1$scores), 15)
  expect_equal(length(cohort2b$gene1$scores), 15)
  expect_equal(length(cohort1b$gene1$maf), 15)
  expect_equal(length(cohort2b$gene1$maf), 15)
  expect_equal(nrow(cohort1b$gene1$cov), 15)
  expect_equal(nrow(cohort2b$gene1$cov), 15)
    
  #### prepScores2 / prepScoresCox equivalent ####
  # TBD
})

