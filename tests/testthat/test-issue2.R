context("Issue 2")
# Monomorphic snps when maf != 0 handled incorrectly
#
# Fix affects prepScores and prepScoresX
#
# Sanity check on prepCox
# TBD
# - prepScores2 

data(seqMetaExample)


###########
# ISSUE 2 #
###########
## Monomophic SNP w/ maf !=0  (CAF == 1  || all hets) handled incorrectly by singlesnpMeta
test_that("Monomophic snps (all 3 cases) handled correctly - prepScores)", {
  si <- SNPInfo[SNPInfo$gene %in% "gene1", ]
  snps_gene1 <- as.character(intersect(colnames(Z1), si$Name))
  Zgene1 <- Z1[ , snps_gene1]
  
  monos <- c("1000011", "1000009", "1000001")
  Zgene1[ , monos[1] ] <- 0
  Zgene1[ , monos[2] ] <- 1
  Zgene1[ , monos[3] ] <- 2
  
  cohort1 <- prepScores(Z=Zgene1, y~sex+bmi, SNPInfo = si, data=pheno1)

  # check maf
  expect_equal(cohort1$gene1$maf[monos[1]], 0, check.attributes=FALSE)
  expect_equal(cohort1$gene1$maf[monos[2]], 0.5, check.attributes=FALSE)
  expect_equal(cohort1$gene1$maf[monos[3]], 1, check.attributes=FALSE)
  
  # check scores
  expect_true(cohort1$gene1$scores[monos[1]] == 0)
  expect_true(cohort1$gene1$scores[monos[2]] == 0)
  expect_true(cohort1$gene1$scores[monos[3]] == 0)

  # check cov
  expect_true(all(cohort1$gene1$cov[ , monos[1]] == 0))
  expect_true(all(cohort1$gene1$cov[ monos[1], ] == 0))
  expect_true(all(cohort1$gene1$cov[ , monos[2]] == 0))
  expect_true(all(cohort1$gene1$cov[ monos[2], ] == 0))
  expect_true(all(cohort1$gene1$cov[ , monos[3]] == 0))
  expect_true(all(cohort1$gene1$cov[ monos[3], ] == 0))
  
  
  out <- singlesnpMeta(cohort1, SNPInfo=si)
  expect_true(out[out$Name == monos[1], "maf"] == 0)
  expect_true(out[out$Name == monos[1], "caf"] == 0)
  expect_true(is.na(out[out$Name == monos[1], "p"]))
  expect_true(is.na(out[out$Name == monos[1], "beta"]))
  expect_true(is.na(out[out$Name == monos[1], "beta.cohort1"]))
  expect_true(is.infinite(out[out$Name == monos[1], "se"]))
  expect_true(is.infinite(out[out$Name == monos[1], "se.cohort1"])) 
  
  expect_true(out[out$Name == monos[2], "maf"] == 0.5)
  expect_true(out[out$Name == monos[2], "caf"] == 0.5)
  expect_true(is.na(out[out$Name == monos[2], "p"]))
  expect_true(is.na(out[out$Name == monos[2], "beta"]))
  expect_true(is.na(out[out$Name == monos[2], "beta.cohort1"]))
  expect_true(is.infinite(out[out$Name == monos[2], "se"]))
  expect_true(is.infinite(out[out$Name == monos[2], "se.cohort1"])) 
  
  expect_true(out[out$Name == monos[3], "maf"] == 0)
  expect_true(out[out$Name == monos[3], "caf"] == 1)
  expect_true(is.na(out[out$Name == monos[3], "p"]))
  expect_true(is.na(out[out$Name == monos[3], "beta"]))
  expect_true(is.na(out[out$Name == monos[3], "beta.cohort1"]))
  expect_true(is.infinite(out[out$Name == monos[3], "se"]))
  expect_true(is.infinite(out[out$Name == monos[3], "se.cohort1"])) 
  
})


test_that("Monomophic snps (all 3 cases) handled correctly - prepCox)", {
  si <- SNPInfo[SNPInfo$gene %in% "gene1", ]
  snps_gene1 <- as.character(intersect(colnames(Z1), si$Name))
  Zgene1 <- Z1[ , snps_gene1]
  
  monos <- c("1000011", "1000009", "1000001")
  Zgene1[ , monos[1] ] <- 0
  Zgene1[ , monos[2] ] <- 1
  Zgene1[ , monos[3] ] <- 2
  
  cohort1 <- prepCox(Z=Zgene1, Surv(time,status)~strata(sex)+bmi, SNPInfo = si, data =pheno1)
  
  # check maf
  expect_equal(cohort1$gene1$maf[monos[1]], 0, check.attributes=FALSE)
  expect_equal(cohort1$gene1$maf[monos[2]], 0.5, check.attributes=FALSE)
  expect_equal(cohort1$gene1$maf[monos[3]], 1, check.attributes=FALSE)
  
  # check scores
  expect_true(cohort1$gene1$scores[monos[1]] == 0)
  expect_true(cohort1$gene1$scores[monos[2]] == 0)
  expect_true(cohort1$gene1$scores[monos[3]] == 0)
  
  # check cov
  expect_true(all(cohort1$gene1$cov[ , monos[1]] == 0))
  expect_true(all(cohort1$gene1$cov[ monos[1], ] == 0))
  expect_true(all(cohort1$gene1$cov[ , monos[2]] == 0))
  expect_true(all(cohort1$gene1$cov[ monos[2], ] == 0))
  expect_true(all(cohort1$gene1$cov[ , monos[3]] == 0))
  expect_true(all(cohort1$gene1$cov[ monos[3], ] == 0))
  
  
  out <- singlesnpMeta(cohort1, SNPInfo=si)
  expect_true(out[out$Name == monos[1], "maf"] == 0)
  expect_true(out[out$Name == monos[1], "caf"] == 0)
  expect_true(is.na(out[out$Name == monos[1], "p"]))
  expect_true(is.na(out[out$Name == monos[1], "beta"]))
  expect_true(is.na(out[out$Name == monos[1], "beta.cohort1"]))
  expect_true(is.infinite(out[out$Name == monos[1], "se"]))
  expect_true(is.infinite(out[out$Name == monos[1], "se.cohort1"])) 
  
  expect_true(out[out$Name == monos[2], "maf"] == 0.5)
  expect_true(out[out$Name == monos[2], "caf"] == 0.5)
  expect_true(is.na(out[out$Name == monos[2], "p"]))
  expect_true(is.na(out[out$Name == monos[2], "beta"]))
  expect_true(is.na(out[out$Name == monos[2], "beta.cohort1"]))
  expect_true(is.infinite(out[out$Name == monos[2], "se"]))
  expect_true(is.infinite(out[out$Name == monos[2], "se.cohort1"])) 
  
  expect_true(out[out$Name == monos[3], "maf"] == 0)
  expect_true(out[out$Name == monos[3], "caf"] == 1)
  expect_true(is.na(out[out$Name == monos[3], "p"]))
  expect_true(is.na(out[out$Name == monos[3], "beta"]))
  expect_true(is.na(out[out$Name == monos[3], "beta.cohort1"]))
  expect_true(is.infinite(out[out$Name == monos[3], "se"]))
  expect_true(is.infinite(out[out$Name == monos[3], "se.cohort1"])) 
  
})

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
