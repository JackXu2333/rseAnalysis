context("gene expression and RNA distance analysis")
library(rseAnalysis)

test_that("Analysis input validation", {

  #None input
  expect_error(RNA <- Analysis.DISEXP())

  #Input type mismatch
  expect_error(RNA <- Analysis.DISEXP(dis.name = 1, dis.distance = as.integer(c(1)),
                                      exp.tumor = c(1), exp.sample = c(2), method = "linear"))
  expect_error(RNA <- Analysis.DISEXP(dis.name = c("1"), dis.distance = c(1),
                                      exp.tumor = c(1), exp.sample = c(2), method = "linear"))
  expect_error(RNA <- Analysis.DISEXP(dis.name = c("1"), dis.distance = as.integer(c(1)),
                                      exp.tumor = c("1"), exp.sample = c(2), method = "linear"))

  #Input out of bound
  expect_error(RNA <- Analysis.DISEXP(dis.name = c("1"), dis.distance = as.integer(c(1)),
                                      exp.tumor = c(1), exp.sample = c(2), method = "linlin"))

  #Input size does not aligned
  expect_error(RNA <- Analysis.DISEXP(dis.name = c("1", "1"), dis.distance = as.integer(c(1)),
                                      exp.tumor = c(1), exp.sample = c(2), method = "linear"))
  expect_error(RNA <- Analysis.DISEXP(dis.name = c("1"), dis.distance = as.integer(c(1)),
                                      exp.tumor = c(1, 2), exp.sample = c(2), method = "linear"))

})

test_that("Analysis working", {

  #Load sample files
  filePath <- system.file("extdata", "test.csv", package = "rseAnalysis")
  expression <- read.csv(filePath, header = TRUE)
  expression <- subset(expression, Read.Type == "reads_per_million_miRNA_mapped")[1:10, ]

  dis.name <- expression$mRNA
  dis.distance <- as.integer(c(12, 1, 34, 19, 103, 18, 45, 83, 49, 23))
  exp.tumor <- expression$Sample
  exp.sample <- expression$Normal

  #Test linear option
  RNA.lm <- Analysis.DISEXP(dis.name = dis.name, dis.distance = dis.distance,
                                    exp.tumor = exp.tumor, exp.sample = exp.sample, method = "linear")

  expect_equal(round(RNA.lm$stats$Correlation, digits = 3), -0.465)
  expect_equal(round(RNA.lm$stats$PValue, digits = 3), 0.176)


  #Test log option
  RNA.log <- Analysis.DISEXP(dis.name = dis.name, dis.distance = dis.distance,
                         exp.tumor = exp.tumor, exp.sample = exp.sample, method = "log")

  expect_equal(round(RNA.log$stats$Correlation, digits = 3), -0.013)
  expect_equal(round(RNA.log$stats$PValue, digits = 3), 0.203)


})
