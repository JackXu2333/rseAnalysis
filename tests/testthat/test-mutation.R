context("generate proper mutated RNA sequence")
library(rseAnalysis)

test_that("DNA2RNA is working", {

  #Empty Input
  expect_error(RNA <- DNA2RNA(""))

  #Incorrect input type
  expect_error(RNA <- DNA2RNA(1))

  #Corrupted DNA
  expect_error(RNA <- DNA2RNA("ABC"))

  #check for result
  expect_equal(DNA2RNA("ATCG"), "AUCG")
  expect_equal(DNA2RNA("CCC"), "CCC")
  expect_equal(DNA2RNA("T"), "U")

})

test_that("Mutation and validation is working", {

  #No Input
  expect_error(mutate <- RNA.validate())

  #Load requirment file
  vcf <- vcf2df("../source/test.vcf")
  fasta <- fasta2df("../source/test.fasta")
  bed <- bed2df("../source/test.bed")

  #Input format error
  expect_error(mutate <- RNA.validate(fasta = data.frame(NAME = c("hsa", "hsa"), SEQQ = c("GGG", "GGG")), bed = bed, vcf = vcf))
  expect_error(mutate <- RNA.validate(fasta = fasta, bed = data.frame(NAME = c("has", "hsa")), vcf = vcf))
  expect_error(mutate <- RNA.validate(fasta = fasta, bed = bed, vcf = data.frame(NAME = c("has", "hsa"))))

  #Input size 0
  expect_error(mutate <- RNA.validate(fasta = fasta[0,], bed = bed, vcf = vcf))
  expect_error(mutate <- RNA.validate(fasta = fasta, bed = bed[0,], vcf = vcf))
  expect_error(mutate <- RNA.validate(fasta = fasta, bed = bed, vcf = vcf[0,]))

  (fasta)

  #Test correct running
  mutate <- RNA.validate(fasta = fasta, bed = bed, vcf = vcf)

  #Test direct output message
  expect_output(mutate <- RNA.validate(fasta = fasta, bed = bed, vcf = vcf), "The matching rate of the dataset is 0.111111111111111")

  expect_equal(mutate$MATCH, c(FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE))


})
