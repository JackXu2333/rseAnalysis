context("import function for file import")
library(dseAnalysis)

test_that("vcf2df file loading", {

  #non present folder
  expect_error(vcf <- vcf2df("./foo/test.vcf"))

  #corrupted file path
  expect_error(vcf <- vcf2df("\test.vcf"))

  #corrupted file type
  expect_error(vcf <- vcf2df("./tests/source/test.vff"))

  vcf <- vcf2df("../source/test.vcf")

  #check for result type
  expect_match(typeof(vcf$CHROM), "character")
  expect_match(typeof(vcf$POS), "integer")
  expect_match(typeof(vcf$ID), "character")
  expect_match(typeof(vcf$REF), "character")
  expect_match(typeof(vcf$ALT), "character")
  expect_match(typeof(vcf$CONSEQUENCE), "character")
  expect_match(typeof(vcf$OCCURRENCE), "character")
  expect_match(typeof(vcf$affected_donors), "integer")
  expect_match(typeof(vcf$project_count), "integer")

})

test_that("fasta2df file loading", {

  #non present folder
  expect_error(fasta <- fasta2df("./foo/test.fasta"))

  #corrupted file path
  expect_error(fasta <- fasta2df("\test.fasta"))

  #corrupted file type
  expect_error(fasta <- fasta2df("./tests/source/test.fastaaa"))

  fasta <- fasta2df("../source/test.fasta")

  #check for result type
  expect_match(typeof(fasta$NAME), "character")
  expect_match(typeof(fasta$SEQ), "character")

})


test_that("bed2df file loading", {

  #non present folder
  expect_error(bed <- bed2df("./foo/test.bed"))

  #corrupted file path
  expect_error(bed <- bed2df("\test.bed"))

  #corrupted file type
  expect_error(bed <- bed2df("./tests/source/test.beeed"))

  bed <- bed2df("../source/test.bed")

  #check for result type
  expect_match(typeof(bed$CHROM), "character")
  expect_match(typeof(bed$STAPOS), "integer")
  expect_match(typeof(bed$ENDPOS), "integer")
  expect_match(typeof(bed$DIR), "character")
  expect_match(typeof(bed$TYPE), "character")
  expect_match(typeof(bed$ID), "character")
  expect_match(typeof(bed$ALIAS), "character")
  expect_match(typeof(bed$NAME), "character")

})
