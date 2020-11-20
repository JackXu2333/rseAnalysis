context("generate proper mutated RNA sequence")
library(dseAnalysis)

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

test_that("Mutation is working", {

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
