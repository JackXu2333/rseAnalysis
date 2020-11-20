context("secondary structure and distance prediction")
library(rseAnalysis)

test_that("File input for secondary structure is correct", {

  #executable path error
  expect_error(mutate <- predict.Structure(executable.path = "./foo", fasta.file = "../source/test.fasta"))

  #fasta file path error
  expect_error(mutate <- predict.Structure(executable.path = "../source", fasta.file = "./foo/datsa"))

  #fasta file unavailable with missing name or seq
  expect_error(mutate <- predict.Structure(executable.path = "../source", rna.name = c("heyhey", "blueblue")))

  #name or seq length mismatched
  expect_error(mutate <- predict.Structure(executable.path = "../source"
                                           , rna.name = c("hsa-1", "hsa-2"), rna.seq = c("GGG", "GGG", "AAA")))

  #name or seq type mismatched
  expect_error(mutate <- predict.Structure(executable.path = "../source"
                                           , rna.name = c(1, 2), rna.seq = c("GGG", "AAA")))

})

if (TRUE) skip("Unable to test under current envirment")

test_that("Secondary structure is working", {

  #Load datapath
  #setwd("./test/source")

  system(paste0("export DATAPATH=", gsub(" ", "\\ ", getwd(), fixed = TRUE), "/tests/source/data_tables"))

  #Load via file
  mutate.file <- predict.Structure(executable.path = "../source/exe", fasta.file = "../source/test.fasta")

  #Load via file
  mutate.input <- predict.Structure(executable.path = "../source/exe",
                                    rna.name = c("hsa-let-7a-1", "hsa-let-7a-2", "hsa-let-7a-3", "hsa-let-7b"),
                                    rna.seq = c("UGGGAUGAGGUAGUAGGUUGUAUAGUUUUAGGGUCACACCCACCACUGGGAGAUAACUAUACAAUCUACUGUCUUUCCUA",
                                            "AGGUUGAGGUAGUAGGUUGUAUAGUUUAGAAUUACAUCAAGGGAGAUAACUGUACAGCCUCCUAGCUUUCCU",
                                            "GGGUGAGGUAGUAGGUUGUAUAGUUUGGGGCUCUGCCCUGCUAUGGGAUAACUAUACAAUCUACUGUCUUUCCU",
                                            "CGGGGUGAGGUAGUAGGUUGUGUGGUUUCAGGGCAGUGAUGUUGCCCCUCGGAAGAUAACUAUACAACCUACUGCCUUCCCUG"))
  #make sure work as expected
  expect_identical(mutate.file, c("(((((.(((((((((((((((((((((.....(((...((((....)))).)))))))))))))))))))))))))))))",
                              "(((..(((.(((.(((((((((((((.........(((......)))))))))))))))).))).))).)))",
                              "(((.(((((((((((((((((((((((((((...)))))).........))))))))))))))))))))).)))",
                              "(((((.((((((((((((((((((((((((((((((.....))).)))).))).....)))))))))))))))))))))))))"))

  #make sure work as expected
  expect_identical(mutate.file, mutate.input)

})


test_that("File input for distance is correct", {

  #executable path error
  expect_error(mutate <- predict.distance(executable.path = ".\foo"))

  #missing attribute
  expect_error(mutate <- predict.distance(executable.path = ""))

  #Input file size misalignment
  expect_error(mutate <- predict.distance(executable.path = "", name = c("heyhey", "blueblue")
                                          , seq.ori = c("JJJ"), seq.alt = c("AAA")))

  #Input file type missed
  expect_error(mutate <- predict.distance(executable.path = "", name = c("heyhey", "blueblue")
                                          , seq.ori = c(1, 2), seq.alt = c("AAA", "BBB")))


})

if (TRUE) skip("Unable to test under current envirment")

test_that("RNA distance prediction is working", {

  #Load file
  distance <- predict.distance(executable.path = "",
                                      name = c("hsa-let-7a-1", "hsa-let-7a-2"),
                               seq.ori = c("(((((.(((((((((((((((((((((.....(((...((((....)))).)))))))))))))))))))))))))))))",
                                               "(((..(((.(((.(((((((((((((.........(((......)))))))))))))))).))).))).)))"),
                               seq.alt = c("(((((.(((((((((((((((((((((.....(((...((((....)))).)))))))))))))))))))))))))))))",
                                              "(((..(((.(((.(((((((((((.(((..................)))))))))))))).))).))).)))"))

  #make sure work as expected
  expect_equal(distance[1], 0)
  expect_equal(distance[2], 12)


})
