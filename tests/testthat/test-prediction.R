context("secondary structure and distance prediction")
library(rseAnalysis)

test_that("File input for secondary structure is correct", {

  #executable path error
  expect_error(mutate <- predict.Structure(executable.path = "./foo", fasta.file = "../../inst/extdata/test.fasta"))

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
  mutate.file <- predict.Structure(executable.path = "../../inst/extdata/exe", fasta.file = "../../inst/extdata/test.fasta")

  #Load via file
  mutate.input <- predict.Structure(executable.path = "../../inst/extdata/exe",
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
                                          , struct.ori = c("JJJ"), struct.alt = c("AAA")))

  #Input file type missed
  expect_error(mutate <- predict.distance(executable.path = "", name = c("heyhey", "blueblue")
                                          , struct.ori = c(1, 2), struct.alt = c("AAA", "BBB")))


})


test_that("RNA distance prediction is working", {

  #Load file
  distance <- predict.distance(executable.path = "",
                                      name = c("hsa-let-7a-1", "hsa-let-7a-2"),
                               struct.ori = c("(((((.((((((((((((((((((((((((((((((.....))).)))).))).....)))))))))))))))))))))))))",
                                               "(((((.((((((((((((((((((((((((((((((.....))).)))).))).....)))))))))))))))))))))))))"),
                               struct.alt = c("(((((.(((((((((((((((((((..(((((((((.....)))))).).....))...))))))))))))))))))))))))",
                                              "(((((.((((((((.(((((((((((((((((((((.....)))))).).))).....))))))))))).)))))))))))))"))

  #make sure work as expected
  expect_equal(distance[1], 10)
  expect_equal(distance[2], 4)


})
