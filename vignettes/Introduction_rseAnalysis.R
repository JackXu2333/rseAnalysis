## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  out.width = "100%"
)

## ----input,  warning=FALSE----------------------------------------------------

#Source library
library(rseAnalysis)
library(ggplot2)

#Load sample data file
vcf <- rseAnalysis::vcf2df(system.file("extdata", "hsa_GRCh37.vcf", package = "rseAnalysis"))
fasta <- rseAnalysis::fasta2df(system.file("extdata", "hsa_GRCh37.fasta", package = "rseAnalysis"))
bed <- rseAnalysis::bed2df(system.file("extdata", "hsa_GRCh37.bed", package = "rseAnalysis"))

#Inspect the imported file
head(vcf)
head(fasta)
head(bed)


## ----mutate,  warning=FALSE---------------------------------------------------

#Mutate RNA using mutation from vcf files
RNA.mutated <- RNA.validate(fasta = fasta, 
             vcf = vcf, 
             bed = bed)


## ----structure----------------------------------------------------------------

# ================== Sample code for RNA secondary structure prediction ==========================
#
#  struct.ori <- suppressMessages(predictStructure(executable.path = "../inst/extdata/exe"
#                                   , rna.name = RNA.mutated$NAME, rna.seq = RNA.mutated$SEQ))
#  struct.alt <- suppressMessages(predictStructure(executable.path = "../inst/extdata/exe"
#                                   , rna.name = RNA.mutated$NAME, rna.seq = RNA.mutated$MUT.SEQ))

# Read prerun result from the predictStructure
RNA.mutated <- subset(RNA.mutated, MATCH)[1:200,]
struct.ori <- read.csv(system.file("extdata", "vignetteSampleORI.csv", package = "rseAnalysis"))
struct.alt <- read.csv(system.file("extdata", "vignetteSampleALT.csv", package = "rseAnalysis"))

head(struct.ori)
head(struct.alt)


## ----distance, message=FALSE, warning=FALSE-----------------------------------

#Run prediction
RNA.distance <- predictDistance(name = RNA.mutated$NAME
                                 , struct.ori = struct.ori$struct.ori
                                 , struct.alt = struct.alt$struct.alt
                                 , method = "gsc")


## ----analysis,  warning=FALSE-------------------------------------------------

#Load expression data

expression <- read.csv(system.file("extdata", "test.csv", package = "rseAnalysis"), header = TRUE)

#Use only standardize read
expression <- subset(expression, Read.Type == "reads_per_million_miRNA_mapped")[1:200, ]

result <- Analysis.DISEXP(dis.name = RNA.mutated$NAME, dis.distance = RNA.distance, 
                exp.tumor = expression$Sample, exp.sample = expression$Normal, method = "linear", showPlot = FALSE)

#Display statistical result
result$stats

#Display images from result
#result$plots


## -----------------------------------------------------------------------------
sessionInfo()

