## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----input,  warning=FALSE----------------------------------------------------

#Source library
library(rseAnalysis)

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


## ----structure, message=FALSE, warning=FALSE----------------------------------

#Filter out the one that is misaligned, choose the first 200, or it takes forever to run
RNA.mutated <- subset(RNA.mutated, MATCH)[1:200,]

#Perform structure prediction on orginal sequence
struct.ori <- suppressMessages(predict.Structure(executable.path = "../inst/extdata/exe"
                                   , rna.name = RNA.mutated$NAME, rna.seq = RNA.mutated$SEQ))

#Perform structure prediction on mutated sequence
struct.alt <- suppressMessages(predict.Structure(executable.path = "../inst/extdata/exe"
                                   , rna.name = RNA.mutated$NAME, rna.seq = RNA.mutated$MUT.SEQ))


## ----distance, message=FALSE, warning=FALSE-----------------------------------

#Run prediction
RNA.distance <- predict.distance(executable.path = "", name = RNA.mutated$NAME, 
                  struct.ori = struct.ori, struct.alt = struct.alt)


## ----analysis-----------------------------------------------------------------

#Load expression data

expression <- read.csv(system.file("extdata", "test.csv", package = "rseAnalysis"), header = TRUE)

#Use only standardize read
expression <- subset(expression, Read.Type == "reads_per_million_miRNA_mapped")[1:200, ]

Analysis.DISEXP(dis.name = RNA.mutated$NAME, dis.distance = RNA.distance, 
                exp.tumor <- expression$Sample, exp.sample <- expression$Normal, method = "linear", seperate = TRUE)


## -----------------------------------------------------------------------------
sessionInfo()

