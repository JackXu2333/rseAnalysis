# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Cmd + Shift + B'
#   Check Package:             'Cmd + Shift + E'
#   Test Package:              'Cmd + Shift + T'

# Return a dataframe object with gene mutation information
vcf2df <- function(filepath) {

  #Validate file path
  if (!file.exists(filepath) && dir.exists(filepath)){ #path point to directory
    stop(sprintf("File path %s is a directory, please name path to a file", filepath))

  } else if (!file.exists(filepath)) { #path not found
    stop(sprintf("File path %s does not exist", filepath))

  } else if (file_ext(filepath) != "vcf"){ #file not in vcf
    stop("Select a validate vcf file")
  }

  #Import the vcf file
  vcf_source <- vcfR::read.vcfR(filepath)

  #Combine and output vcf base and mutation info
  vcf <- cbind(as.data.frame(vcfR::getFIX(vcf_source)), vcfR::INFO2df(vcf_source))

  return(vcf)

}

# Return a dataframe object with fasta information
fasta2df <- function(filepath) {

  #Validate file path
  if (!file.exists(filepath) && dir.exists(filepath)){ #path point to directory
    stop(sprintf("File path %s is a directory, please name path to a file", filepath))

  } else if (!file.exists(filepath)) { #path not found
    stop(sprintf("File path %s does not exist", filepath))

  } else if (file_ext(filepath) != "fasta"){ #file not in vcf
    stop("Select a validate fasta file")
  }

  #Import the fasta file
  fasta_source <- Biostrings::readRNAStringSet("./GRCh38_mRNA.fa")

  #Format FASTA file
  fasta <- data.frame(name = names(fasta_source), sequence = paste(fasta_source))

  return(fasta)

}
