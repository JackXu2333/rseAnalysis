#' Return a dataframe object with gene mutation information
#'
#' Perform vcf file retrieving using read.vcfR function, result stored
#' in standarlize data.frame for easy usage.
#'
#' @param filepath The path directed to vcf file on current working directory
#'
#' @return Return a data.frame object contain vcf document loaded
#' \itemize{
#'   \item CHROM - Located in chromosome CHROM
#'   \item POS - The start position of the mutation
#'   \item ID - ID of the mutation (if available)
#'   \item REF - The reference base(s) on the mutation point
#'   \item ALT - The alternative base(s) appeared in during mutation
#'   \item CONSEQUENCE - Record the effect of mutation under biological studies
#'   \item OCCURRENCE - Record the associated cancer and its distribution
#'   \item affected_donors - Record the total number of donor affected
#'   \item project_count - Number of project associated with current studies
#' }
#'
#'
#' @examples
#' filePath <- system.file("extdata", "test.vcf", package = "rseAnalysis")
#' vcf <- vcf2df(filePath)
#'
#' @author Sijie Xu, \email{sijie.xu@mail.utoronto.ca}
#'
#' @references
#' Knaus BJ, Grünwald NJ (2017). “VCFR: a package to manipulate and visualize variant call format data in R." _Molecular Ecology Resources_, *17*(1), 44-53. ISSN 757,
#' <URL: http://dx.doi.org/10.1111/1755-0998.12549>.
#'
#' @export
#' @importFrom tools file_ext
#' @importFrom vcfR read.vcfR

vcf2df <- function(filepath) {

  #Validate file path
  if (!file.exists(filepath) && dir.exists(filepath)){ #path point to directory
    stop(sprintf("File path %s is a directory, please name path to a file", filepath))

  } else if (!file.exists(filepath)) { #path not found
    stop(sprintf("File path %s does not exist", filepath))

  } else if (tools::file_ext(filepath) != "vcf"){ #file not in vcf
    stop("Select a validate vcf file")
  }

  #Import the vcf file
  vcf_source <- vcfR::read.vcfR(filepath)

  #Combine and output vcf base and mutation info
  vcf <- cbind(as.data.frame(vcfR::getFIX(vcf_source)), vcfR::INFO2df(vcf_source))

  vcf$POS <- strtoi(vcf$POS)

  return(subset(vcf, select = c("CHROM", "POS", "ID", "REF", "ALT", "CONSEQUENCE", "OCCURRENCE", "affected_donors", "project_count")))

}



########################################################################################



#' Return a dataframe object with fasta information
#'
#' Perform fasta file retrieving using readRNAStringSet function from Biostrings, result stored
#' in standarlize data.frame for easy usage.
#'
#' @param filepath The path directed to fasta file on current working directory
#'
#' @return Return a data.frame object contain fasta document loaded
#' \itemize{
#'   \item NAME - Corresponding of the RNA sequence
#'   \item SEQ - The original RNA sequence
#' }
#'
#'
#' @examples
#' filePath <- system.file("extdata", "test.fasta", package = "rseAnalysis")
#' fasta <- fasta2df(filePath)
#'
#' @author Sijie Xu, \email{sijie.xu@mail.utoronto.ca}
#'
#' @references H. Pagès, P. Aboyoun, R. Gentleman and S. DebRoy (2020). Biostrings: Efficient manipulation of biological strings. R package version 2.56.0.
#'
#' @export
#' @importFrom tools file_ext
#' @importFrom Biostrings readRNAStringSet

fasta2df <- function(filepath) {

  #Validate file path
  if (!file.exists(filepath) && dir.exists(filepath)){ #path point to directory
    stop(sprintf("File path %s is a directory, please name path to a file", filepath))

  } else if (!file.exists(filepath)) { #path not found
    stop(sprintf("File path %s does not exist", filepath))

  } else if (tools::file_ext(filepath) != "fasta"){ #file not in vcf
    stop("Select a validate fasta file")
  }

  #Import the fasta file
  fasta_source <- Biostrings::readRNAStringSet(filepath)

  #Format FASTA file
  fasta <- data.frame(NAME = sub("\\ .*", "", names(fasta_source)), SEQ = paste(fasta_source))

  return(fasta)

}



########################################################################################



#' Return a dataframe object with RNA location information
#'
#' Perform bed file retrieving using result stored
#' in standardized list for easy usage.
#'
#' @param filepath The path directed to bed file on current working directory
#'
#' @return Return a data.frame object contain bed document loaded
#' \itemize{
#'   \item CHROM - Located in chromosome CHROM
#'   \item STAPOS - The start location of the sequence
#'   \item ENDPOS - The end location of the sequence
#'   \item DIR - The direction of the sequence
#'   \item TYPE - The type of the sequence
#'   \item ID - ID of the RNA (if avilable)
#'   \item ALIAS - The alias of the sequence
#'   \item NAME - The name of the sequence
#' }
#'
#'
#' @examples
#' filePath <- system.file("extdata", "test.bed", package = "rseAnalysis")
#' bed <- bed2df(filePath)
#'
#' @author Sijie Xu, \email{sijie.xu@mail.utoronto.ca}
#'
#' @export
#' @importFrom tidyr separate
#' @importFrom utils read.table

bed2df <- function(filepath){

  #Validate file path
  if (!file.exists(filepath) && dir.exists(filepath)){ #path point to directory
    stop(sprintf("File path %s is a directory, please name path to a file", filepath))

  } else if (!file.exists(filepath)) { #path not found
    stop(sprintf("File path %s does not exist", filepath))

  } else if (tools::file_ext(filepath) != "bed"){ #file not in vcf
    stop("Select a validate bed file")
  }

  #Read entire file, separate the other info column
  bed <- utils::read.table(filepath, header = FALSE, sep="\t",stringsAsFactors=FALSE, quote="")
  bed <- tidyr::separate(data = bed, col = "V10", into = c("ID", "Alias", "Name"), sep = ";")

  #Rename the tables
  names(bed) <- c("CHROM", "STAPOS", "ENDPOS", "ID2", "NULL", "DIR", "NULL2", "TYPE", "NULL3", "ID", "ALIAS", "NAME")

  #Remove the header info among the lines
  bed$ID <- substr(bed$ID, 4, nchar(bed$ID))
  bed$ALIAS <- substr(bed$ALIAS, 7, nchar(bed$ALIAS))
  bed$NAME <- substr(bed$NAME, 6, nchar(bed$NAME))
  bed$CHROM <- as.character(bed$CHROM)

  #Return only the useful lines
  return(subset(bed, select = c("CHROM", "STAPOS", "ENDPOS", "DIR", "TYPE", "ID", "ALIAS", "NAME")))
}
