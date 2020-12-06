#' Return the transcript of input DNA
#'
#' Utility function that return the complimentary
#' sequence of transcript RNA of input DNA sequence.
#' This function will validate the DNA sequence.
#'
#' @param DNA.Seq A deoxyribonucleic acid sequence
#'
#' @return Return the complementary sequence of DNA.Seq's transcript RNA
#'
#' @examples
#' (RNA <- DNA2RNA("TGGGATGAGGTGGATGTTTCCTA"))
#'
#' @author Sijie Xu, \email{sijie.xu@mail.utoronto.ca}
#'
#' @export

DNA2RNA <- function(DNA.Seq){

  DNAbase <- c("T", "A", "C", "G")

  if (typeof(DNA.Seq) != "character"){
    stop("Invalidated input type")
  }

  #Validate DNASeq
  DNASeqList <- seqinr::s2c(DNA.Seq)

  if (length(DNASeqList) == 0){
    stop("Invalidate DNA sequence found")
  }

  for (nucltide in DNASeqList) {
    if (!nucltide %in% DNAbase){
      stop("Invalidate DNA sequence found")
    }
  }

  #Transcript the DNA
  DNASeqList[which(DNASeqList == "T")] <- "U"

  return(gsub(", ", "", toString(DNASeqList)))

}



########################################################################################



#' Apply DNA mutation to RNA and validate the file
#'
#' Calculate the mutated series of RNA from RNA fasta, bed file
#' , and DNA vcf file, both should be result file from vcf2df
#' and fasta2df file, function gives validation based on fasta
#' and vcf comparisons, user should user validation result to make
#' sure their vcf and fasta file are aligned under same version.
#' Algorithm inspired by the mRNA mutation function, a utility function
#' created by Sijie Xu and Zhiwen Tan, the updated method here is a more generalized
#' and input validated method, more checking are embeded to insure the
#' mutation is applied correctly.
#'
#' @param fasta RNA sequence to be matched read using fasta2df
#' @param vcf mutation information read using vcf2df function
#' @param bed bed information read using bed2df function
#'
#' @return Return the RNA.mutated data frame that includes all information
#' about the RNA and its information
#' \itemize{
#'   \item NAME - Corresponding of the RNA sequence
#'   \item MATCH - Boolean value indicating whether the mutation can be applied,
#'         FALSE value indicate the sequence given has inconsistent bases comparing to the
#'         reference base in VCF file
#'   \item STAPOS - The start location of the sequence
#'   \item ENDPOS - The end location of the sequence
#'   \item DIR - The direction of the sequence
#'   \item SEQ - The original RNA sequence
#'   \item CUR.REF - The current RNA base(s) on the mutation point
#'   \item MUT.REF - The reference base(s) on the mutation point
#'   \item MUT.ALT - The alternative base(s) appeared in during mutation
#'   \item MUT.SEQ - The mutated RNA sequence when MATCH is TRUE, otherwise empty string
#'   \item MUT.POS - The start position of the mutation
#' }
#'
#' @examples
#' mutated <- RNA.validate(
#'               fasta = fasta2df(system.file("extdata", "test.fasta", package = "rseAnalysis")),
#'               vcf = vcf2df(system.file("extdata", "test.vcf", package = "rseAnalysis")),
#'               bed = bed2df(system.file("extdata", "test.bed", package = "rseAnalysis")))
#'
#' @author Sijie Xu, \email{sijie.xu@mail.utoronto.ca}
#'
#' @references Xu, S. Tan, Z., BCB430 "Analysis.zip/File Validation.rmd" <https://github.com/Deemolotus/BCB330Y-and-BCB430Y>
#'
#' @export

RNA.validate <- function(fasta, vcf, bed){

  #validate input type
  if (typeof(fasta) != "list" || typeof(vcf) != "list" || typeof(bed) != "list"){
    stop("invalidated fasta type")
  }

  #Validate fasta input
  if (any(!c("NAME", "SEQ") %in% names(fasta))) {
    stop("corrupted fasta input")
  }

  #Validate fasta input
  if (any(!c("CHROM", "STAPOS", "ENDPOS", "DIR", "NAME") %in% names(bed))) {
    stop("corrupted bed input")
  }

  #Validate vcf input
  if (any(!c("CHROM", "POS", "REF", "ALT") %in% names(vcf))) {
    stop("corrupted vcf input")
  }

  #Validate input size
  if (nrow(fasta) == 0 || nrow(bed) == 0  || nrow(vcf) == 0 ){
    stop("empty data found")
  }

  #Merge bed and fasta file to retrieve RNA sequence and location
  RNA <- merge(bed, fasta, by="NAME")

  #Private function to find the reference sequence on the target location
  RNA.find <- function(sequence, startPos, pointPos, orgin){

    #Calculate the relative distance of position of the mutation from start
    mutatePos <- pointPos - startPos

    return(substr(sequence, mutatePos, mutatePos + nchar(orgin) - 1))

  }

  #Private function return mutated RNA using fasta and vcf files
  RNA.mutate <- function (sequence, startPos, pointPos, orgin, mutation){

    LHS <- substr(sequence, 0, pointPos - startPos - 1)
    RHS <- substr(sequence, pointPos - startPos + nchar(orgin), nchar(sequence))

    return (paste0(LHS, mutation, RHS))

  }

  #Counnting the number of matches and the matched itemed
  match.count = 0
  error.count = 0

  #Initialized the empty list for data input
  RNA.matched <- data.frame(NAME=character(), MATCH=logical()
                            , STAPOS=integer() , ENDPOS=integer()
                            , DIR=character(), SEQ=character()
                            , CUR.REF=character(), MUT.REF=character()
                            , MUT.ALT=character(), MUT.SEQ=character()
                            , MUT.POS=character(), MUT.ID=character())

  for (i in seq(nrow(RNA))){#Iterate through mRNA

        #Locate mutation in mRNA
        mutation <- vcf[which(RNA$ENDPOS[i] >= vcf$POS & vcf$POS >= RNA$STAPOS[i] & vcf$CHROM == RNA$CHROM[i]), ]

    if (nrow(mutation) >= 1){#If mutation is found

      for (p in seq(nrow(mutation))){# Iterate through each mutation

        #Find sequence according to direction
        if(RNA$DIR[i] == "-") {current.seq <- Biostrings::reverse(RNA$SEQ[i])}
        else {current.seq <- RNA$SEQ[i]}

        #Find out RNA sequence locating at given mutation location
        current.ref <- RNA.find(RNA$SEQ[i]
                      , RNA$STAPOS[i]
                      , mutation$POS[p]
                      , mutation$REF[p])

        #Record the matching data
        newMatch <- data.frame(NAME=RNA$NAME[i], MATCH=TRUE
                                  , STAPOS=RNA$STAPOS[i], ENDPOS=RNA$ENDPOS[i]
                                  , DIR=RNA$DIR[i], SEQ=current.seq
                                  , CUR.REF=current.ref, MUT.REF=DNA2RNA(mutation$REF[p])
                                  , MUT.ALT=DNA2RNA(mutation$ALT[p]), MUT.SEQ=""
                                  , MUT.POS=mutation$POS[p], MUT.ID=mutation$ID[p])

        #Compare against the RNA mutation reference
        if(current.ref == newMatch$MUT.REF){ #If Match is found

          #Increment the match count
          match.count=match.count + 1

          #Update the newMatch
          newMatch$MATCH = TRUE
          newMatch$MUT.SEQ = RNA.mutate(newMatch$SEQ, newMatch$STAPOS
                                             , newMatch$MUT.POS, current.ref
                                             , newMatch$MUT.ALT)
        } else {

          #Increment the error count
          error.count=error.count + 1

          #Update the newMatch
          newMatch$MATCH = FALSE
        }

        RNA.matched <- rbind(RNA.matched, newMatch)

      }
    }
  }

  #Print the success rate of the matching
  cat(sprintf("The matching rate of the dataset is %s \n", match.count/(match.count + error.count)))

  return(unique(RNA.matched))

}

