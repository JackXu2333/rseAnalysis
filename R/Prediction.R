#' Return list of predicted secondary structure
#'
#' Perform RNA structural prediction based on fasta file or
#' input of RNA name and sequence. Secondary structure prediction
#' algorithm based on command line program RNAStructure from
#' Mathews Lab, program require current path to the executable,
#' MacOS or Unix based user installed entire package is allowed to
#' omit.
#'
#' @param executable.path The path referred to RNAStructure executable (.../RNAStructure/exe) used by program
#' @param fasta.file The file path with mature sequences in fasta format
#' @param rna.name List consisted with name of the RNA sequence in seq
#' @param rna.seq List consisted with RNA sequence
#'
#' @return Returns a list of dot bracket form, ordered according to the input name sequence
#'
#' @examples
#' \dontrun{dot <- predict.Structure(executable.path = "./exe"
#'                                  , rna.name = c("hsa-let-7b", "hsa-let-7a-2")
#'                                  , rna.seq = c("UGGGAUGAGGUGGAUGUCUUUCCUA", "GAUAACUAUACAAUCUAC"))}
#'
#' @author Sijie Xu, \email{sijie.xu@mail.utoronto.ca}
#'
#' @references
#'
#' Duan, S., Mathews, D.H. and Turner, D.H. (2006).
#' Interpreting oligonucleotide microarray data to determine RNA secondary structure: application to the 3' end of Bombyx mori R2 RNA.
#' Biochemistry, 45:9819-9832.
#'
#' Wuchty, S., Fontana, W., Hofacker, I.L. and Schuster, P. (1999).
#' Complete suboptimal folding of RNA and the stability of secondary structures.
#' Biopolymers, 49:145-165.
predict.Structure <- function (executable.path = "", fasta.file = "", rna.name = c(), rna.seq = c()) {

  AllSub.path <- paste(executable.path, "/Allsub", sep = "")
  ct2dot.path <- paste(executable.path, "/ct2dot", sep = "")

  #Validate Input
  if (!file.exists(AllSub.path) || !file.exists(ct2dot.path)){
    stop("Structural prediction method unavailable")
  } else if (!missing(fasta.file) && !file.exists(fasta.file)){
    stop("fasta file unavailable")
  } else if (missing(fasta.file) && missing(rna.name) && missing(rna.seq)){
    stop("Input file unavailable")
  } else if (missing(fasta.file) && (missing(rna.name) || missing(rna.seq))){
    stop("Must input both name and sequence")
  }

  if (missing(fasta.file)){
    if (length(rna.name) != length(rna.seq)){
      stop("Name and Sequence does not share the same length")
    } else if (typeof(rna.name) != "character" ||  typeof(rna.seq) != "character"){
      stop("Invalidate name or seq type")
    }
  }

  #Create new folder to store temperary files
  dir.create("temp")

  #If fasta file exist, read them, otherwise load them
  if (!missing(fasta.file)){
    fasta <- fasta2df(fasta.file)
  } else {
    fasta <- data.frame(NAME = rna.name, SEQ = rna.seq)
  }

  #Create empty list to store the outcome
  dot.list <- rep("", nrow(fasta))

  #Put fasta file into individual files
  cat("Loading fasta files..")
  for (i in seq(nrow(fasta))){

    #Prepare file names for command line execution
    fa.filename <- paste("./temp/", i, ".fasta", sep="")
    ct.filename <- paste("./temp/", i, ".ct", sep="")
    dot.filename <- paste("./temp/", i, ".dot", sep="")

    #Generate single line fasta file
    cat(">", fasta$NAME[i], "\n", fasta$SEQ[i], "\n", sep="", file=fa.filename)

    #Predict Secondary Structure
    system(paste(AllSub.path, fa.filename, ct.filename, sep=" "))

    #Convert into the dot file, default chosen by free energy
    system(paste(ct2dot.path, ct.filename, "1", dot.filename, sep=" "))

    #Read result from the dot file
    dot <- scan(dot.filename, what="", sep="\n")[3]

    #Load dot into dot list
    dot.list[i] <- dot
  }

  unlink("./temp",recursive = TRUE)

  return(dot.list)

}



########################################################################################




#' Return list of calculated RNA distance
#'
#' Perform RNA distance calculation based on imported file or
#' input of RNA orginal sequence and alternative sequence. The RNA
#' distance calculation is based on the RNAdistance command line
#' algorithm by Walter Fontana, Ivo L Hofacker, Peter F Stadler.
#' RNA distance calculation is based on -XF parameter, which
#' indicates single alignment among each comparisons.
#' MacOS or Unix user can ignore executable.path, when RNADistance
#' is already installed in the system tools.
#'
#' @param executable.path The path referred to RNADistance executable (.../RNADistance) used by program
#' @param name List consisted with name of the RNA sequence in seq
#' @param seq.ori List consisted with RNA orginal sequence
#' @param seq.alt List consisted with RNA alternative sequence
#'
#' @return Returns a list of RNA distance, ordered according to the input name sequence
#'
#' @examples
#' \dontrun{dot <- calculate.distance(name = c("hsa-let-7b", "hsa-let-7a-2")
#' , seq.ori = c("(((((.((((((((((((((((((((((((((((((.....))).)))).))).....)))))))))))))))))))))))))",
#'           "(((((.((((((((((((((((((((((((((((((.....))).)))).))).....)))))))))))))))))))))))))")
#' , seq.alt = c("(((((.(((((((((((((((((((..(((((((((.....)))))).).....))...))))))))))))))))))))))))",
#'            "(((((.((((((((.(((((((((((((((((((((.....)))))).).))).....))))))))))).)))))))))))))"))}
#'
#' @author Sijie Xu, \email{sijie.xu@mail.utoronto.ca}
#'
#' @references
#'
#' Lorenz, Ronny and Bernhart, Stephan H. and Höner zu Siederdissen, Christian and Tafer, Hakim and Flamm, Christoph and Stadler, Peter F. and Hofacker, Ivo L.
#' ViennaRNA Package 2.0 Algorithms for Molecular Biology, 6:1 26, 2011, doi:10.1186/1748-7188-6-26
predict.distance <- function(executable.path = "", name = c(), seq.ori = c(), seq.alt = c()){

  #Validate Input
  if (missing(name) || missing(seq.ori) || missing(seq.alt)){
    stop("Input file unavailable")
  } else if (Sys.which("RNADistance") == "" &&  executable.path == ""){
    stop("Unable to find RNADistance installation")
  } else if (executable.path == ""){
    executable.path = "RNADistance"
  }

  if (length(name) != length(seq.ori) || length(name) != length(seq.alt)){
    stop("Name and Sequence does not share the same length")
  } else if (typeof(name) != "character" ||  typeof(seq.ori) != "character" ||  typeof(seq.alt) != "character"){
    stop("Invalidate input type")
  }

  #Create new folder to store temperary files
  dir.create("temp")

  #Constructe Object to store input
  fasta <- data.frame(name = name, orginal = seq.ori, alternative = seq.alt)

  #Create empty list to store the outcome
  dis.list <- rep(0, nrow(fasta))

  #Put fasta file into individual files
  cat("Loading files..")
  for (i in seq(nrow(fasta))){

    #Prepare file names for command line execution
    fa.filename <- paste("./temp/", i, ".fasta", sep="")
    dis.filename <- paste("./temp/", i, ".txt", sep="")

    #Generate single line fasta file
    cat(">", fasta$name[i], "\n", fasta$orginal[i], "\n", fasta$alternative[i], "\n", sep="", file=fa.filename)

    #Predict Secondary Structure
    system(paste(executable.path, "-Xf < ", fa.filename, ">", dis.filename))

    #Read result from the dot file
    dis <- scan(dis.filename, what="", sep="\n")[2]

    #Load dot into dot list
    dis.list[i] <- substr(dis, 4, nchar(dis) - 2)
  }

  unlink("./temp",recursive = TRUE)

  return(strtoi(dis.list))
}


