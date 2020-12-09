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
#' \dontrun{dot <- predictStructure(executable.path = "./exe"
#'                                  , rna.name = c("hsa-let-7b", "hsa-let-7a-2")
#'                                  , rna.seq = c("UGGGAUGAGGUGGAUGUCUUUCCUA", "GAUAACUAUACAAUCUAC"))}
#'
#' @author Sijie Xu, \email{sijie.xu@mail.utoronto.ca}
#'
#' @references
#' {Duan, S., Mathews, D.H. and Turner, D.H. (2006).
#' Interpreting oligonucleotide microarray data to determine RNA secondary structure: application to the 3' end of Bombyx mori R2 RNA.
#' Biochemistry, 45:9819-9832.}
#'
#' {Wuchty, S., Fontana, W., Hofacker, I.L. and Schuster, P. (1999).
#' Complete suboptimal folding of RNA and the stability of secondary structures.
#' Biopolymers, 49:145-165.}
#'
#' @export

predictStructure <- function (executable.path = "", fasta.file = "", rna.name = c(), rna.seq = c()) {

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

  # Add name back to the list
  names(dot.list) <- rna.name

  return(dot.list)

}



########################################################################################



#' (Internal Function) Return list of calculated RNA distance using RNADistance
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
#' @param struct.ori List consisted with RNA orginal sequence
#' @param struct.alt List consisted with RNA alternative sequence
#'
#' @return Returns a list of RNA distance, ordered according to the input name sequence
#'
#' @author Sijie Xu, \email{sijie.xu@mail.utoronto.ca}
#'
#' @references
#' {Lorenz, Ronny and Bernhart, Stephan H. and Höner zu Siederdissen, Christian and Tafer, Hakim and Flamm, Christoph and Stadler, Peter F. and Hofacker, Ivo L.
#' ViennaRNA Package 2.0 Algorithms for Molecular Biology, 6:1 26, 2011, doi:10.1186/1748-7188-6-26}
#' {Steipe B., ABC R Project, A Bioinformatics Course: Applied Bioinformatics http://steipe.biochemistry.utoronto.ca/abc/index.php/Bioinformatics_Main_Page}
#'

predictDistance_RNADis <- function(executable.path = "", name = c(), struct.ori = c(), struct.alt = c()){

  #Helper Function that draw a progress bar in the console
  #Referenced from ABC project (.utility 4.07) by Professor Steipe B. at University of Toronto
  pBar <- function(i, l, nCh = 50) {
    # i: the current iteration
    # l: the total number of iterations
    # nCh: width of the progress bar
    ticks <- round(seq(1, l-1, length.out = nCh))
    if (i < l) {
      if (any(i == ticks)) {
        p <- which(i == ticks)[1]  # use only first, in case there are ties
        p1 <- paste(rep("#", p), collapse = "")
        p2 <- paste(rep("-", nCh - p), collapse = "")
        cat(sprintf("\r|%s%s|", p1, p2))
        flush.console()
      }
    }
    else { # done
      cat("\n")
    }
  }

  #Create new folder to store temperary files
  dir.create("temp")

  #Constructe Object to store input
  fasta <- data.frame(name = name, orginal = struct.ori, alternative = struct.alt)

  #Create empty list to store the outcome
  dis.list <- rep(0, nrow(fasta))

  #Put fasta file into individual files
  cat("Loading files..")
  for (i in seq(nrow(fasta))){

    pBar(i, nrow(fasta))

    #Prepare file names for command line execution
    fa.filename <- paste("./temp/", i, ".fasta", sep="")
    dis.filename <- paste("./temp/", i, ".txt", sep="")

    #Generate single line fasta file
    cat(">", fasta$name[i], "\n", fasta$orginal[i], "\n", fasta$alternative[i], "\n", sep="", file=fa.filename)

    #Predict Secondary Structure
    system(paste(executable.path, "-Xf < ", fa.filename, ">", dis.filename))

    #Read result from the dot file
    dis <- suppressWarnings(suppressMessages(scan(dis.filename, what="", sep="\n", quiet = TRUE)))[2]

    #Load dot into dot list
    dis.list[i] <- substr(dis, 4, nchar(dis) - 2)
  }

  unlink("./temp",recursive = TRUE)

  return(strtoi(dis.list))
}



########################################################################################



#' (Internal Function) Return list of calculated RNA distance using gscVisualizer
#'
#' Perform RNA distance calculation based on imported file or
#' input of RNA original sequence and alternative sequence. The RNA
#' distance calculation is based on the gscVisualizer R package
#' algorithm by Zhiwen Tan, a BCB student from UofT.
#'
#' @param name List consisted with name of the RNA sequence in seq
#' @param struct.ori List consisted with RNA orginal sequence
#' @param struct.alt List consisted with RNA alternative sequence
#'
#' @return Returns a list of RNA distance, ordered according to the input name sequence
#'
#' @author Sijie Xu, \email{sijie.xu@mail.utoronto.ca}
#'
#' @references
#' {Z Tan (2020). gscVisualizer: Visualize The Difference among Gene Sequences. R package version 0.1.0.}
#'
#' @importFrom gscVisualizer dotComp

predictDistance_gsc <- function(name = c(), struct.ori = c(), struct.alt = c()){

  #Constructe Object to store input
  fasta <- data.frame(name = name, orginal = struct.ori, alternative = struct.alt)

  #Create empty list to store the outcome
  dis.list <- rep(0, nrow(fasta))

  #Put fasta file into individual files
  cat("Loading files..")
  for (i in seq(nrow(fasta))){

    dis.list[i] <- gscVisualizer::dotComp(fasta$orginal[i], fasta$alternative[i])

  }

  return(strtoi(dis.list))
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
#' @param struct.ori List consisted with RNA orginal sequence
#' @param struct.alt List consisted with RNA alternative sequence
#' @param method Options to use "RNADis" (RNADistance) or "gsc" gscVisualizer package for RNAdistance calculation, default value is "gsc"
#'
#' @return Returns a list of RNA distance, ordered according to the input name sequence
#'
#' @examples
#' dot <- predictDistance(name = c("hsa-let-7b", "hsa-let-7a-2")
#' , struct.ori = c("(((((.((((((((((((((((((((((((((((((.....))).)))).))).....)))))))))))))))))))))))))",
#'           "(((((.((((((((((((((((((((((((((((((.....))).)))).))).....)))))))))))))))))))))))))")
#' , struct.alt = c("(((((.(((((((((((((((((((..(((((((((.....)))))).).....))...))))))))))))))))))))))))",
#'            "(((((.((((((((.(((((((((((((((((((((.....)))))).).))).....))))))))))).)))))))))))))"))
#'
#' @author Sijie Xu, \email{sijie.xu@mail.utoronto.ca}
#'
#' @references
#' {Lorenz, Ronny and Bernhart, Stephan H. and Höner zu Siederdissen, Christian and Tafer, Hakim and Flamm, Christoph and Stadler, Peter F. and Hofacker, Ivo L.
#' ViennaRNA Package 2.0 Algorithms for Molecular Biology, 6:1 26, 2011, doi:10.1186/1748-7188-6-26}
#'
#' @export
#'
#' @importFrom gscVisualizer dotComp

predictDistance <- function(executable.path = "", name = c(), struct.ori = c(), struct.alt = c(), method = "gsc"){

  #Validate RNA input
  if (missing(name) || missing(struct.ori) || missing(struct.alt)){
    stop("Input file unavailable")
  } else if (typeof(name) != "character" || typeof(struct.ori) != "character" || typeof(struct.alt) != "character"){
    stop("Input file should be list of strings")
  } else if (length(name) != length(struct.ori) || length(name) != length(struct.alt)){
    stop("Name and Sequence does not share the same length")
  }

  #Validate method selection
  if (!(method %in% c("gsc", "RNADis"))){
    stop("Please select validate prediction algorithm by input gsc or RNADis")

  } else if (method == "RNADis"){ #If method selection is RNADistance

    #Validate RNADistance
    if (Sys.which("RNADistance") == "" &&  executable.path == ""){
      stop("Unable to find RNADistance executable")
    } else if (executable.path == ""){
      executable.path = "RNADistance"
    }

    return(predictDistance_RNADis(executable.path, name, struct.ori, struct.alt))

  } else { #If method selection is gscVisualizer

    return(predictDistance_gsc(name, struct.ori, struct.alt))

  }
}


# [END]
