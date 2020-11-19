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

#Perform secondary structure prediction for sequence
predict.Structure <- function (execuable.path, fasta.file, name, seq) {

  AllSub.path <- paste(execuable.path, "Allsub", sep = "")
  ct2dot.path <- paste(execuable.path, "ct2dot", sep = "")

  #Validate Input
  if (!file.exists(AllSub.path) || !file.exists(ct2dot.path)){
    stop("Structural prediction method unavailable")
  } else if (!file.exists(fasta.file)){
    stop("fasta file unavailable")
  } else if (missing(fasta.file) && missing(name) && missing(seq)){
    stop("Input file unavailable")
  } else if (missing(fasta.file) && (missing(name) || missing(seq))){
    stop("Must input both name and sequence")
  }

  if (length(name) != length(seq)){
    stop("Name and Sequence does not share the same length")
  }

  #Create new folder to store temperary files
  dir.create("temp")

  #If fasta file exist, read them, otherwise load them
  if (!missing(fasta.file)){
    fasta <- fasta2df(fasta.file)
  } else {
    fasta <- data.frame(name = name, sequence = seq)
  }

  #Create empty list to store the outcome
  dot.list <- vector("list", length = nrow(fasta))

  #Put fasta file into individual files
  cat("Loading fasta files..")
  for (i in seq(nrow(fasta))){

    #Prepare file names for command line execution
    fa.filename <- paste("./temp/", fasta$name[1], ".fasta", sep="")
    ct.filename <- paste("./temp/", fasta$name[1], ".ct", sep="")
    dot.filename <- paste("./temp/", fasta$name[1], ".dot", sep="")

    #Generate single line fasta file
    cat(">", fasta$name[1], "\n", fasta$sequence[1], "\n", sep="", file=fa.filename)

    #Predict Secondary Structure
    invisible(capture.output(system(paste(AllSub.path, "1", fa.filename, ct.filename, sep=" "))))

    #Convert into the dot file, default chosen by free energy
    system(paste(ct2dot.path, ct.filename, dot.filename, sep=" "))

    #Read result from the dot file
    dot <- read.table(file=dot.filename, header=FALSE, sep=" ")[3]

    #Load dot into dot list
    dot.list[i] <- dot
  }

  return(dot.list)

}

#Perform RND distance calculation
calculate.distance <- function(){

}


