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


#Transcript DNA into RNA
DNA2RNA <- function(DNA.Seq){

  #Break into small pieces
  seq <- seqinr::s2c(DNA.Seq)

  #Transcript the DNA
  seq[which(seq == "T")] <- "U"

  return(seq)

}

#Preform validation checking on fasta and vcf files
mRNA.validate <- function(fasta, vcf){

  #Counnting the number of matches and the matched itemed
  match.count=0
  error.count=0

  fasta.matched <- data.frame(name=character()
                                , sequence=character()
                                , startPos=integer()
                                , pointPos=integer()
                                , mutationRNA=character()
                                , referenceRNA=character()
                                , realization=character()
                                , doesMatch=boolean())

  for (i in seq(nrow(mRNA))){#Iterate through mRNA

    #Locate mutation in mRNA
    mutation <- vcf[ mRNA$ENDPOS[i] >= vcf$POS & vcf$POS >= mRNA$STAPOS[i], ]

    if (nrow(mutation) > 1){#If mutation is found

      for (p in seq(nrow(mutation))){# Iterate through each mutation

        #Find out sequence locating at given mutation location
        orgin <- find(mRNA$SEQ[i]
                      , mRNA$STAPOS[i]
                      , mutation$POS[p]
                      , mutation$REF[p])

        #Record the matching data
        newMatch <-data.frame(name=paste0(mRNA$NAME[i], ',', mutation$OCCURRENCE[p])
                              , sequence=mRNA$SEQ[i]
                              , startPos=mRNA$STAPOS[i]
                              , pointPos=mutation$POS[p]
                              , mutationRNA=DNA_Transcript(mutation$ALT[p])
                              , referenceRNA=DNA_Transcript(mutation$REF[p])
                              , realization=orgin)

        #Compare against the mutation reference
        if(orgin == newMatch$referenceRNA){
          match_count=match_count + 1

          newMatch$doesMatch = TRUE
          newMatch$mutated.sequence = mutate(newMatch$sequence, newMatch$startPos
                                             , newMatch$pointPos, newMatch$referenceRNA
                                             , newMatch$mutationRNA)
        } else {
          error_count=error_count + 1

          newMatch$doesMatch = FALSE
          newMatch$mutated.sequence = 'NO MATCH'
        }

        fasta.matched <- as.data.frame(rbind(fasta.matched, newMatch))

      }
    }
  }

  cat(sprintf("The matching rate of the dataset is %s", match_count/(match_count + error_count)))

  return()

}


#Preform mutation on mRNA using fasta and vcf files
mRNA.mutate <- function(){

}
