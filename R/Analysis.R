#' Create a scatter plot from RNA distance
#'
#' Utility function that return the complimentary
#' sequence of transcript RNA of input DNA sequence.
#' This function will validate the DNA sequence.
#'
#' @param DNA.Seq A deoxyribonucleic acid sequence
#'
#' @return Return the complementary sequence of DNA.Seq's transcript RNA
#'
#' @exampless
#' \dontrun{RNA <- DNA2RNA("TGGGATGAGGTGGATGTTTCCTA")}
#'
#' @author Sijie Xu, \email{sijie.xu@mail.utoronto.ca}
#'
#' @import dylpr
#' @import ggplot
Analysis.DISEXP <- function(dis.name, dis.distance, exp.name, exp.tumor, exp.sample, method = "linear"){

  #Validate input
  if (typeof(dis.name) != "list" || typeof(dis.distance) != "list" || typeof(exp.name) != "list"
      || typeof(exp.tumor) != "list" || typeof(exp.sample) != "list"){
    stop("Input must be packaged in list")
  } else if (length(dis.name) != length(dis.distance)){
    stop("RNA Distance data does not share the same size")
  } else if (length(exp.name) != length(exp.tumor) || length(exp.tumor) != length(exp.sample) ||
             length(exp.name) != length(exp.sample)){
    stop("Gene expression data does not share the same size")
  } else if (!method %in% c("linear", "poisson")){
    stop("select method from linear, possion")
  }

  #Constructs the distance and expression data frame
  dis <- data.frame(name = dis.name, distance = dis.distance)
  exp <- data.frame(name = exp.name, read_tumor = exp.tumor, read_sample = exp.sample)

  disexp.df <- merge(dis, exp, by="name")

  #Check to see of data point matches each another
  if (nrow(disexp.df) == 0){
    stop("RNA distance and gene expression file can not match")
  }

  #Calculate the difference in read from sample to tumor
  disexp.df$read_diff <- abs(disexp.df$read_sample - disexp.df$read_tumor)

  #Perform correlation analysis based on the difference gene expression and RNA distance
  #Plot the diagonise plots
  par(mfrow=c(1,2))

  #Scatter plot with distance as x and difference in read being y
  ggplot2::ggplot(disexp.df, ggplot2::aes(x=distance, y=read_diff)) +
    ggplot2::geom_point() +
    ggplot2::geom_smooth(method=lm)

  #Boxplot checking outliars on read_diff
  ggplot2::ggplot(disexp.df, ggplot2::aes(x=name, y=read_diff)) +
    ggplot2::geom_boxplot(outlier.colour="red", outlier.shape=8,
                 outlier.size=4)

  #Boxplot checking outliars on distance
  ggplot2::ggplot(disexp.df, ggplot2::aes(x=name, y=distance)) +
    ggplot2::geom_boxplot(outlier.colour="red", outlier.shape=8,
                 outlier.size=4)

  #Density plot on distance by different RNA
  ggplot2::ggplot(disexp.df, ggplot2::aes(x=distance, color=name)) +
    ggplot2::geom_density()

  if (method = "linear"){
    #Perform linear modeling on the object
    lm <- lm(read_diff ~ distance, data=disexp.df)

    #Calculate result from linear modeling
    disexp.attr <- data.frame(corrolation = cor(disexp.df$distance, disexp.df$read_diff),
                            p_value = summary(lm)$coefficients[,4][2])
  } else {

    #Perfom glm modeling on the object
    lm <- glm(read_diff ~ distance, data=disexp.df, family=poisson())

    #Calculate result from linear modeling
    disexp.attr <- data.frame(corrolation = summary(lm)$coefficients[,1][2],
                              p_value = summary(lm)$coefficients[,4][2])

  }

  return(disexp.attr)

}
