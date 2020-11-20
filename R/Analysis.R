#' Perform correlation analysis based on RNA distance and Gene expression
#'
#' Analysis.DISEXP is a complete analysis based on user selection of
#' linear or log regression. The gene expression is calculated as the
#' absolute differences between sampled and normal gene expression data.
#' Analysis also export sets of graphs to facilitate in model selection
#' and analysis result validation.
#'
#' @param dis.name Set of name of RNA distance
#' @param dis.distance Set of RNA distance between mutate and original data
#' @param exp.tumor Set of reads from gene expression from tumor samples
#' @param exp.sample Set of reads from gene expression from normal samples (usually blood sample)
#' @param method Selection of linear or possion (log link) function for regression
#'
#' @return Show of four graph for validation, includes
#' \itemize{
#'   \item Scatter plot of gene expression data on RNA distance data
#'   \item Boxplot of gene expression data by each RNA
#'   \item Boxplot of RNA distance data by each RNA
#'   \item Density plot of RNA distance data, labelled by each RNA
#' }
#' list of output indicating the corrolation\itemize{
#'   \item correlation - Beta value of the regression model based on data
#'   \item p_value - P value calculated associated with the correlation
#' }
#'
#' @examples \dontrun{disexp <- Analysis.DISEXP(
#'              dis.name = c("hsa-let-7a-1", "hsa-let-7a-1", "hsa-let-7a-3", "hsa-let-7a-3", "hsa-let-7a-3"),
#'              dis.distance = c(10, 35, 91, 100, 92, 5),
#'              exp.tumor = c(98691, 49201, 57540, 148702, 97721),
#'              exp.sample = c(23495, 23310, 13274, 19337, 14389))}
#'
#' @author Sijie Xu, \email{sijie.xu@mail.utoronto.ca}
#'
#' @import dplyr
#' @import ggplot2
Analysis.DISEXP <- function(dis.name, dis.distance, exp.tumor, exp.sample, method = "linear"){

  #Validate input
  if (typeof(dis.name) != "list" || typeof(dis.distance) != "list"
      || typeof(exp.tumor) != "list" || typeof(exp.sample) != "list"){
    stop("Input must be packaged in list")
  } else if (length(dis.name) != length(dis.distance)){
    stop("RNA Distance data does not share the same size")
  } else if (length(exp.tumor) != length(exp.sample)){
    stop("Gene expression data does not share the same size")
  } else if (!method %in% c("linear", "poisson")){
    stop("select method from linear, possion")
  }

  #Constructs the distance and expression data frame
  dis <- data.frame(name = dis.name, distance = dis.distance, read_tumor = exp.tumor, read_sample = exp.sample)

  disexp.df <- merge(dis, exp, by="name")

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

  if (method == "linear"){
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
