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
#' @param method Selection of linear or gaussian log link function for regression
#' @param seperate TRUE and FALSE varible if TRUE the output image will be shown separately
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
#'              dis.name = c("hsa-let-7a-1", "hsa-let-7a-1",
#'              "hsa-let-7a-3", "hsa-let-7a-3", "hsa-let-7a-3"),
#'              dis.distance = as.integer(c(10, 35, 91, 100, 92, 5)),
#'              exp.tumor = c(98691, 49201, 57540, 148702, 97721),
#'              exp.sample = c(23495, 23310, 13274, 19337, 14389))}
#'
#' @author Sijie Xu, \email{sijie.xu@mail.utoronto.ca}
#'
#' @export
#' @importFrom ggplot2 ggplot geom_point aes_string geom_boxplot geom_density geom_point geom_smooth
#' @importFrom gridExtra grid.arrange
#' @importFrom stats gaussian glm

Analysis.DISEXP <- function(dis.name, dis.distance, exp.tumor, exp.sample, method = "linear", seperate = TRUE){

  #Validate input
  if (typeof(dis.name) != "character" || typeof(dis.distance) != "integer"
      || typeof(exp.tumor) != "double" || typeof(exp.sample) != "double"){
    stop("Input must be packaged in list")
  } else if (length(dis.name) != length(dis.distance)){
    stop("RNA Distance data does not share the same size")
  } else if (length(exp.tumor) != length(exp.sample)){
    stop("Gene expression data does not share the same size")
  } else if (length(dis.name) != length(exp.tumor)){
    stop("RNA Distance and Gene expression data is misaligned")
  } else if (!method %in% c("linear", "log")){
    stop("select method from linear, log")
  }

  #Constructs the distance and expression data frame
  disexp.df <- data.frame(name = dis.name, distance = dis.distance, read_tumor = exp.tumor, read_sample = exp.sample)

  #Calculate the difference in read from sample to tumor
  disexp.df$read_diff <- abs(disexp.df$read_sample - disexp.df$read_tumor)

  #Perform correlation analysis based on the difference gene expression and RNA distance
  #Plot the diagonise plots

  if (method == "linear"){
    #Scatter plot with distance as x and difference in read being y
    scatter <- ggplot2::ggplot(disexp.df, ggplot2::aes_string(y="read_diff", x="distance")) +
      ggplot2::geom_point() +
      ggplot2::geom_smooth(formula = y ~ x, method="lm")
  } else {
    #Scatter plot with glm and log linkage
    scatter <- ggplot2::ggplot(disexp.df, ggplot2::aes_string(y="read_diff", x="distance")) +
      ggplot2::geom_point() +
      ggplot2::geom_smooth(method = "glm", formula = y~x, method.args = list(family = stats::gaussian(link = 'log')))
  }

  #Boxplot checking outliars on read_diff
  box.genexp <- ggplot2::ggplot(disexp.df, ggplot2::aes_string(x="name", y="read_diff")) +
    ggplot2::geom_boxplot(outlier.colour="red", outlier.shape=8,
                 outlier.size=4)

  #Boxplot checking outliars on distance
  box.dis <- ggplot2::ggplot(disexp.df, ggplot2::aes_string(x="name", y="distance")) +
    ggplot2::geom_boxplot(outlier.colour="red", outlier.shape=8,
                 outlier.size=4)

  #Density plot on distance by different RNA
  density.dis <- ggplot2::ggplot(disexp.df, ggplot2::aes_string(x="distance", color="name")) +
    ggplot2::geom_density()

  #Density plot on read_diff by different RNA
  density.genexp <- ggplot2::ggplot(disexp.df, ggplot2::aes_string(x="read_diff", color="name")) +
    ggplot2::geom_density()

  #Generate the plots according to selection of separation
  if(!seperate) {
    plot(gridExtra::grid.arrange(scatter, density.dis, box.genexp, box.dis, ncol=2, nrow=2))
  } else {
    plot(scatter)
    plot(density.dis)
    plot(box.genexp)
    plot(box.dis)
  }

  if (method == "linear"){
    #Perform linear modeling on the object
    lm <- lm(read_diff ~ distance, data=disexp.df)

    #Calculate result from linear modeling
    disexp.attr <- data.frame(correlation = stats::cor(disexp.df$distance, disexp.df$read_diff),
                            p_value = summary(lm)$coefficients[,4][2])
  } else {

    #Perfom glm modeling on the object
    lm <- stats::glm(read_diff ~ distance, data=disexp.df, stats::gaussian(link="log"))

    #Calculate result from linear modeling
    disexp.attr <- data.frame(correlation = summary(lm)$coefficients[,1][2],
                              p_value = summary(lm)$coefficients[,4][2])

  }

  return(disexp.attr)

}
