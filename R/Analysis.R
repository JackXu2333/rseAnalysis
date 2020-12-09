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
#' @param method Selection of linear or gaussian log link function for regression (linear or log)
#' @param showPlot TRUE and FALSE variable if TRUE the output image will be shown on the run
#'
#' @return Returns an S3 object of class DISEXP with results.
#' list of output stats from the model
#' \itemize{
#' \item stats
#' \itemize{
#'   \item Correlation - Beta value of the regression model based on data
#'   \item PValue - P value calculated associated with the correlation
#' }
#' Show of four graph for validation
#' \item plots
#' \itemize{
#'   \item ScatterPlot plots of gene expression data on RNA distance data with fitted line
#'   \item GeneExpressionBoxPlot Boxplot checking outliars on gene expression
#'   \item RNADistanceBoxPlot Boxplot checking outliars on distance
#'   \item RNADistanceDensityPlot  Density plot on distance by different RNA
#'   \item GeneExpressionDensityPlot Density plot on gene expression by different RNA
#' }
#' }
#'
#' @examples
#' disexp <- Analysis.DISEXP(
#'              dis.name = c("hsa-let-7a-1", "hsa-let-7a-1",
#'              "hsa-let-7a-3", "hsa-let-7a-3", "hsa-let-7a-3"),
#'              dis.distance = as.integer(c(10, 35, 91, 100, 92)),
#'              exp.tumor = c(98691, 49201, 57540, 148702, 97721),
#'              exp.sample = c(23495, 23310, 13274, 19337, 14389))
#'
#' @author Sijie Xu, \email{sijie.xu@mail.utoronto.ca}
#'
#' @export
#' @importFrom ggplot2 ggplot geom_point aes_string geom_boxplot geom_density geom_point geom_smooth theme element_text ggtitle xlab ylab
#' @importFrom gridExtra grid.arrange
#' @importFrom stats gaussian glm

Analysis.DISEXP <- function(dis.name, dis.distance, exp.tumor, exp.sample, method = "linear", showPlot = FALSE){

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

  if (method == "linear"){
    #Scatter plot with distance as x and difference in read being y
    scatter <- ggplot2::ggplot(disexp.df, ggplot2::aes_string(y="read_diff", x="distance")) +
      ggplot2::geom_point() +
      ggplot2::geom_smooth(formula = y ~ x, method="lm") +
      ggplot2::ggtitle("Differential normalize gene expression vs RNA Gene Distance \n Linear Model") +
      ggplot2::xlab("RNA Gene distance") + ggplot2::ylab("Differential gene expresson (read)")
  } else {
    #Scatter plot with glm and log linkage
    scatter <- ggplot2::ggplot(disexp.df, ggplot2::aes_string(y="read_diff", x="distance")) +
      ggplot2::geom_point() +
      ggplot2::geom_smooth(method = "glm", formula = y~x, method.args = list(family = stats::gaussian(link = 'log'))) +
      ggplot2::ggtitle("Differential normalize gene expression vs RNA Gene Distance \n Log Link Applied") +
      ggplot2::xlab("RNA Gene distance") + ggplot2::ylab("Differential gene expresson (read)")
  }

  #Boxplot checking outliars on read_diff
  box.genexp <- ggplot2::ggplot(disexp.df, ggplot2::aes_string(x="name", y="read_diff")) +
    ggplot2::geom_boxplot(outlier.colour="red", outlier.shape=8,
                 outlier.size=4) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust=1)) +
    ggplot2::ggtitle("Boxplot of differential normalize gene expression read \n by RNA") +
    ggplot2::xlab("RNA") + ggplot2::ylab("Differential gene expresson (read)")

  #Boxplot checking outliars on distance
  box.dis <- ggplot2::ggplot(disexp.df, ggplot2::aes_string(x="name", y="distance")) +
    ggplot2::geom_boxplot(outlier.colour="red", outlier.shape=8,
                 outlier.size=4) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust=1)) +
    ggplot2::ggtitle("Boxplot of RNA gene distance read \n by RNA") +
    ggplot2::xlab("RNA") + ggplot2::ylab("RNA distance")

  #Density plot on distance by different RNA
  density.dis <- ggplot2::ggplot(disexp.df, ggplot2::aes_string(x="distance", color="name")) +
    ggplot2::geom_density() +
    ggplot2::ggtitle("Density plot of RNA gene distance read \n by RNA") +
    ggplot2::xlab("RNA distance") + ggplot2::ylab("Density")

  #Density plot on read_diff by different RNA
  density.genexp <- ggplot2::ggplot(disexp.df, ggplot2::aes_string(x="read_diff", color="name")) +
    ggplot2::geom_density() +
    ggplot2::ggtitle("Density plot of differential normalize gene expression read \n by RNA") +
    ggplot2::xlab("RNA distance") + ggplot2::ylab("Density")

  #Generate the plots according to selection of separation
  if(showPlot) {
    plot(gridExtra::grid.arrange(scatter, density.dis, box.genexp, box.dis, ncol=2, nrow=2))
  }

  if (method == "linear"){

    #Perform linear modeling on the object
    lm <- lm(read_diff ~ distance, data=disexp.df)

    #Calculate result from linear modeling
    correlation <- stats::cor(disexp.df$distance, disexp.df$read_diff)
    p_value <- summary(lm)$coefficients[,4][2][[1]]

  } else {

    #Perfom glm modeling on the object
    lm <- stats::glm(read_diff ~ distance, data=disexp.df, stats::gaussian(link="log"))

    #Calculate result from linear modeling
    correlation <- summary(lm)$coefficients[,1][2][[1]]
    p_value <- summary(lm)$coefficients[,4][2][[1]]

  }

  #Create DISEXP object for export
  AnalysisResults <- list(
    stats = list(
      Correlation = correlation,
      PValue = p_value
    ),
    plots = list(
      ScatterPlot = scatter,
      GeneExpressionBoxPlot = box.genexp,
      RNADistanceBoxPlot = box.dis,
      RNADistanceDensityPlot = density.dis,
      GeneExpressionDensityPlot= density.genexp)
    )

  class(AnalysisResults) <- "DISEXP"

  return(AnalysisResults)

}

#' Launch Shiny App For Package rseAnalysis
#'
#' A function that launches the shiny app for this package.
#' Analysis.DISEXP is a complete analysis based on user selection of
#' linear or log regression. The gene expression is calculated as the
#' absolute differences between sampled and normal gene expression data.
#' Analysis also export sets of graphs to facilitate in model selection
#' and analysis result validation. Inputs:
#'
#' FILE NormalRNAStructureFile .csv file with two column, RNA name and orginal structure in dot/bracket
#' FILE MutatedRNAStructureFile .csv file with two column, RNA name and mutated structure in dot/bracket
#' FILE geneExpressonFile .csv file with at least three column, $Read.Type indicates whether
#' it is "reads_per_million_miRNA_mapped" or "read_count", $Sample and $Normal stores the differential gene
#' expression information.
#' SELCETOR ModelMethod Selection of linear or gaussian log link function for regression (linear or log)
#' SELECTOR PredictionMethod Selection of gscVistualizer or RNADistance to predict secondary structure of RNA
#' CHECKBOX Run example to run shiny app under example files
#'
#' @return No return value but open up a shiny page.
#'
#' @examples
#' \dontrun{
#' runRSEAnalysis()
#' }
#'
#' @author Sijie Xu, \email{sijie.xu@mail.utoronto.ca}
#'
#' @export
#' @importFrom shiny runApp

runRSEAnalysis <- function() {
  appDir <- system.file("shiny-scripts",
                        package = "rseAnalysis")
  shiny::runApp(appDir, display.mode = "normal")
  return()
}

# [END]
