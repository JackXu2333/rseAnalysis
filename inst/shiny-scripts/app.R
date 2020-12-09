library(shiny)
library(rseAnalysis)

ui <- fluidPage(
  # Define UI for random distribution app ----
  ui <- fluidPage(

    # App title ----
    titlePanel("RNA Distance & Gene Expression Analysis"),

    # Sidebar layout with input and output definitions ----
    sidebarLayout(

      # Sidebar panel for inputs ----
      sidebarPanel(

        # Input: Upload files
        fileInput(inputId = "normalDot", "Select Normal RNA Structure file", accept = ".csv"),
        checkboxInput(inputId = "normalDotHeader", "Header", value = TRUE, width = NULL),

        fileInput(inputId = "mutatedDot", "Select Mutated RNA Structure file", accept = ".csv"),
        checkboxInput(inputId = "mutatedDotHeader", "Header", value = TRUE, width = NULL),

        fileInput(inputId = "geneExp", "Select gene expression file", accept = ".csv"),
        checkboxInput(inputId = "geneExpHeader", "Header", value = TRUE, width = NULL),

        # Input: Select the random distribution type ----
        radioButtons("modelMethod", "Modeling method",
                     c("Linear model" = "linear",
                       "Log Transformation" = "log")),

        # Input: Select the distance calculation method ----
        radioButtons("distaneMethod", "RNA Distance calculation method",
                     c("gscVisualizer (Zhiwen T)" = "gsc",
                       "RNA Distance (Walter F)"  = "RNADis")),

        # br() element to introduce extra vertical spacing ----
        br(),

        # Button
        actionButton("Compute", "Compute", icon("refresh")),

        # Checkbox for running with example
        checkboxInput(inputId = "sampleCompute", "Run Example", value = FALSE, width = NULL),


      ),

      # Main panel for displaying outputs ----
      mainPanel(

        verbatimTextOutput("value"),

        # Output: Tabset w/ plot, summary, and table ----
        tabsetPanel(type = "tabs",
                    tabPanel("Corrolation",
                             plotOutput("scatterPlot"), verbatimTextOutput("summary")),
                    tabPanel("Gene Expression Visualization",
                             plotOutput("GeneExpressionBoxPlot"), plotOutput("GeneExpressionDensityPlot")),
                    tabPanel("RNA Distance Visualization",
                             plotOutput("RNADistanceBoxPlot"), plotOutput("RNADistanceDensityPlot")),
                    tabPanel("Gene Expression Data", tableOutput("geneExpData")),
                    tabPanel("RNA Distance Data", tableOutput("rnaDisData"))
        )
      )
    )
  )
)

server <- function(input, output) {

  # Reactive expression to generate the requested distribution ----
  # This is called whenever the inputs change. The output functions
  # defined below then use the value computed from this expression

  recomputed <- eventReactive(input$Compute, {

    if (input$sampleCompute) {

      # Display header notation
      output$value <- renderText({"Running Example from rseAnalysis package"})

      # Read structure data
      struct.ori <- read.csv(system.file("extdata", "vignetteSampleORI.csv", package = "rseAnalysis"))
      struct.alt <- read.csv(system.file("extdata", "vignetteSampleALT.csv", package = "rseAnalysis"))

      # Read expression Data
      expression <- read.csv(system.file("extdata", "test.csv", package = "rseAnalysis"), header = TRUE)


    } else {


      # Read structure data
      struct.ori <- read.csv(input$normalDot$datapath, header = input$normalDotHeader)
      struct.alt <- read.csv(input$mutatedDot$datapath, header = input$mutatedDotHeader)

      # Read expression Data
      expression <- read.csv(input$geneExp$datapath, header = input$geneExpHeader)

    }

    # Run prediction
    RNA.distance <- predictDistance(name = struct.ori[,1]
                                    , struct.ori = struct.ori[,2]
                                    , struct.alt = struct.alt[,2]
                                    , method = input$distaneMethod)


    expression <- subset(expression, Read.Type == "reads_per_million_miRNA_mapped")[1:200, ]

    result <- Analysis.DISEXP(dis.name = struct.alt$X, dis.distance = RNA.distance,
                              exp.tumor = expression$Sample, exp.sample = expression$Normal, method = input$modelMethod, showPlot = FALSE)

    result <- list(stats = result$stats, plots = result$plots
                   , data = list( RNADistance = cbind(name = struct.ori$X, distance = RNA.distance), geneExpression = expression))

  })

  # Plot main result
  output$scatterPlot <- renderPlot ({
    if ((is.null(input$normalDot) || is.null(input$mutatedDot) || is.null(input$geneExp)) & !input$sampleCompute) return(NULL)
    recomputed()$plots$ScatterPlot })

  output$summary <- renderPrint ({
    if ((is.null(input$normalDot) || is.null(input$mutatedDot) || is.null(input$geneExp)) & !input$sampleCompute) return(NULL)
    recomputed()$stats })

  # Plot associated plots
  output$GeneExpressionBoxPlot <- renderPlot ({
    if ((is.null(input$normalDot) || is.null(input$mutatedDot) || is.null(input$geneExp)) & !input$sampleCompute) return(NULL)
    recomputed()$plots$GeneExpressionBoxPlot })

  output$GeneExpressionDensityPlot <- renderPlot ({
    if ((is.null(input$normalDot) || is.null(input$mutatedDot) || is.null(input$geneExp)) & !input$sampleCompute) return(NULL)
    recomputed()$plots$GeneExpressionDensityPlot })

  output$RNADistanceBoxPlot <- renderPlot ({
    if ((is.null(input$normalDot) || is.null(input$mutatedDot) || is.null(input$geneExp)) & !input$sampleCompute) return(NULL)
    recomputed()$plots$RNADistanceBoxPlot })

  output$RNADistanceDensityPlot <- renderPlot ({
    if ((is.null(input$normalDot) || is.null(input$mutatedDot) || is.null(input$geneExp)) & !input$sampleCompute) return(NULL)
    recomputed()$plots$RNADistanceDensityPlot })

  # Display data
  output$geneExpData <- renderTable ({
    if ((is.null(input$normalDot) || is.null(input$mutatedDot) || is.null(input$geneExp)) & !input$sampleCompute) return(NULL)
    recomputed()$data$RNADistance })

  output$rnaDisData <- renderTable ({
    if ((is.null(input$normalDot) || is.null(input$mutatedDot) || is.null(input$geneExp)) & !input$sampleCompute) return(NULL)
    recomputed()$data$geneExpression })

}

shinyApp(ui = ui, server = server)

# [END]
