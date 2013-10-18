shinyUI(
  pageWithSidebar(
  headerPanel("simViewer - 0.1"),
  
  # Input
  sidebarPanel(
    p(strong("Select input data")),
    fileInput(inputId= "file", label="",
              accept=c(".Rdata", "Rdata containig simulation results")),
    selectInput("qualityMeasureSelect", "Select a Quality Measure", 
                choices = c("BIAS", "RBIAS", "ABIAS", "MSE", "RRMSE")),
    wellPanel(p(strong("Filter Data")),
              uiOutput("estimatorSelect"),
              uiOutput("scenarioSelect")),
    wellPanel(p(strong("Graphic Parameter")),
              checkboxInput("figLine", "Zero-Line", value = FALSE),
              numericInput("figMin", "Axis-Min QM", NA),
              numericInput("figMax", "Axis-Max QM", NA))
  ),
  
  # Output
  mainPanel(
    tabsetPanel(
      tabPanel("Plot",
               img(plotOutput(outputId = "qmPlot"))
               ),
      tabPanel(title= "Tab: Estimated Parameters",
               tableOutput("parameterData")),
      tabPanel(title= "Tab: Data and Quality Measures",
               h3('Data Summary'),
               verbatimTextOutput("estimatorSummary"),
               h3('Quality Measure Summary'),
               verbatimTextOutput("qmSummary"))
    )
  )
  )
)