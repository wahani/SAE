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
    wellPanel(p(strong("Graphic Parameter Quality Measures")),
              checkboxInput("figLine", "Zero-Line", value = FALSE),
              numericInput("figMin", "Axis-Min QM", NA),
              numericInput("figMax", "Axis-Max QM", NA)),
    wellPanel(p(strong("Select Domains of Interest")),
              selectInput("selectDomain", "", 
                          c("Single Domain" = 1, "Min-to-Value" = 2, "Value-to-Max" = 3)),
              numericInput("selectDomainValue", "Domain", NA)
              ),
    checkboxInput("debug", "Debug", value = FALSE)
  ),
  
  # Output
  mainPanel(
    tabsetPanel(
      tabPanel("Plot - Quality Measures",
               div(plotOutput(outputId = "qmPlot"))),
      tabPanel("Plot - Bland-Altman",
               div(plotOutput(outputId = "baPlot"))),
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