disable <- function(x) {
  if (inherits(x, 'shiny.tag')) {
    if (x$name %in% c('input', 'select'))
      x$attribs$disabled <- 'disabled'
    x$children <- disable(x$children)
  }
  else if (is.list(x) && length(x) > 0) {
    for (i in 1:length(x))
      x[[i]] <- disable(x[[i]])
  }
  x
}

shinyUI(
  pageWithSidebar(
  headerPanel("simViewer - 0.1"),
  
  # Input
  sidebarPanel(
    p(strong("Choose Eval Data")),
    fileInput(inputId= "file", label="",
              accept=c(".Rdata", "Rdata containig simulation results")),

    wellPanel(
      p(strong("Type")),
      selectInput(inputId = "SomeThingToSelect", label = "",
                  choices = c("1" = "passedToServer1",
                              "2" = "passedToServer2")),      
      br(),
      p(strong("Something to check")),
      checkboxInput(inputId = "something", label = "something", value = TRUE)),
    
    
#     #     wellPanel(
#     p(strong("Playing minute")),
#     checkboxInput(inputId='checkSlider', 'Slider', TRUE),
#     conditionalPanel(condition = "input.checkSlider == true",
#                     sliderInput(inputId = "minute1",
#                                 label = "",
#                                 min = 1, max = 90, step = 1, value = 45,
#                                 ticks = TRUE
# #                               ,animate = list(interval = 666)
#                                 )
#                     ),
#     conditionalPanel(condition = "input.checkSlider != true",
#                      numericInput(inputId = "minute2", label="", value=0 , min=0, max=90)
#                      ),        
    #     br(),
    checkboxInput(inputId = "logY", label = "Logarithmic scale", value = FALSE),
    checkboxInput(inputId = "grid", label = "Grid", value = FALSE),
    uiOutput("varSelect")
    
  )
  #   )
  ,
  
  
  
  # Output
  mainPanel(
    tabsetPanel(
      tabPanel("Plot", 
               conditionalPanel(condition = "input.something",
                                # br(),
                                # div(plotOutput(outputId = "slidingArrow")),
                                div(plotOutput(outputId = "somePlot"))
               )
      ),
      
      tabPanel(title= "Tab: Estimated Parameters",
               tableOutput("parameterData")),
      tabPanel(title= "Tab: Eval Data",
               verbatimTextOutput("evalSummary"))
    )
  )
  )
)