require(SAE)
require(ggplot2)

shinyServer(function(input, output) {
  
  ################################### Data-Reactives ###########################
  
  # Load data
  getData <- reactive({
    
    # return NULL if nothing is selected
    if(is.null(input$file)) return(NULL)
    
    # else:
    file <- input$file$datapath
    load(file, envir = environment()) #object in env is named: simData
    
    # Aufbereitung:
    # ...
    
    return(simData)
  })
  
  evalData <- reactive({
    
    # Daten aufrufen:
    simData <- getData()
    
    # Filtern auf Domain
    # ...
    
    # Filtern auf Gütemaße zum Anzeigen - Spalten
    # ...
    
    # Temp
    simData <- data.frame(x = 1:10, y = rnorm(10))
    
    return(simData)
      
  })
  
  
  # Erstelle Reactives fuer die gewuenschten Ausgabe-Dataframes
  parameterData <- reactive({
    
    #Daten aufrufen:
    simData <- getData()
    
    # ...
    
    # Temp
    simData <- data.frame(x = 1:10, y = rnorm(10))
    
    return(simData)
    
  })
    
  somePlotReactive <- reactive({
    x <- 0
    if(input$checkSlider==FALSE){
      x<-input$minute2 
      print(input$minute2 )
    }
    
    if(input$checkSlider==TRUE){
      x<-input$minute1  
      print(input$minute1 )
    }
    x
  })

  logY <- reactive( input$logY )
  grid <- reactive( input$grid)
  referenceLine <- reactive( input$betfair )
  
  ################################### PlotFunctions ############################
  
  somePlotFunction <- function(){
    
    # Zugriff auf die Reactives
    dat <- evalData()
    
    isLog <- if(logY()) "y" else ""
    isGrid <- grid()
    
    somePlot <- ggplot(dat) + geom_point(aes(x, y))
    
    return(somePlot)
  }
  
  ################################### OUTPUT ###################################
  
  # Grafiken
  output$somePlot <- renderPlot({ print(somePlotFunction()) })  
  
  # Output Tables
  output$evalData  <- renderTable({
    evalData()                               
  }, booktabs = TRUE, include.rownames = FALSE)
  
  # Output Tables
  output$parameterData  <- renderTable({
    parameterData()                               
  }, booktabs = TRUE, include.rownames = FALSE)
  
})