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
    estData <- getData()$estData
    
    # Berechnung der Gütemaße - BIAS, RBIAS, etc.
    
    aggregate(as.matrix(estData[,grepl("est", names(estData))]) ~ domain, data = estData, FUN = mean)
    
    # Filtern auf Gütemaße zum Anzeigen - Spalten
    # ...
    
    # Temp
        
    return(estData)
      
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
    
  logY <- reactive( input$logY )
  grid <- reactive( input$grid)
  referenceLine <- reactive( input$betfair )
  
  ################################### PlotFunctions ############################
  
  somePlotFunction <- function(){
    
    # Zugriff auf die Reactives
    dat <- parameterData()
    
    isLog <- if(logY()) "y" else ""
    isGrid <- grid()
    
    somePlot <- ggplot(dat) + geom_point(aes(x, y))
    
    return(somePlot)
  }
  
  ################################### OUTPUT ###################################
  
  # Grafiken
  output$somePlot <- renderPlot({ print(somePlotFunction()) })  
  
  # Output Tables
  output$evalSummary  <- renderPrint({
    dat <- evalData()
    dat <- dat[grepl("est.", names(dat))]
    summary(dat)
  })
  
  # Output Tables
  output$parameterData  <- renderTable({
    parameterData()                               
  }, booktabs = TRUE, include.rownames = FALSE)
  
  output$varSelect <- renderUI({
    dat <- evalData()
    checkboxGroupInput(inputId = "varSelect", label = "Select variables",
                names(dat))
  })
})