require(SAE)
require(ggplot2)
require(reshape2)

shinyServer(function(input, output) {
  
  ################################### Data-Reactives ###########################
  
  # Load data
  getData <- reactive({
    
    # return NULL if nothing is selected
    if(is.null(input$file)) return(list(estData = data.frame(),
                                        parData = data.frame()))
    
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
    
    #aggregate(as.matrix(estData[,grepl("est", names(estData))]) ~ domain, data = estData, FUN = mean)
    
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
  
  
  # Erstelle Reactives fuer die gewuenschten Ausgabe-Dataframes
  qmData <- reactive({
    #browser()
    
    #Daten aufrufen:
    estData <- getData()$estData
    
    # Quality Measures for each scenario and each domain:
    dataList <- split(estData, list(estData$scenario, estData$domain), drop = TRUE)
    
    qmData <- list(BIAS = function(mu, est) 
                     mean(est - mu),
                   RBIAS = function(mu, est) 
                     mean((est - mu) / mu),
                   ABIAS = function(mu, est) 
                     mean(abs(est - mu)),
                   MSE = function(mu, est)
                     mean((est - mu)^2),
                   RRMSE = function(mu, est)
                     mean(sqrt(mean(((est - mu)/mu)^2))))
    
    qmData <- lapply(qmData, 
           function(fun, dataList) {
             do.call(rbind, lapply(dataList, function(dat) {
               estimators <- names(dat)[grepl("est.", names(dat))]
               for (est in estimators) {
                 dat[, est] <- fun(dat[, "y"], dat[, est])
               }
               rownames(dat) <- NULL
               dat[1, c("domain", "scenario", estimators)]
             }))
           }, dataList = dataList)
    
    return(qmData)
    
  })
  
  ################################### input-Reactives###########################
  
  inputEstimator <- reactive(input$estimatorSelect)
  inputScenario <- reactive(input$scenarioSelect)
  inputQualityMeasure <- reactive(input$qualityMeasureSelect)
  inputFigLine <- reactive(input$figLine)
  inputFigMinMax <- reactive({
    c(input$figMin, input$figMax)
  })
  
  ################################### PlotFunctions ############################
  
  somePlotFunction <- reactive({
    browser()
    # Zugriff auf die Reactives
        
    dat <- qmData()[[inputQualityMeasure()]]
    
    somePlot <- ggplot(dat) + geom_point(aes(get("est.FH"), get("y")))
    
    return(somePlot)
  })
  
  ################################### OUTPUT ###################################
  
  # Grafiken
  output$qmPlot <- renderPlot({
    #browser()
    
    # Zugriff auf die Reactives
    dat <- qmData()[[inputQualityMeasure()]]
    scenario <- inputScenario()
    estimator <- inputEstimator()
    qm <- inputQualityMeasure()
    minMax <- inputFigMinMax()
    
    # Filter:
    dat <- dat[which(dat$scenario %in% scenario), c("domain", "scenario", estimator)]
    dat <- melt(dat, id.vars = c("domain", "scenario"), measure.vars = estimator, 
         variable.name = "Estimator", value.name = qm)
    levels(dat[, "Estimator"]) <- gsub("est.", "", estimator)
    
    # Plot:
    somePlot <- ggplot(dat, aes_string(x = "Estimator", y = qm)) + geom_boxplot() +
      coord_flip() + facet_grid(scenario ~.)
    
    if(inputFigLine()) somePlot <- somePlot + geom_hline(aes(yintercept = 0), colour = "red")
    if(!all(is.na(minMax))) 
      somePlot <- somePlot + 
      coord_flip(ylim = c(if(is.na(minMax[1])) min(dat[, qm])*1.1 else minMax[1],
                          if(is.na(minMax[2])) max(dat[, qm])*1.1 else minMax[2]))
    
    print(somePlot)
    
  })
  
  # Output summarys
  output$estimatorSummary  <- renderPrint({
    dat <- evalData()
    if(is.null(inputEstimator())) {
      summary(dat[, c("domain", "r")])
    } else {
      dat <- dat[which(dat$scenario %in% inputScenario()), inputEstimator()]
      summary(dat)
    }
  })
  
  output$qmSummary  <- renderPrint({
    dat <- qmData()[[inputQualityMeasure()]]
    
    tmp <- lapply(as.list(inputScenario()), 
                  function(scenario) {
                    cat(scenario, "\n")
                    dat <- dat[which(dat$scenario == scenario), inputEstimator()]
                    print(summary(dat))
                  })
      
#     }
  })
  
  # Output Tables
  output$parameterData  <- renderTable({
    qmData()
    parameterData()                               
  }, booktabs = TRUE, include.rownames = FALSE)
  
  
  ############################ OUTPUT - Dynamic UI #############################
  
  # Filter Data
  output$estimatorSelect <- renderUI({
    dat <- evalData()
    estimator <- names(dat)[grepl("est.", names(dat))]
    names(estimator) <- gsub("est.", "", estimator)
    checkboxGroupInput(inputId = "estimatorSelect", label = "Select Estimators:",
                       estimator, selected = names(dat)[grepl("est.", names(dat))])
  })
  
  output$scenarioSelect <- renderUI({
    dat <- evalData()
    scenario <- levels(dat$scenario)
    checkboxGroupInput(inputId = "scenarioSelect", label = "Select Scenarios:",
                       scenario, selected = scenario[1])
  })
})