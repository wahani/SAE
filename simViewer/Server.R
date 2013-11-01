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
    return(estData)
  })
  
  
  # Erstelle Reactives fuer die gewuenschten Ausgabe-Dataframes
  parameterData <- reactive({
    browser()
    #Daten aufrufen:
    simData <- getData()
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
  
  ############################### input-Reactives ##############################
   
  inputEstimator <- reactive(input$estimatorSelect)
  inputScenario <- reactive(input$scenarioSelect)
  inputQualityMeasure <- reactive(input$qualityMeasureSelect)
  inputFigLine <- reactive(input$figLine)
  inputFigMinMax <- reactive({c(input$figMin, input$figMax)})
  nDomains <- reactive({max(as.numeric(levels(evalData()$domain)))})
  inputSelectDomainRange <- reactive({
    if(input$selectDomain == 1) 
      return(input$selectDomainValue) else if(input$selectDomain == 2) 
        return(1:input$selectDomainValue) else 
          return(input$selectDomainValue:nDomains())
  })
  
  debugger <- reactive({
    if(input$debug) browser() else NULL
  })
  
  ################################### OUTPUT ###################################
  
  ############################### OUTPUT - Plots ###############################
  
  output$qmPlot <- renderPlot({
    # Boxplots for Quality Measures
    
    # Zugriff auf die Reactives
    dat <- qmData()[[inputQualityMeasure()]]
    scenario <- inputScenario()
    estimator <- inputEstimator()
    qm <- inputQualityMeasure()
    minMax <- inputFigMinMax()
    selectDomains <- inputSelectDomainRange()
    
    # Filter:
    dat <- dat[which(dat$scenario %in% scenario), c("domain", "scenario", estimator)]
    dat <- melt(dat, id.vars = c("domain", "scenario"), measure.vars = estimator, 
                variable.name = "Estimator", value.name = qm)
    levels(dat[, "Estimator"]) <- gsub("est.", "", estimator)
    
    # Plot:
    somePlot <- ggplot(dat, aes_string(x = "Estimator", y = qm)) + 
      geom_boxplot(outlier.size = 0) +
      coord_flip() + facet_grid(scenario ~.) 

    if(!all(is.na(selectDomains)))
      somePlot <- somePlot + 
      geom_point(data = subset(dat, domain %in% as.character(selectDomains)), 
                 aes_string(x = "Estimator", y = qm), 
                 colour = "dodgerblue4")
    
    if(inputFigLine()) somePlot <- somePlot + geom_hline(aes(yintercept = 0), colour = "red")
    if(!all(is.na(minMax))) 
      somePlot <- somePlot + 
      coord_flip(ylim = c(if(is.na(minMax[1])) min(dat[, qm])*1.1 else minMax[1],
                          if(is.na(minMax[2])) max(dat[, qm])*1.1 else minMax[2]))
    

    print(somePlot)
    
  })
  
  output$baPlot <- renderPlot({
    # Bland-Altman-Plot
    
    # Used reactives
    dat <- evalData()
    scenario <- inputScenario()
    estimator <- inputEstimator()
    
    # Which estimators are compared - always 2!
    if(length(estimator) == 0) 
      stop("Select one or more estimators!") 
    if (length(estimator) == 1) 
      estimator <- c("y", estimator)
    if (length(estimator) > 2) 
      estimator <- estimator[1:2]
    
    # Filter:
    dat <- dat[which(dat$scenario %in% scenario), c("domain", "scenario", estimator)]
    
    # Labels:
    xLabel <- gsub("est.", "", paste("Average of", estimator[1], "+", estimator[2]))
    yLabel <- gsub("est.", "", paste("Difference of", estimator[1], "-", estimator[2]))
    
    # Bland-Altman Statistics
    dat$Avrg <- (dat[, estimator[1]] + dat[, estimator[2]])/2
    dat$Diff <- dat[, estimator[1]] - dat[, estimator[2]]
    
    # Plot:
    somePlot <- ggplot(dat) + geom_point(aes(x = Avrg, y = Diff), alpha = 0.1) +
      labs(y = yLabel, x = xLabel) + 
      geom_hline(aes(yintercept = mean(Diff)), colour = "dodgerblue4") +
      geom_hline(aes(yintercept = 1.96 * sd(Diff)), colour = "red", linetype = 2) +
      geom_hline(aes(yintercept = -1.96 * sd(Diff)), colour = "red", linetype = 2) +
      geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = 1.96 * sd(Diff), ymax = Inf), 
                alpha = 0.01, fill = "bisque") +
      geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = -1.96 * sd(Diff)), 
                alpha = 0.01, fill = "bisque")
    
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
    debugger()
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
                       estimator, selected = names(estimator))
  })
  
  output$scenarioSelect <- renderUI({
    dat <- evalData()
    scenario <- levels(dat$scenario)
    checkboxGroupInput(inputId = "scenarioSelect", label = "Select Scenarios:",
                       scenario, selected = scenario[1])
  })
})