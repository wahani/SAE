#' evaluate simulation results
#' 
#' @param simResults a simSetup object with simulation results
#' @param critFunctionName a function name. Function with two arguments. First is
#' the vector of true values, second the vector of predicitons.
#' @export
getEvalCrit <- function(simResults, critFunctionName = "calcRRMSE", scenario = "") {
  require(reshape2)
  critFunction <- match.fun(critFunctionName)
  data <- subset(do.call("rbind", simResults@data), Time == simResults@nTime)
  vars <- names(data)[grepl("fit", names(data))]
  newNames <- vars
  newNames[vars == "yHat.fitEB"] <- "FH"
  newNames[vars == "yHat.fitSTREBLUP"] <- "STRFH"
  newNames[vars == "yHat.fitSTEBLUP"] <- "STFH"
  
  dataList <- split(data, as.factor(data$Domain))
    
  evalList <- lapply(dataList, 
                     function(dat) {
                       evalData <- data.frame(Direct = critFunction(dat$trueY, dat$y))
                       evalData$nDirect <- simResults@n
                       #browser()
                       for (i in seq_along(vars)) {
                         evalData[newNames[i]] <- critFunction(as.numeric(dat$trueY), as.numeric(dat[[vars[i]]]))
                         evalData[paste("n", newNames[i], sep = "")] <- sum(as.numeric(dat[[vars[i]]]) != 0)
                       }
                       evalData
                     }
  )
  
  result <- do.call("rbind", evalList)
  
  result1 <- melt(result[, !grepl("^n", names(result))], variable.name="model", value.name = sub("calc", "", critFunctionName))
  result1$Domain <- rep(1:simResults@nDomains, length(vars)+1)
  #browser()
  result1$n <- melt(result[, grepl("^n", names(result))], variable.name="model", value.name = "n")$n
  result1$Scenario <- scenario
  result1
}

calcAABIAS <- function(trueValues, estimates) {
  estimates[estimates == 0] <- NA
  mean(abs(as.numeric(trueValues-estimates)), na.rm=T)
}

calcRBIAS <- function(trueValues, estimates) {
  estimates[estimates == 0] <- NA
  mean(as.numeric((trueValues-estimates)/trueValues), na.rm=T)
}

calcRRMSE <- function(trueValues, estimates) {
  estimates[estimates == 0] <- NA
  sqrt(mean(as.numeric((trueValues-estimates)/trueValues)^2, na.rm = T))
}

calcMSE <- function(trueValues, estimates) {
  estimates[estimates == 0] <- NA
  mean(as.numeric(trueValues-estimates)^2, na.rm = T)
}