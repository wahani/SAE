combineSimResults <- function(simResults, functionNames) {
  simSetup <- Reduce("+", x=simResults)
  datList <- slot(simSetup, "data")
  datList <- lapply(datList, 
                    function(dat, functionNames) {
                      names(dat)[grepl("yHat", names(dat))] <- paste("yHat", functionNames, sep = ".")
                      dat
                    }
                    ,
                    functionNames = functionNames)
  
  slot(simSetup, "data") <- datList
  simSetup
}

consoleOutput <- function(consoleOutput) if (consoleOutput) cat(".") else NULL

omega2Diag <- function(Ome2, nDomains) {
  out <- matrix(0, ncol = ncol(Ome2) * nDomains, nrow = ncol(Ome2) * nDomains)
  for (i in 1:nDomains) {
    #browser()
    ind <- 1:ncol(Ome2) + (i-1) * ncol(Ome2)
    
    out[ind, ind] <- Ome2
  }
  out
}


#' makeXY
#' 
#' @description extract respone vector and design matrix from data
#' 
#' @param formula formula object
#' @param data data.frame
#' 
makeXY <- function(formula, data){
  mf <- model.frame(formula=formula, data=data)
  x <- model.matrix(attr(mf, "terms"), data=mf)
  y <- model.response(mf)
  
  list(y = y,
       x = x)
}

getNDomains <- function(dat) length(unique(dat$Domain))
getNTime <- function(dat) length(unique(dat$Time))

#' psiOne
#' 
#' @description Psi-function
#' 
#' @param u, k, deriv see \code{\link{MASS::psi.huber}}
#' 
psiOne <- function(u,k = 1.345,deriv = 0){
  var.weights = rep(1, length(u))
  sm<-median(abs(u/sqrt(var.weights)))/0.6745
  w <- MASS::psi.huber(u/(sm * sqrt(var.weights)),k, deriv)
  if (!deriv) return(w*u) else return(w)
}