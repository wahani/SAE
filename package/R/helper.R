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