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