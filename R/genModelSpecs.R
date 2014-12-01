genModelSpecs <- function(optsRobust = genOptsRobust(), optsOptim = genOptsOptim(), type) {
  out <- new.env()
  mapply("assign", names(optsRobust), optsRobust, MoreArgs = list(envir = out))
  mapply("assign", names(optsOptim), optsOptim, MoreArgs = list(envir = out))
  out$psiFunction <- psiOne
  out$fitparam <- data.frame(param = character(), m = numeric(), stepIterations = numeric(), 
                             stepParam = list(), returnStatus = numeric(),
                             timeElapsed = numeric(), stringsAsFactors=FALSE)
  out$fitparam$stepParam <- list()  
  class(out) <- c(class(out), paste("MS", type, sep = ""))
  out
}

#' genOptsRobust
#' 
#' Constructs a list with necessary elements for \code{fitfh} - only a convinience function.
#' 
#' @param k 
#' @param K 
#' 
#' @return A named list with \code{k} and \code{K} as elements.
#' @export
genOptsRobust <- function(k = 1.345, K = 2 * pnorm(k) - 1 - 2 * k * dnorm(k) + 2 * (k^2) * (1-pnorm(k))) {
  list(k = k, K = K)
}

#' genOptsOptim
#' 
#' Constructs a list with necessary elements for \code{fitfh} - only a convinience function.
#' 
#' @param tol tolerance criterion for algorithms
#' @param maxIter maximal number of iterations, used as a parameter in each optimization step (not cummulative)
#' @param progress logical - print progress to the console
#' 
#' @return A named list with \code{tol}, \code{maxIter} and \code{progress} as elements.
#' @export
genOptsOptim <- function(tol = 1e-6, maxIter = 100, progress = FALSE) {
  list(tol = tol, maxIter = maxIter, progress = progress)
}