#' Fit Fay-Herriot Model
#' 
#' @param formula for fixed effects
#' @param vardir variable name as \code{character} of known sampling variances
#' @param idName variable name as \code{character} of the domain ID
#' @param data area-level data as \code{data.frame}
#' @param optsRobust a \code{list} containing parameters for robust estimation. See \code{\link{genOptsRobust}} for available options
#' @param optsOptim a \code{list} containing parameters controlling the algorithm. See \code{\link{genOptsOptim}} for available options
#' @param type The type of model, see details
#' @param y \code{logical}, \code{TRUE} if the vector of direct estimates (response) is part of the return value. Necessary for (postestimation) computing residuals etc.
#' 
#' @details The type controls the model to be estimated. At this time only the robust Fay-Herriot model ('RFH') is supported.
#' 
#' @examples 
#' data(milk)
#' milk$SD <- milk$SD^2
#' fit <- fitfh(formula = yi ~ as.factor(MajorArea), vardir="SD", idName="SmallArea", data = milk)
#' 
fitfh <- function(formula, vardir, idName, data, optsRobust = genOptsRobust(), 
                  optsOptim = genOptsOptim(), type = "RFH", y=TRUE) {
  # Check Input:
  availableTypes <- c("RFH")
  if(!(type %in% availableTypes)) 
    stop("The type '", type, "' is not supported. Choose one in: ", 
         paste(availableTypes, collapse = ", "))
  
  # Fit Fay-Herriot model
  modelSpecs <- genModelSpecs(optsRobust, optsOptim, type)
  modelSpecs <- addModelFrame(modelSpecs, formula, vardir, idName, data)
  modelSpecs <- addStartValues(modelSpecs)
  modelSpecs <- optimizeParam(modelSpecs)
  modelSpecs <- optimizeRE(modelSpecs)
  modelSpecs <- addPrediction(modelSpecs)
  out <- list(call = match.call(), prediction = modelSpecs$prediction, fitparam = modelSpecs$fitparam, 
              fitre = modelSpecs$fitre)
  if (y) out$y <- modelSpecs$y
  class(out) <- type
  out
}