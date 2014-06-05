calcABIAS <- function(trueValues, estimates) {
  estimates[estimates == 0] <- NA
  mean(abs(as.numeric(trueValues-estimates)), na.rm=T)
}

calcRBIAS <- function(trueValues, estimates) {
  estimates[estimates == 0] <- NA
  trueValues[is.na(estimates)] <- NA
  bias <- trueValues-estimates
  median(bias / trueValues, na.rm=T)
}

calcBIAS <- function(trueValues, estimates) {
  estimates[estimates == 0] <- NA
  bias <- trueValues-estimates
  mean(bias, na.rm=T)
}

calcRRMSE <- function(trueValues, estimates) {
  estimates[estimates == 0] <- NA
  sqrt(mean(as.numeric((trueValues-estimates)/trueValues)^2, na.rm = T))
}

calcMSE <- function(trueValues, estimates) {
  estimates[estimates == 0] <- NA
  mean(as.numeric(trueValues-estimates)^2, na.rm = T)
}