summary.fitSTREBLUP <- function(fit) {
  out <- matrix(c(fit$beta, fit$sigma, fit$rho), ncol = 1)
  rownames(out) <- c(rownames(fit$beta), "sigmaSquaredSAR", "sigmaSquaredAR", "rhoSAR", "rhoAR")
  colnames(out) <- "estimated coefficient"
  as.table(out)
}

"+.simSetup" <- function(x, y) {
  datListX <- slot(x, "data")
  datListY <- slot(y, "data")
  
  datList <- mapply(function(dat1, dat2) data.frame(dat1, dat2$yHat),
                    datListX, datListY, SIMPLIFY = FALSE)
  
  slot(x, "data") <- datList
  x
}

summary.RFH <- function(fit, ...) {
  out <- list(summaryPred = summary(as.numeric(fit$prediction)))
  out$coeffFixed <- fit$fitparam[fit$fitparam$m == max(fit$fitparam$m) & fit$fitparam$param == "beta", "stepParam"][[1]]
  colnames(out$coeffFixed) <- "Estimate"
  out$coeffRandom <- matrix(fit$fitparam[fit$fitparam$m == max(fit$fitparam$m) & fit$fitparam$param == "reVar", "stepParam"][[1]])
  colnames(out$coeffRandom) <- "Estimate"
  rownames(out$coeffRandom) <- "Variance Parameter"
  
  cat("\nCall:\n")
  print(fit$call)
  cat("\n\nCoefficients (fixed-effects):\n")
  print(out$coeffFixed)
  cat("\nCoefficients (random-effects):\n")
  print(out$coeffRandom)
  cat("\n\nSummary for predicted dependent variable:\n")
  cat("\t", out$summaryPred)
  cat("\n\n")
  invisible(out)
}

