summary.fitSTREBLUP <- function(fit) {
  out <- matrix(c(fit$beta, fit$sigma, fit$rho), ncol = 1)
  rownames(out) <- c(rownames(fit$beta), "sigmaSquaredSAR", "sigmaSquaredAR", "rhoSAR", "rhoAR")
  colnames(out) <- "estimated coefficient"
  as.table(out)
}

#predict.fitSTREBLUP <- function(fit) eval(expression(xb + u1 + u2), envir=fit$data)