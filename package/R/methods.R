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

#predict.fitSTREBLUP <- function(fit) eval(expression(xb + u1 + u2), envir=fit$data)
