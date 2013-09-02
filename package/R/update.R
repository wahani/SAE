#' updateOmega1
#' 
#' @description Marhuenda et. al. (2013): page 310
#' 
updateOmega1 <- function(sarCorr, w0) {
  reSar1Var(w0, sarCorr=sarCorr)
}

#' updateOmega2
#' 
#' @description Marhuenda et. al. (2013): page 310
#'
updateOmega2 <- function(arCorr, nTime) reAr1Var(arCorr, nTime)

#' updateP
#' 
#' @description Marhuenda et. al. (2013): page 311
#'
updateP <- function(solvedV, X) {
  tmp1 <- solve(crossprod(X, solvedV) %*% X)
  solvedV - solvedV %*% X %*% tmp1 %*% crossprod(X, solvedV)
}

#' updateP
#' 
#' @description Marhuenda et. al. (2013): page 311 Remark 1
#'
updateV <- function(sigma1, Ome1, A, Z1) sigma1 * Z1 %*% Ome1 %*% t(Z1) + A


#' updateP
#' 
#' @description Marhuenda et. al. (2013): page 311 Remark 1
#'
updateA <- function(sigma2, Ome2, nDomains, nTime, sigmaSamplingError) {
  
  diagTemp <- sigma2 * omega2Diag(Ome2, nDomains)
  #samplingErrors <- seSigmaClosure(nDomains, nTime)()
  diag(diagTemp) <- diag(diagTemp) + sigmaSamplingError
  diagTemp
}

#' updateSolvedV
#' 
#' @description Marhuenda et. al. (2013): page 311 Remark 1
#'
updateSolvedV <- function(sarCorr, sigma1, arCorr, sigma2, A, Ome1, Z1) {
  solvedA <- solve(A)
  solvedA - solvedA %*% Z1 %*% solve(sigma1^(-1) * solve(Ome1) + crossprod(Z1, solvedA) %*% Z1) %*% crossprod(Z1, solvedA)
}

#' updateDerivativeClosure
#' 
#' @description Marhuenda et. al. (2013): page 312 Remark 1
#'
updateDerivativesClosure <- function(nDomains, nTime, Z1, W) {
  force(nDomains); force(nTime); force(Z1); force(W)
  
  #Some Functions - Make only sense in this environment!!!
  derOmega2ArCorr <- function (arCorr, Ome2) {
    #Derivative for Omega2 with respect to arcorr
    elements <- lapply(as.list((nTime - 1):1), function(num) arCorr^(1:num)/arCorr * (1:num))
    omega2 <- diag(1, nrow = nTime)
    omega2 <- as.vector(omega2)
    ind <- which(omega2 == 1)
    for (i in 1:(nTime - 1)) {
      omega2[(ind[i] + 1):(ind[i] + length(elements[[i]]))] <- elements[[i]]
    }
    omega2 <- matrix(omega2, nTime)
    omega2 <- omega2 + t(omega2)
    diag(omega2) <- 0
    (1 - arCorr^2)^(-1) * omega2 + 2 * arCorr * Ome2 * (1 - arCorr^2)^(-1)
  }
  
  derVSigma1 <- function(Ome1) Z1 %*% Ome1 %*% t(Z1)
  derVSarCorr <- function(sigma1, Ome1, sarCorr) (-1) * sigma1 * Z1 %*% Ome1 %*% (-W-t(W) + 2 * sarCorr * t(W) %*% W) %*% Ome1 %*% t(Z1)
  derVSigma2 <- function(Ome2, nDomains) omega2Diag(Ome2, nDomains)
  derVArCorr <- function(sigma2, arCorr, Ome2, nDomains) sigma2 * omega2Diag(derOmega2ArCorr(arCorr, Ome2), nDomains)
  
  #####
  
  function(sarCorr, sigma1, arCorr, sigma2, Ome1, Ome2, parSet = c("sigma", "rho")) {
    #Return a list with the updated derivatives of V
    if (parSet == "sigma") return(
      list(derVSigma1 = derVSigma1(Ome1),         
           derVSigma2 = derVSigma2(Ome2, nDomains))) else return(         
             list(derVSarCorr = derVSarCorr(sigma1, Ome1, sarCorr),
                  derVArCorr = derVArCorr(sigma2, arCorr, Ome2, nDomains)))
  }
}

#' updateSqrtU
#' 
#' @description Sinha & Rao (2009): page 386
#'
updateSqrtU <- function(V) {
  U <- diag(diag(V), dim(V)[1])
  # SquareRoot of a matrix with single value decomposition
  #   uDecomp <- svd(U)
  #   t(uDecomp$v %*% (t(uDecomp$u)*sqrt(uDecomp$d)))
  sqrt(U)
}

rFun <- MASS::psi.huber