# in R
#######
rm(list=ls())
#beachte: Eingabe 'by row'

#pos-def matrix:
k   <- 2000
rho <- .3
S       <- matrix(rep(rho, k*k), nrow=k)
diag(S) <- 1
dat <- mvrnorm(10000, mu=rep(0,k), Sigma=S)
R <- cor(dat)

# Inverse:
# Standard:
system.time(A <- solve(R))
# Was auch immer die Überlegung hierbei ist, es ist nicht schneller:
system.time(B <- qr.solve(R))
# Inverse für Cholesky-Decomposition scheint am schnellsten zu sein: (immer noch langsam)
system.time(C <- chol2inv(chol(R))) 



######################################
## In Schleifen:
######################################

#pos-def matrix: etwas kleiner!
k   <- 800
S       <- matrix(rep(rho, k*k), nrow=k)
diag(S) <- 1
dat <- mvrnorm(10000, mu=rep(0,k), Sigma=S)
R <- cor(dat)

invChol <- function(R, iter = 10) {
  system.time({
    Rfor <- list()
    for (i in 1:iter) {
      Rfor[[i]] <- chol2inv(chol(R), LINPACK = FALSE)
    }
  })
}

invSolve <- function(R, iter = 10) {
  system.time({
    Rfor <- list()
    for (i in 1:iter) {
      Rfor[[i]] <- solve(R)
    }
  })
}

iter = 100

# Keine überraschungen, die Berechnungszeit steigt linear:
invChol(R, iter)
invSolve(R, iter)


# Bringt hierbei anscheinen nichts.
library(compiler)
invCholCmp <- cmpfun(invChol)
invSolveCmp <- cmpfun(invSolve)

invCholCmp(R, iter)
invSolveCmp(R, iter)


# In Kombination mit matrixMultiplikation:
# berechne viel unsinn:
# (R'R)^-1 * (R'R)

k <- 3000
R <- matrix(rnorm(k^2), ncol = k)

# Alles was ich bisher an Performance-Einsparungen lernen konnte:
system.time({
  tmp1 <- crossprod(R)
  tmp2 <- chol2inv(chol(tmp1))
  check <- tmp2 %*% tmp1
})

# Dauert zu lange - gerne ausprobieren:
# system.time({
#   chol2inv(chol(t(R) %*% R)) %*% t(R) %*% R
# })
# 
# # Zum Vergleich:
# system.time({
#   solve(t(R) %*% R) %*% t(R) %*% R
# })


############################# Vergleich mit RCpp ###############################
library(SAE)

# Berechnung ohne spezielle Funktionen etc.
system.time({
  cppArmFunc(R)
  })

# Nutzt die Info, dass R'R positiv-definit ist
system.time({
  cppArmFuncOpt(R)
})

# # Macht irgendetwas fürchterlich langsames
# system.time({
#   cppArmFuncPseudo(R)
# })

# Rechnet intern die inverse mit der Cholesky-Decomposition - Die decomposition
# ist implementiert, die Inverse dafür nicht.
system.time({
  cppArmFuncChol(R)
})





