test_that("omega2Diag", {
  Ome2 <- matrix(1:9, ncol = 3)
  nDomains <- 2
  expOut <- matrix(0, ncol = ncol(Ome2) * nDomains, nrow = ncol(Ome2) * nDomains)
  expOut[1:3, 1:3] <- Ome2
  expOut[4:6, 4:6] <- Ome2
  expect_that(omega2Diag(Ome2, nDomains), equals(expOut))
})