library(SAE)
library(microbenchmark)
library(ggplot2)


benchi <- microbenchmark(
  reSAR1(wMatrix(1000), 0.5),
  reSar1Var(w0Matrix(1000), sarCorr=0.5))

autoplot(benchi)