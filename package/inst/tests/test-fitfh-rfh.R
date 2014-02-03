context(desc="Testing robust FH:")

test_that(desc="fitfh is producing 'good'/'expected' results:", {
  library(saedevel)
  library(sae)
  tmp <- genModelSpecs(type = "test")
  expect_error(addModelFrame(tmp, y ~ x, "x", "x", data.frame(y = rnorm(100), x = rnorm(100))))
  
  modelSpecs <- genModelSpecs(type = "RFH")
  
  data(milk, envir=environment())
  milk$SD <- milk$SD^2
  
  modelSpecs <- addModelFrame(modelSpecs, formula = yi ~ MajorArea, vardir="SD", idName="SmallArea", data = milk)
  modelSpecs <- addStartValues(modelSpecs)
  
  modelSpecs <- optimizeParam(modelSpecs)
  modelSpecs <- optimizeRE(modelSpecs)
    
  fh <- eblupFH(yi ~ MajorArea, SD, method = "REML", MAXITER = 100, PRECISION = 0.0001, data = milk)$eblup
  
  rfh <- modelSpecs$X %*% modelSpecs$beta + modelSpecs$fitre$x
  expect_that(fh, equals(rfh, tolerance = 0.015))
  
  fit <- fitfh(formula = yi ~ MajorArea, vardir="SD", idName="SmallArea", data = milk, type = "RFH",
               optsRobust = genOptsRobust(k = 10000))
  
  expect_that(fh, equals(fit$prediction, tolerance = 0.015))
  expect_that(1 > confint(lm(fh ~ fit$prediction))[2, 1] & 
                1 < confint(lm(fh ~ fit$prediction))[2, 2], is_true())
})


