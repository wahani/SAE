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
    
  fh <- eblupFH(yi ~ MajorArea, SD, method = "REML", MAXITER = 100, PRECISION = 1e-06, data = milk)$eblup
  
  rfh <- modelSpecs$X %*% modelSpecs$beta + modelSpecs$fitre$x
  expect_that(fh, equals(rfh, tolerance = 0.015))
  
  fit <- fitfh(formula = yi ~ MajorArea, vardir="SD", idName="SmallArea", data = milk, type = "RFH",
               optsRobust = genOptsRobust(k = 10000))
  
  
  expect_that(fh, equals(fit$prediction, tolerance = 0.015))
  expect_that(1 > confint(lm(fh ~ fit$prediction))[2, 1] & 
                1 < confint(lm(fh ~ fit$prediction))[2, 2], is_true())
})

test_that("RFH is handling 0 variances", {
  dat <- structure(list(idD = 1:20, 
                        y = c(1.91334788421053, 1.31219963636364, 0.465175530756437, 
                              3.02856546, 0.676592046, 0.3379209, 0.120637807017544, 
                              0.470554980362538, 0.508989952, 0.452383672961443, 
                              0.543936295833333, 0.758410128920415, 0.505185970149254, 
                              0.41766147754491, 0.324333567447743, 0.0957985490196078, 
                              0.243834842931937, 0.122924552023121, 0.151441988970877, 
                              0.135980818181818), 
                        x = c(9, 4.81818181818182, 3.32769830949285, 11.34, 
                              4.15533980582524, 3.2, 1.57894736842105, 2.88821752265861, 
                              3.496, 3.13930348258706, 2.54166666666667, 
                              4.05536332179931, 3.29850746268657, 2.70209580838323, 
                              2.67102137767221, 1.62352941176471, 1.37696335078534, 
                              1.63583815028902, 1.63561643835616, 1.70087976539589), 
                        dirVar = c(0.0691726074203466, 0.00377145639910195, 
                                   0.0216503883117238, 0.278008289005939, 
                                   0.0166537591225976, 1.37298169483161, 
                                   0.698959846705716, 2.92265293884939, 
                                   0.956783706938724, 0.422775156250175, 
                                   5.45007382273668, 0.244958830050519, 
                                   0.177044995542312, 0.229011362363644, 
                                   0.0626689312819137, 0.220433595440947, 
                                   0.251129347662288, 0.357922374351885, 
                                   0.582439697566219, 0.298090953862522)), 
                   .Names = c("idD", "y", "x", "dirVar"), class = "data.frame", 
                   row.names = c(NA, -20L))
  
  expect_true(all(fitfh(y ~ x, vardir = "dirVar", data=dat, idName="idD")$fitre$x == 0))
  
})
