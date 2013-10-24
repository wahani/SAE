rm(list = ls(all = TRUE))

require(shiny)

runApp("simViewer")

#http://glimmer.rstudio.com/brechtdv/report/


dat <- data.frame(x = rnorm(100), y = rnorm(100))


ggplot(dat) + geom_point(aes(x = x, y = y)) +
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = 1.96 * sd(dat$y), ymax = Inf), 
            fill = "bisque", alpha = 0.02)
  


