\name{genOptsRobust}
\alias{genOptsRobust}
\title{genOptsRobust}
\usage{
  genOptsRobust(k = 1.345,
    K = 2 * pnorm(k) - 1 - 2 * k * dnorm(k) + 2 * (k^2) * (1 - pnorm(k)))
}
\arguments{
  \item{k}{}

  \item{K}{}
}
\value{
  A named list with \code{k} and \code{K} as elements.
}
\description{
  Constructs a list with necessary elements for
  \code{fitfh} - only a convinience function.
}

