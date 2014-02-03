\name{genOptsOptim}
\alias{genOptsOptim}
\title{genOptsOptim}
\usage{
  genOptsOptim(tol = 1e-06, maxIter = 100,
    progress = FALSE)
}
\arguments{
  \item{tol}{tolerance criterion for algorithms}

  \item{maxIter}{maximal number of iterations, used as a
  parameter in each optimization step (not cummulative)}

  \item{progress}{logical - print progress to the console}
}
\value{
  A named list with \code{tol}, \code{maxIter} and
  \code{progress} as elements.
}
\description{
  Constructs a list with necessary elements for
  \code{fitfh} - only a convinience function.
}

