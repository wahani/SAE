\name{fitfh}
\alias{fitfh}
\title{Fit Fay-Herriot Model}
\usage{
  fitfh(formula, vardir, idName, data,
    optsRobust = genOptsRobust(),
    optsOptim = genOptsOptim(), type = "RFH", y = TRUE)
}
\arguments{
  \item{formula}{for fixed effects}

  \item{vardir}{variable name as \code{character} of known
  sampling variances}

  \item{idName}{variable name as \code{character} of the
  domain ID}

  \item{data}{area-level data as \code{data.frame}}

  \item{optsRobust}{a \code{list} containing parameters for
  robust estimation. See \code{\link{genOptsRobust}} for
  available options}

  \item{optsOptim}{a \code{list} containing parameters
  controlling the algorithm. See \code{\link{genOptsOptim}}
  for available options}

  \item{type}{The type of model, see details}

  \item{y}{\code{logical}, \code{TRUE} if the vector of
  direct estimates (response) is part of the return value.
  Necessary for (postestimation) computing residuals etc.}
}
\description{
  Fit Fay-Herriot Model
}
\details{
  The type controls the model to be estimated. At this time
  only the robust Fay-Herriot model ('RFH') is supported.
}
\examples{
data(milk)
milk$SD <- milk$SD^2
fit <- fitfh(formula = yi ~ as.factor(MajorArea), vardir="SD", idName="SmallArea", data = milk)
}

