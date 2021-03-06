% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/flow_dimension_utilities.R
\name{flow_dimension}
\alias{flow_dimension}
\title{flow_dimension}
\usage{
flow_dimension(t, s, d = 2, method = "central")
}
\arguments{
\item{t}{Numeric vector with time values}

\item{s}{Numeric vector with drawdown values}

\item{d}{Derivative parameter. if method equals to bourdet then d is equal to
the number of adjacent values used in the derivative calculation. If method is
equal to spline then d is equal to the number of knots used in the interpolation
of the drawdown data. In this case a value of d=20 to d=30 is recommended. If
method is equal to spane then d is equal to the number of points used in the linear
regression approach.}

\item{method}{Method to calculate the derivative. See log_derivative}
}
\value{
This function returns a list with components:
\itemize{
\item x: Numeric vector with times at which the second derivative is evaluated
\item y: Numeric vector with the values of the second derivative of drawdown
\item n: Numeric vector with the flow dimension values evaluated at each time
}
}
\description{
Function to calculate the flow dimension using different approaches
}
\examples{
data(boulton)
t <- boulton$t
s <- boulton$s
boulton.fdim <- flow_dimension(t,s, method = "smoothspline")
plot(boulton.fdim$x, boulton.fdim$n, type = "p", log = "x", 
ylim = c(0, 10))
}
\author{
Oscar Garcia-Cabrejo, \email{khaors@gmail.com}
}
\concept{flow_dimension functions}
