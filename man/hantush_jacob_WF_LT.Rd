% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hantush_jacob_utilities.R
\name{hantush_jacob_WF_LT}
\alias{hantush_jacob_WF_LT}
\title{hantush_jacob_WF_LT}
\usage{
hantush_jacob_WF_LT(p, arg1, ...)
}
\arguments{
\item{p}{Laplace parameter}

\item{arg1}{Beta parameter. The beta parameter is defined as
\deqn{\beta = \frac{r}{4} \sqrt{\frac{k_{h}^{'}S_{s}^{'}}{TS}}}
where \eqn{T} and \eqn{S} are the transmissivity and storage coefficient of the
aquifer, and \eqn{k_{h}^{'}} and \eqn{S_{s}^{'}} are the vertical hydraulic
conductivity and the Storage coefficient of the aquitard respectively.}

\item{...}{Additional parameters}
}
\value{
A numeric vector with the calculated drawdown in the Laplace domain
}
\description{
Calculates the Hantush-Jacob solution (drawdown) in the Laplace domain
}
\examples{
p <- 0.5
beta <- 0.5
res <- hantush_jacob_WF_LT(p, beta)
print(res)
}
\references{
Hantush, M. S. Flow to wells in aquifers separated by a semipervious layer. Journal
of Geophysical Research, 1967, 72, 1709.
}
\seealso{
Other hantush_jacob functions: \code{\link{hantush_jacob_WF_LT_dlogt}},
  \code{\link{hantush_jacob_calculate_parameters}},
  \code{\link{hantush_jacob_solution_dlogt}},
  \code{\link{hantush_jacob_solution_initial}},
  \code{\link{hantush_jacob_solution}},
  \code{\link{hantush_jacob_well_function}}
}
\author{
Oscar Garcia-Cabrejo \email{khaors@gmail.com}
}
\concept{hantush_jacob functions}
