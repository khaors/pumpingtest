% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/agarwal_skin_utilities.R
\name{agarwal_skin_solution_dlogt}
\alias{agarwal_skin_solution_dlogt}
\title{agarwal_skin_solution_dlogt}
\usage{
agarwal_skin_solution_dlogt(ptest, a, t0, sigma, t)
}
\arguments{
\item{ptest}{A pumping_test object}

\item{a}{Slope of the straight line fitted to the drawdown data using the
Cooper-Jacob approach}

\item{t0}{Intercept of the straight line fitted to the drawdown data using the
Cooper-Jacob approach}

\item{sigma}{Value of the skin factor}

\item{t}{Numeric vector with time}
}
\value{
A list with
\itemize{
\item x: time coordinates
\item y: Derivative value
}
}
\description{
Function to calculate the derivative of drawdown with respect to log of time of the Agarwal solution (skin and wellbore storage)
}
\examples{
data(papadopulos_cooper)
}
\references{
Agarwal, R.; Al-Hussainy, R. & Ramey, H. An Investigation of Wellbore Storage and
Skin Effect in Unsteady Liquid Flow: I. Analytical Treatment. Society of Petroleum
Engineers Journal., 1970
}
\seealso{
Other agarwal_skin functions: \code{\link{agarwal_skin_WF_LT}},
  \code{\link{agarwal_skin_calculate_parameters}},
  \code{\link{agarwal_skin_solution_initial}},
  \code{\link{agarwal_skin_solution}},
  \code{\link{agarwal_skin_well_function}}
}
\author{
Oscar Garcia-Cabrejo \email{khaors@gmail.com}
}
\concept{agarwal_skin functions}
