% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/jacob_lohman_utilities.R
\name{jacob_lohman_G_WF_LT}
\alias{jacob_lohman_G_WF_LT}
\title{jacob_lohman_G_WF_LT}
\usage{
jacob_lohman_G_WF_LT(p, ...)
}
\arguments{
\item{p}{Laplace parameter}

\item{...}{Additional parameters}
}
\value{
Value of the dimensionless discharge in the Laplace domain of the Jacob-Lohman solution
}
\description{
Function to calculate the well function of the Jacob-Lohman solution (Discharge) in the
Laplace domain
}
\examples{
p <- 0.1
qp <- jacob_lohman_G_WF_LT(p)
print(qp)
}
\references{
Jacob, C. E. & Lohman, S. W. Nonsteady Flow to a Well Of Constant Drawdown. American
Geophysical Union, 1952, 33, 10.
}
\seealso{
Other jacob_lohman functions: \code{\link{jacob_lohman_F_calculate_parameters}},
  \code{\link{jacob_lohman_F_solution_initial}},
  \code{\link{jacob_lohman_F_solution}},
  \code{\link{jacob_lohman_F_well_function}},
  \code{\link{jacob_lohman_G_calculate_parameters}},
  \code{\link{jacob_lohman_G_solution_initial}},
  \code{\link{jacob_lohman_G_solution}},
  \code{\link{jacob_lohman_G_well_function}},
  \code{\link{jacob_lohman_WF_F_LT}}
}
\author{
Oscar Garcia-Cabrejo \email{khaors@gmail.com}
}
\concept{jacob_lohman functions}
