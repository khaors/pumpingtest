% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/theis_utilities.R
\name{theis_solution_initial}
\alias{theis_solution_initial}
\title{theis_solution_initial}
\usage{
theis_solution_initial(ptest)
}
\arguments{
\item{ptest}{Object of a pumping_test class}
}
\value{
A list with the values of \eqn{a} and \eqn{t0}
}
\description{
Function to estimate the initial parameters of the Theis solution
}
\examples{
# Data from a confined aquifer
data(theis)
ptest <- pumping_test('Well1', Q = 1.3e-3, r = 200, t = theis$t, s = theis$s)
res <- theis_solution_initial(ptest)
print(res)
}
\references{
Theis, C. The relation between the lowering of the piezometric surface and the rate
and duration of discharge of a well using groundwater storage. Transactions of the
American Geophysical Union, 1935, 2, 519-524
}
\seealso{
Other theis functions: \code{\link{exponential_integral}},
  \code{\link{theis_WF_LT_dlogt}},
  \code{\link{theis_WF_LT}},
  \code{\link{theis_calculate_parameters}},
  \code{\link{theis_solution_dlogt}},
  \code{\link{theis_solution}},
  \code{\link{theis_well_function_large}},
  \code{\link{theis_well_function_small}},
  \code{\link{theis_well_function}}
}
\author{
Oscar Garcia-Cabrejo \email{khaors@gmail.com}
}
\concept{theis functions}
