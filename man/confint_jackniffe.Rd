% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pumpingtest.R
\name{confint_jackniffe}
\alias{confint_jackniffe}
\title{confint_jackniffe}
\usage{
confint_jackniffe(ptest, level = c(0.025, 0.975), slug = F)
}
\arguments{
\item{ptest}{A pumping_test object}

\item{level}{A numeric vector with the significance levels of the confidence intervals.}

\item{slug}{A logical flag indicating if a slug test is used in the analysis.}
}
\value{
A list with the following entries:
\itemize{
\item method: A character string specifying the method used in the CI estimation (jackniffe)
\item hydraulic.parameters.names: A character vector with the names of the hydraulic
parameters.
\item hydraulic.parameters.ci: A matrix with the requested confidence of the hydraulic
parameters.
}
}
\description{
Function to calculate the
}
\seealso{
Other base functions: \code{\link{additional.parameters<-}},
  \code{\link{confint.pumping_test}},
  \code{\link{confint_bootstrap}},
  \code{\link{confint_wald}}, \code{\link{estimated<-}},
  \code{\link{evaluate}}, \code{\link{fit.optimization}},
  \code{\link{fit.parameters<-}},
  \code{\link{fit.sampling}}, \code{\link{fit}},
  \code{\link{hydraulic.parameter.names<-}},
  \code{\link{hydraulic.parameters<-}},
  \code{\link{model.parameters}}, \code{\link{model<-}},
  \code{\link{plot.pumping_test}},
  \code{\link{plot_model_diagnostic}},
  \code{\link{plot_sample_influence}},
  \code{\link{plot_uncert}},
  \code{\link{print.pumping_test}},
  \code{\link{pumping_test}}, \code{\link{simulate}},
  \code{\link{summary.pumping_test}}
}
\author{
Oscar Garcia-Cabrejo \email{khaors@gmail.com}
}
\concept{base functions}
