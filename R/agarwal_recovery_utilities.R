#' @section Agarwal functions (Recovery tests):
#'
#' The functions included in this section are:
#'
#' agarwal_recovery_well_function, agarwal_recovery_solution_initial, agarwal_recovery_calculate_parameters, agarwal_recovery_solution, agarwal_recovery_WF_LT, agarwal_recovery_solution_dlogt, equivalent_time
#'
#' @docType package
#' @name pumpingtest
NULL
#' @title
#' agarwal_recovery_well_function
#' @description
#' Wrapper to function theis_well_function (The recovery test is analyzed using the Theis solution)
#' @param td Numeric vector with dimensionless time
#' @param par Additional parameters
#' @return
#' A numeric vector with the dimensionless drawdown
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @export
#' @family agarwal_recovery functions
#' @references
#' Agarwal, R. (1980). A New Method To Account For Producing Time Effects When Drawdown
#' Type Curves Are Used To Analyze Pressure Buildup And Other Test Data. SPE Paper 9289.
#' @examples
#' data(agarwal_recovery)
agarwal_recovery_well_function <- function(td, par = NULL){
  res <- theis_well_function(td)
  return(res)
}
#' @title
#' agarwal_recovery_solution_initial
#' @description
#' Wrapper to function theis_solution_initial (The recovery test is analyzed using the Theis solution)
#' @param ptest A pumping_test object
#' @return
#' A list with:
#' \itemize{
#' \item a: Slope of the straight line fitted to the drawdown data using the
#' Cooper-Jacob approach
#' \item t0: Intercept of the straight line fitted to the drawdown data using the
#' Cooper-Jacob approach
#' }
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @export
#' @family agarwal_recovery functions
#' @references
#' Agarwal, R. (1980). A New Method To Account For Producing Time Effects When Drawdown
#' Type Curves Are Used To Analyze Pressure Buildup And Other Test Data. SPE Paper 9289.
#' @examples
#' data(agarwal_recovery)
#' pump_history <- 14400
#' agarwal_corr <- recovery_equivalent_time(agarwal_recovery$t, agarwal_recovery$s,
#'                                 pump_history)
#' ptest.agarwal <- pumping_test("Well1", Q = 0.02893519, r = 60, t = agarwal_corr$ta,
#' s = agarwal_corr$sa)
#' sol0 <- agarwal_recovery_solution_initial(ptest.agarwal)
agarwal_recovery_solution_initial <- function(ptest){
  res <- theis_solution_initial(ptest)
  return(res)
}
#' @title
#' agarwal_recovery_calculate_parameters
#' @description
#' Wrapper to function theis_calculate_parameters (The recovery test is analyzed using the Theis solution)
#' @param ptest A pumping_test object
#' @param par A list with the values of slope and intercept (a and t0) of the straight
#' line fitted to the drawdown data using the Cooper-Jacob approach
#' @param hydraulic Logical flag to indicate if hydraulic parameters are calculated. If False,
#' the the statistcal parameter (a and t0) are calculated.
#' @return
#' A list with the hydraulic parameters:
#' \itemize{
#' \item Tr: Transmissivity
#' \item Ss: Storage coefficient
#' \item radius_influence
#' }
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @export
#' @family agarwal_recovery functions
#' @references
#' Agarwal, R. (1980). A New Method To Account For Producing Time Effects When Drawdown
#' Type Curves Are Used To Analyze Pressure Buildup And Other Test Data. SPE Paper 9289.
#' @examples
#' data(agarwal_recovery)
#' pump_history <- 14400
#' agarwal_corr <- recovery_equivalent_time(agarwal_recovery$t, agarwal_recovery$s,
#'                                 pump_history)
#' ptest.agarwal <- pumping_test("Well1", Q = 0.02893519, r = 60, t = agarwal_corr$ta,
#' s = agarwal_corr$sa)
#' sol0 <- agarwal_recovery_solution_initial(ptest.agarwal)
#' hydr.par <- agarwal_recovery_calculate_parameters(ptest.agarwal, sol0, TRUE)
agarwal_recovery_calculate_parameters <- function(ptest, par, hydraulic = TRUE){
  res <- theis_calculate_parameters(ptest, par, hydraulic)
  return(res)
}
#' @title
#' agarwal_recovery_solution
#' @description
#' Wrapper to function theis_solution (The recovery test is analyzed using the Theis solution)
#' @param ptest A pumping_test object
#' @param a Slope of the straight line fitted to the drawdown data using the
#' Cooper-Jacob approach
#' @param t0 Intercept of the straight line fitted to the drawdown data using the
#' Cooper-Jacob approach
#' @param t Numeric vector with time
#' @return
#' A numeric vector with the calculate drawdown
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @export
#' @family agarwal_recovery functions
#' @references
#' Agarwal, R. (1980). A New Method To Account For Producing Time Effects When Drawdown
#' Type Curves Are Used To Analyze Pressure Buildup And Other Test Data. SPE Paper 9289.
#' @examples
#' data(agarwal_recovery)
#' pump_history <- 14400
#' agarwal_corr <- recovery_equivalent_time(agarwal_recovery$t, agarwal_recovery$s,
#'                                 pump_history)
#' ptest.agarwal <- pumping_test("Well1", Q = 0.02893519, r = 60, t = agarwal_corr$ta,
#' s = agarwal_corr$sa)
#' sol0 <- agarwal_recovery_solution_initial(ptest.agarwal)
#' scal <- agarwal_recovery_solution(ptest.agarwal, sol0$a, sol0$t, ptest.agarwal$t)
#' print(scal)
agarwal_recovery_solution <- function(ptest, a, t0, t){
  res <- theis_solution(ptest, a, t0, t)
  return(res)
}
#' @title
#' agarwal_recovery_WF_LT
#' @description
#' Wrapper to function theis_WF_LT (The recovery test is analyzed using the Theis solution)
#' @param p Laplace Transform parameter
#' @param ... Additional parameters
#' @return
#' A numeric vector with the calculate drawdown in the Laplace domain
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @export
#' @family agarwal_recovery functions
#' @references
#' Agarwal, R. (1980). A New Method To Account For Producing Time Effects When Drawdown
#' Type Curves Are Used To Analyze Pressure Buildup And Other Test Data. SPE Paper 9289.
#' @examples
#' coeffs <- stehfest_coeffs(8)
#' td <- logseq(-2,4,50)
#' W <- stehfest_inversion(td, coeffs, agarwal_recovery_WF_LT)
#' plot(td, W, log = "xy")
agarwal_recovery_WF_LT <- function(p, ...){
  res <- theis_WF_LT(p,...)
  return(res)
}
#' @title
#' agarwal_recovery_solution_dlogt
#' @description
#' Wrapper to function theis_solution_dlogt (The recovery test is analyzed using the Theis solution)
#' @param ptest A pumping_test object
#' @param a Slope of the straight line fitted to the drawdown data using the
#' Cooper-Jacob approach
#' @param t0 Intercept of the straight line fitted to the drawdown data using the
#' Cooper-Jacob approach
#' @param t Numeric vector with the values of time
#' @return
#' A list with
#' \itemize{
#' \item x: Numeric vector with the time values
#' \item y: Numeric vector with the derivative values
#' }
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @export
#' @family agarwal_recovery functions
#' @references
#' Agarwal, R. (1980). A New Method To Account For Producing Time Effects When Drawdown
#' Type Curves Are Used To Analyze Pressure Buildup And Other Test Data. SPE Paper 9289.
#' @examples
#' data(agarwal_recovery)
agarwal_recovery_solution_dlogt <- function(ptest, a, t0, t){
  res <- theis_solution_dlogt(ptest, a, t0, t)
  return(res)
}