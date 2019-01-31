#' @section Hvorslev functions:
#'
#' These functions are used in the estimation of the aquifer parameters using a
#' slug test.
#'
#' The functions included in this section are:
#'
#' hvorslev_well_function, hvorslev_solution_initial, hvorslev_calculate_parameters, hvorslev_solution, hvorslev_solution_dlogt
#'
#' @docType package
#' @name pumpingtest
NULL
#' @title
#' hvorslev_well_function
#' @description
#' Function calculate the drawdown using the Hvorslev Solution during a slug test.
#' @param td A numeric vector with the dimensionless time
#' @param par A list with the values of intercept t0 of the straight line fitted to the
#' drawdown data.
#' @return
#' A numeric vector with the calculated drawdown
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @export
#' @family hvorslev functions
#' @references
#' Hvorslev, M.J., 1951. Time Lag and Soil Permeability in Ground-Water Observations,
#' Bull. No. 36, Waterways Exper. Sta. Corps of Engrs, U.S. Army, Vicksburg, Mississippi,
#' pp. 1-50.
hvorslev_well_function <- function(td, par){
  t0 <- par$t0
  s <- exp(-t/t0)
  return(s)
}
#' @title
#' hvorslev_initial_solution
#' @description
#' Function calculate the drawdown using the Hvorslev Solution during a slug test.
#' @param ptest A pumpingtest object
#' @return
#' A list with the time t37.
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @export
#' @family hvorslev functions
#' @references
#' Hvorslev, M.J., 1951. Time Lag and Soil Permeability in Ground-Water Observations,
#' Bull. No. 36, Waterways Exper. Sta. Corps of Engrs, U.S. Army, Vicksburg, Mississippi,
#' pp. 1-50.
hvorslev_initial_solution <- function(ptest){
  if(class(ptest) != 'pumping_test'){
    stop('ERROR A pumping_test object is required as input')
  }
  t0 <- pracma::interp1(rev(ptest$s), rev(ptest$t), xi = 0.37, method = 'linear')
  res <- list(t0 = t0)
  return(res)
}
#' @title
#' hvorslev_calculate_parameters
#' @description
#' Function calculate the drawdown using the Hvorslev Solution during a slug test.
#' @param ptest A pumpingtest object
#' @param par A list with the fit parameters
#' @param hydraulic A logical flag indicating if hydraulic parameters are required.
#' @return
#' A list with the time t37.
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @export
#' @family hvorslev functions
#' @references
#' Hvorslev, M.J., 1951. Time Lag and Soil Permeability in Ground-Water Observations,
#' Bull. No. 36, Waterways Exper. Sta. Corps of Engrs, U.S. Army, Vicksburg, Mississippi,
#' pp. 1-50.
hvorslev_calculate_parameters <- function(ptest, par, hydraulic = TRUE){
  if(class(ptest) != 'pumping_test'){
    stop('ERROR A pumping_test object is required as input')
  }
  res <- NULL
  if(hydraulic){
    Le <- ptest$additional_parameters$Le
    rw <- ptest$additional_parameters$rw #Radius well screen
    rc <- ptest$additional_parameters$rc #Radius casing
    t0 <- par$t0
    K <- (rc^2)*log(Le/rw)/(2*Le*t0)
    res <- list(K = K)
  }
  else{
    Le <- ptest$additional_parameters$Le
    rw <- ptest$additional_parameters$rw #Radius well screen
    rc <- ptest$additional_parameters$rc #Radius casing
    K <- par$K
    t0 <- (rc^2)*log(Le/rw)/(2*Le*K)
    res <- list(t0 = t0)
  }
  return(res)
}
#' @title
#' hvorslev_solution
#' @description
#' Function calculate the drawdown using the Hvorslev Solution during a slug test.
#' @param ptest A pumpingtest object
#' @param t0 Time for which the normalized drawdown is equal to 0.37
#' @param t A numeric vector with the time
#' @return
#' A numeric vector with the calculated drawdown
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @export
#' @family hvorslev functions
#' @references
#' Hvorslev, M.J., 1951. Time Lag and Soil Permeability in Ground-Water Observations,
#' Bull. No. 36, Waterways Exper. Sta. Corps of Engrs, U.S. Army, Vicksburg, Mississippi,
#' pp. 1-50.
hvorslev_solution <- function(ptest, t0, t){
  if(class(ptest) != 'pumping_test'){
    stop('ERROR A pumping_test object is required as input')
  }
  Le <- ptest$additional_parameters$Le
  rw <- ptest$additional_parameters$rw #Radius well screen
  rc <- ptest$additional_parameters$rc #Radius casing
  K <- (rc^2)*log(Le/rw)/(2*Le*t0)
  s <- exp(-(2*K*Le*t)/(log(200)*rc^2) )
  return(s)
}