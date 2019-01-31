#' @section Gringarten functions:
#'
#' These functions are used in the estimation of pumping tests in a single
#' fracture aquifers.
#'
#' The functions included in this section are:
#'
#' gringarten_well_function, gringarten_solution_initial, gringarten_calculate_parameters, gringarten_solution, gringarten_solution_dlogt
#'
#' @docType package
#' @name pumpingtest
NULL
#' @title
#' gringarten_well_function
#' @description
#' Function to calculate the drawdown using the Gringarten Solution
#' @param td A numeric vector with the dimensionless time
#' @param par A list with the values of slope and intercept (a and t0) of the straight
#' line fitted to the drawdown data using the Cooper-Jacob approach
#' @return
#' A numeric vector with the calculated drawdown
#' @importFrom pracma erf expint
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @export
#' @references
#' Gringarten A. C., Ramey H. J., & Raghavan R., 1975.
#' Applied Pressure Analysis for Fractured Wells, Petroleum Transactions,
#' AIME follows page 784.
#' @family gringarten functions
gringarten_well_function <- function(td, par){
  a <- par$a
  t0 <- par$t0
  u <- td/t0
  s <- a*(2.0*sqrt(pi*u)*erf(1./(2.0*sqrt(u)))+expint(1./(4.*u)))
  return(s)
}
#' @title
#' gringarten_solution_initial
#' @description
#' Function to calculate the initial parameters of the Gringarten Solution
#' @param ptest A pumping_test object
#' @return
#' A list with
#' \itemize{
#' \item a: Slope of the straight line fitted to the drawdown data using the
#' Cooper-Jacob approach
#' \item t0: Intercept  of the straight line fitted to the drawdown data using the
#' Cooper-Jacob approach
#' }
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @export
#' @references
#' Gringarten A. C., Ramey H. J., & Raghavan R., 1975.
#' Applied Pressure Analysis for Fractured Wells, Petroleum Transactions,
#' AIME follows page 784.
#' @family gringarten functions
gringarten_solution_initial <- function(ptest){
  if(class(ptest) != 'pumping_test'){
    stop('A pumping_test object is required as input')
  }
  endp <- length(ptest$t)
  nt <- round(endp/4)
  t <- ptest$t[nt:endp]
  s <- ptest$s[nt:endp]
  ptest1 <- pumping_test(ptest$id, Q = ptest$Q, r = ptest$r, t = t, s = s)
  res <- cooper_jacob_solution_initial(ptest1)
  return(res)
}
#' @title
#' gringarten_calculate_parameters
#' @description
#' Calculate the hydraulic parameters of the Gringarten solution
#' @param ptest A pumping_test object
#' @param par A list with the values of slope and intercept (a and t0) of the straight
#' line fitted to the drawdown data using the Cooper-Jacob approach
#' @param hydraulic Logical flag to indicate if hydraulic parameters are calculated.
#' If False, the the statistcal parameter (a and t0) are calculated.
#' @return
#' A list with
#' \itemize{
#' \item Tr: Transmissivity
#' \item Sxf2: ????
#' }
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @export
#' @references
#' Gringarten A. C., Ramey H. J., & Raghavan R., 1975.
#' Applied Pressure Analysis for Fractured Wells, Petroleum Transactions,
#' AIME follows page 784.
#' @family gringarten functions
gringarten_calculate_parameters <- function(ptest, par, hydraulic = TRUE){
  if(class(ptest) != 'pumping_test'){
    stop('A pumping_test object is required as input')
  }
  if(class(par) != 'list'){
    stop('A list is required as input')
  }
  Q <- ptest$Q
  if(hydraulic){
    a <- par$a
    t0 <- par$t0
    pos_valid_a <- abs(a) < 1.0e-15
    pos_valid_t0 <- abs(t0) < 1.0e-15
    if( sum(pos_valid_a) > 0 | sum(pos_valid_t0) > 0){
      stop('a or t0 or both are close to numerical precision')
    }
    Tr <- Q/(4*pi*a)
    Sxf2 <- t0/Tr
    res <- list(Tr = Tr, Sxf2 = Sxf2)
  }
  else {
    Tr <- par$Tr
    Sxf2 <- par$Sxf2
    a <- Q/(4*pi*Tr)
    t0 <- Sxf2*Tr
    res <- list(a = a, t0 = t0)
  }
  return(res)
}
#' @title
#' gringarten_solution
#' @description
#' Calculate the drawdown using the Gringarten Solution
#' @param ptest A pumping_test object
#' @param a Slope  of the straight line fitted to the drawdown data using
#' the Cooper-Jacob approach
#' @param t0 Intercept of the straight line fitted to the drawdown data using
#' the Cooper-Jacob approach
#' @param t A numeric vector with the time values
#' @return
#' A numeric vector with the calculated drawdown
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @references
#' Gringarten A. C., Ramey H. J., & Raghavan R., 1975.
#' Applied Pressure Analysis for Fractured Wells, Petroleum Transactions,
#' AIME follows page 784.
#' @family gringarten functions
gringarten_solution <- function(ptest, a, t0, t){
  if(class(ptest) != 'pumping_test'){
    stop('A pumping_test object is required as input')
  }
  par <- list(a = a, t0 = t0)
  res <- gringarten_well_function(t, par)
  return(res)
}
#' @title
#' gringarten_solution_dlogt
#' @description
#' Function to calculate the derivative of drawdown with respect to logarithm of time using the
#' Gringarten Solution
#' @param ptest A pumping_test object
#' @param a Slope  of the straight line fitted to the drawdown data using
#' the Cooper-Jacob approach
#' @param t0 Intercept of the straight line fitted to the drawdown data using
#' the Cooper-Jacob approach
#' @param t A numeric vector with the time values
#' @return
#' A numeric vector with the derivative of the calculated drawdown
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @export
#' @references
#' Gringarten A. C., Ramey H. J., & Raghavan R., 1975.
#' Applied Pressure Analysis for Fractured Wells, Petroleum Transactions,
#' AIME follows page 784.
#' @family gringarten functions
gringarten_solution_dlogt <- function(ptest, a, t0, t){
  if(class(ptest) != 'pumping_test'){
    stop('A pumping_test object is required as input')
  }
  par <- list(a = a, t0 = t0)
  s <- gringarten_well_function(t, par)
  s_dlogt <- log_derivative_central(t, s)
  return(s_dlogt$y)
}