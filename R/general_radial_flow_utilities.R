#' @section General_radial_flow functions:
#'
#' The functions included in this section are:
#'
#' general_radial_flow_well_function, general_radial_flow_solution_initial, general_radial_flow_calculate_parameters, general_radial_flow_solution, general_radial_flow_WF_LT, general_radial_flow_dlogt
#'
#' @docType package
#' @name pumpingtest
NULL
#' @title
#' general_radial_flow_well_function
#' @description
#' Calculates the well function of the general radial flow model
#' @param td Numeric vector with the dimensionless time
#' @param par A list with the parameters flow dimension(n), dimensionless radius (rd), and
#' the vector of stehfest coefficients (coeffs)
#' @return
#' A numeric vector with the dimensionless drawdown
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @family general_radial_flow functions
#' @export
#' @references
#' Barker, J. A. A generalized radial flow model for hydraulic tests in fractured rock
#' Water Resources Research, 1988, 24, 1796-1804.
#' @examples
#' td <- logseq(-1, 4, 50)
#' par <- list(n = 2, rd = 10, coeffs = stehfest_coeffs(8))
#' W <- general_radial_flow_well_function(td, par)
#' plot(td, W, type = "l", log = "xy", main = "GRAL RADIAL FLOW: Well function")
general_radial_flow_well_function <- function(td, par){
  #coeffs <- par$coeffs
  arg1 <- par$n
  arg2 <- par$rd
  #sd <- stehfest_inversion(td, coeffs, general_radial_flow_WF_LT, arg1 = arg1, arg2 = arg2)
  sd <- general_radial_flow_well_function_cpp(td, arg1, arg2, 0.0)
  return(sd)
}
#' @title
#' general_radial_flow_solution_initial
#' @description
#' This function calculate the initial parameters
#' @param ptest A pumping_test object
#' @return
#' A list with
#' \itemize{
#' \item a: Slope of the straight line fitted to the drawdown data using the
#' Cooper-Jacob approach
#' \item t0: Intercept of the straight line fitted to the drawdown data using the
#' Cooper-Jacob approach
#' \item n: Flow dimension (n=2)
#' }
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @export
#' @family general_radial_flow functions
#' @references
#' Barker, J. A. A generalized radial flow model for hydraulic tests in fractured rock
#' Water Resources Research, 1988, 24, 1796-1804.
#' @examples
#' data("general_radial_flow1")
#' ptest.grf1 <- pumping_test("Well1", Q= 0.02322, r = 26.2, t = general_radial_flow1$t,
#'                            s = general_radial_flow1$s)
#' par <- list(rw = 0.1, rc = 0.1, rd = 1)
#' ptest.grf1$additional_parameters <- par
#' sol0 <-  general_radial_flow_solution_initial(ptest.grf1)
#' print(sol0)
general_radial_flow_solution_initial <- function(ptest){
  if(class(ptest) != 'pumping_test'){
    stop('A pumping_test object is required as input')
  }
  n1 <- floor(length(ptest$t)/3)
  endp <- length(ptest$t)
  ptest1 <- as.data.frame(cbind(ptest$t[n1:endp],ptest$s[n1:endp]))
  names(ptest1) <- c('t', 's')
  initial <- lm(s ~ log10(t), data = ptest1)
  t0 <- 10^(-initial$coefficients[1]/initial$coefficients[2]);
  t0 <- t0*(ptest$additional_parameters$rw/ptest$r)^2;
  res <- list(a = unname(initial$coefficients[2]), t0 = unname(t0), n = 2)
  return(res)
}
#' @title
#' general_radial_flow_calculate_parameters
#' @description
#' Calculates the hydraulic parameters form general_radial_flow model
#' @param ptest A pumping_test object
#' @param par A list with the slope and intercept (a and t0) of the straight line
#' fitted to the drawdown data using the Cooper-Jacob approach
#' @param hydraulic Logical flag to indicate if hydraulic parameters are calculated.
#' If False, the the statistcal parameter (a and t0) are calculated.
#' @return
#' A list with
#' \itemize{
#' \item Tr: Transmissivity
#' \item Ss: Storage coefficient
#' }
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @family general_radial_flow functions
#' @export
#' @references
#' Barker, J. A. A generalized radial flow model for hydraulic tests in fractured rock
#' Water Resources Research, 1988, 24, 1796-1804.
#' @examples
#' data("general_radial_flow1")
#' ptest.grf1 <- pumping_test("Well1", Q= 0.02322, r = 26.2, t = general_radial_flow1$t,
#'                            s = general_radial_flow1$s)
#' par <- list(rw = 0.1, rc = 0.1, rd = 1)
#' ptest.grf1$additional_parameters <- par
#' sol0 <-  general_radial_flow_solution_initial(ptest.grf1)
#' grf1.par <- general_radial_flow_calculate_parameters(ptest.grf1, sol0)
#' print(grf1.par)
general_radial_flow_calculate_parameters <- function(ptest, par, hydraulic = TRUE){
  if(class(ptest) != 'pumping_test'){
    stop('A pumping_test object is required as input')
  }
  Q <- ptest$Q
  r <- ptest$r
  if(is.null(ptest$additional_parameters)){
    stop('The additional_parameters list is empty')
  }
  rw <- ptest$additional_parameters$rw
  if(hydraulic){
    a <- par$a
    t0 <- par$t0
    pos_valid_a <- abs(a) < 1.0e-15
    pos_valid_t0 <- abs(t0) < 1.0e-15
    if( sum(pos_valid_a) > 0 | sum(pos_valid_t0) > 0){
      stop('a or t0 or both are close to numerical precision')
    }
    # Equivalent Cylindrical Transmissivity
    Tr <- log(10)*Q/4/pi/a
    Ss <- 2.2458394*Tr*t0/rw^2
    res <- list(Tr = Tr, Ss = Ss, n = par$n)
  }
  else {
    Tr <- par$Tr
    Ss <- par$Ss
    a <- log(10)*Q/4/pi/Tr
    t0 <- (Ss*rw^2)/(2.2458394*Tr)
    res <- list(a = a, t0 = t0, n = par$n)
  }
  return(res)
}
#' @title
#' general_radial_flow_solution
#' @description
#' Function to calculate the drawdown using a general_radial_flow solution
#' @param ptest A pumping_test object
#' @param a Slope of the straight line fitted to the drawdown data using the
#' Cooper-Jacob approach
#' @param t0 Intercept of the straight line fitted to the drawdown data using the
#' Cooper-Jacob approach
#' @param n Flow dimension
#' @param t Numeric vector with time
#' @return
#' A numeric vector with the calculated drawdown
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @export
#' @family general_radial_flow functions
#' @references
#' Barker, J. A. A generalized radial flow model for hydraulic tests in fractured rock
#' Water Resources Research, 1988, 24, 1796-1804.
#' @examples
#' data("general_radial_flow1")
#' ptest.grf1 <- pumping_test("Well1", Q= 0.02322, r = 26.2, t = general_radial_flow1$t,
#'                            s = general_radial_flow1$s)
#' par <- list(rw = 0.1, rc = 0.1, rd = 1)
#' ptest.grf1$additional_parameters <- par
#' sol0 <-  general_radial_flow_solution_initial(ptest.grf1)
#' sol <- general_radial_flow_solution(ptest.grf1, sol0$a, sol0$t0, sol0$n, ptest.grf1$t)
#' print(sol)
general_radial_flow_solution <- function(ptest, a, t0, n, t){
  Q <- ptest$Q
  r <- ptest$r
  rw <- ptest$additional_parameters$rw
  dimensionless_r <- r/rw
  coeffs <- stehfest_coeffs(8)
  td <- t/2.2458/abs(t0)
  par <- list(rd = dimensionless_r, n = n, coeffs = coeffs)
  grf_well_function <- general_radial_flow_well_function(td, par)
  res <- 0.868588963806504*a*grf_well_function
  return(res)
}
#' @title
#' general_radial_flow_WF_LT
#' @description
#' Function to calculate the general_radial_flow well function in the Laplace domain
#' @param p Laplace Transform parameter
#' @param arg1 Flow dimension (n)
#' @param arg2 Dimensionless radius (rd)
#' @return
#' A numeric vector with the calculated drawdown in the Laplace domain
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @export
#' @family general_radial_flow functions
#' @references
#' Barker, J. A. A generalized radial flow model for hydraulic tests in fractured rock
#' Water Resources Research, 1988, 24, 1796-1804.
#' @examples
#' p <- 0.5
#' sd <- general_radial_flow_WF_LT(p, 2.0, 10)
#' print(sd)
general_radial_flow_WF_LT <- function(p, arg1, arg2){
  n <- arg1
  rd <- arg2
  sp <- sqrt(p)
  sd <- rd^(2-n)*(rd^2.*p/4)^(n/4-0.5)*besselK(rd*sp,(n/2-1))/p/gamma(n/2)
  return(sd)
}
#' @title
#' general_radial_flow_solution_dlogt
#' @description
#' This function calculates the derivative of the drawdown with respect to the
#' derivative of the log of time
#' @param ptest A pumping_test object
#' @param a Slope of the straight line fitted to the drawdown data using the
#' Cooper-Jacob approach
#' @param t0 Intercept of the straight line fitted to the drawdown data using the
#' Cooper-Jacob approach
#' @param n Flow dimension
#' @param t Numeric vector with time
#' @return
#' A numeric vector with the derivative of drawdown
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @export
#' @family general_radial_flow functions
#' @references
#' Barker, J. A. A generalized radial flow model for hydraulic tests in fractured rock
#' Water Resources Research, 1988, 24, 1796-1804.
#' @examples
#' data("general_radial_flow1")
#' ptest.grf1 <- pumping_test("Well1", Q= 0.02322, r = 26.2, t = general_radial_flow1$t,
#'                            s = general_radial_flow1$s)
#' par <- list(rw = 0.1, rc = 0.1, rd = 1)
#' ptest.grf1$additional_parameters <- par
#' sol0 <-  general_radial_flow_solution_initial(ptest.grf1)
#' sol_dlogt <- general_radial_flow_solution_dlogt(ptest.grf1, sol0$a, sol0$t0, sol0$n, ptest.grf1$t)
#' print(sol_dlogt)
general_radial_flow_solution_dlogt <- function(ptest, a, t0, n, t){
  sb <- general_radial_flow_solution(ptest, a, t0, n, t)
  dl_sb <- log_derivative_central(t, sb)
  res <- dl_sb$y
  return(res)
}