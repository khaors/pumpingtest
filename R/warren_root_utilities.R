#' @section Warren and Root functions:
#'
#' These functions are used in the estimation of pumping tests in dual porosity
#' aquifers.
#'
#' The functions included in this section are:
#'
#' warren_root_well_function, warren_root_solution_initial, warren_root_calculate_parameters, warren_root_solution, warren_root_WF_LT, warren_root_solution_dlogt
#'
#' @docType package
#' @name pumpingtest
NULL
#' @title
#' warren_root_well_function
#' @description
#' Calculates the well function of the Warrent and Root model
#' @param td A numeric vector with dimensionless time
#' @param par List with parameters sigma, interflow porosity (lambda) and Stehfest coefficients
#' @return
#' A numeric vector with the dimensionless drawdown
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @export
#' @family warren_root functions
#' @references
#' Warren, J. & Root, P. The Behavior of Naturally Fractured Reservoirs. Society of
#' Petroleum Engineers Journal., 1963, 3.
#' @examples
#' td <- logseq(-1, 4, 50)
#' par <- list(sigma = 1, lambda = 1, coeffs = stehfest_coeffs(8))
#' W <- warren_root_well_function(td, par)
warren_root_well_function <- function(td, par){
  sigma <- par$sigma
  lambda <- par$lambda
  #coeffs <- par$coeffs
  #sd <- stehfest_inversion(td, coeffs, warren_root_WF_LT, arg1 = sigma, arg2 = lambda)
  sd <- warren_root_well_function_cpp(td, sigma, lambda, 0.0)
  return(sd)
}
#' @title
#' warren_root_solution_initial
#' @description
#' Calculates the initial values of parameters of Warren and Root solution
#' @param ptest A pumping_test object
#' @return
#' A list with
#' \itemize{
#' \item a: Slope
#' \item t0: Intercept of the early time asymptote
#' \item t1: Intercept of the late time asymptote
#' \item tm: time of the minimum derivative
#' }
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @export
#' @family warren_root functions
#' @references
#' Warren, J. & Root, P. The Behavior of Naturally Fractured Reservoirs. Society of
#' Petroleum Engineers Journal., 1963, 3.
warren_root_solution_initial <- function(ptest){
  if(class(ptest) != 'pumping_test'){
    stop('A pumping_test object is required as input')
  }
  t <- ptest$t
  s <- ptest$s
  s_dlogt <- log_derivative_spline(t, s, n = 40)
  endp <- length(s_dlogt$y)
  dd <- mean(s_dlogt$y[(endp-3):endp]);
  a <- log(10)*dd;
  t0 <- t[1]/exp(s[1]/dd)
  t1 <- t[endp]/exp(s[endp]/dd)
  pos_min <- which.min(s_dlogt$y)
  tm <- s_dlogt$x[pos_min]
  res <- list(a = a, t0 = t0, t1 = t1, tm = tm)
  return(res)
}
#' @title
#' warren_root_calculate_parameters
#' @description
#' A function to calculates the hydraulic parameters for a Warren and Root model
#' @param ptest A pumping_test object
#' @param par A list with the parameters a, t0, t1, tm
#' @param hydraulic Logical flag to indicate if hydraulic parameters are calculated.
#' If False, the the statistcal parameter (a and t0) are calculated.
#' @return
#' A list with the parameters:
#' \itemize{
#' \item Tr: Transmissivity Tf
#' \item Ss: Storativity Sf
#' \item Sm: Storativity Sm
#' \item lambda: Interflow porosity
#' }
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @export
#' @references
#' Warren, J. & Root, P. The Behavior of Naturally Fractured Reservoirs. Society of
#' Petroleum Engineers Journal., 1963, 3.
warren_root_calculate_parameters <- function(ptest, par, hydraulic = TRUE){
  if(class(ptest) != 'pumping_test'){
    stop('A pumping_test object is required as input')
  }
  if(class(par) != 'list'){
    stop('a list is required as input')
  }
  Q <- ptest$Q
  r <- ptest$r
  if(hydraulic){
    a <- par$a
    t0 <- par$t0
    t1 <- par$t1
    tm <- par$tm
    pos_valid_a <- abs(a) < 1.0e-15
    pos_valid_t0 <- abs(t0) < 1.0e-15
    pos_valid_t1 <- abs(t1) < 1.0e-15
    pos_valid_tm <- abs(tm) < 1.0e-15
    if( sum(pos_valid_a) > 0 | sum(pos_valid_t0) > 0 |
        sum(pos_valid_t1) > 0| sum(pos_valid_tm) >0 ){
      stop('a or t0 or t1 or tm or all of them are close to numerical precision')
    }
    Tf <- 0.1832339*Q/a
    Sf <- 2.245839*Tf*t0/r^2
    Sm <- 2.245839*Tf*t1/r^2-Sf
    sigma <- (t1-t0)/t0
    lambda <- 2.2458394*t0*log(t1/t0)/tm
    res <- list(Tr = Tf, Ss = Sf, Sm = Sm, lambda = lambda)
  }
  else{
    Tr <- par$Tr
    Ss <- par$Ss
    Sm <- par$Sm
    lambda <- par$lambda
    a <- 0.1832339*Q/Tr
    t0 <- (Ss*r^2)/(2.245839*Tr)
    t1 <- (r^2)*(Sm+Ss)/(2.245839*Tr)
    tm <- 2.2458394*t0*log(t1/t0)/lambda
    res <- list(a = a, t0 = t0, t1 = t1, tm = tm)
  }
  return(res)
}
#' @title
#' warren_root_solution
#' @description
#' Function to calculate the drawdown using the Warren Root model
#' @param ptest A pumping_test object
#' @param a Slope of the Cooper-Jacob solution r
#' @param t0 Intercept with the time axis for the early time asymptote
#' @param t1 Intercept with the time axis for the late time asymptote
#' @param tm Time of the minimum derivative
#' @param t A numeric vector with the time values
#' @return
#' A numeric vector with the calculated drawdown
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @export
#' @family warren_root functions
#' @references
#' Warren, J. & Root, P. The Behavior of Naturally Fractured Reservoirs. Society of
#' Petroleum Engineers Journal., 1963, 3.
warren_root_solution <- function(ptest, a, t0, t1, tm, t){
  if(class(ptest) != 'pumping_test'){
    stop('A pumping_test object is required as input')
  }
  td <- 0.445268*t/t0
  sigma <- (t1-t0)/t0
  lambda <- 2.2458394*t0*log(t1/t0)/tm
  par <- list(sigma = sigma, lambda = lambda, coeffs = ptest$coeffs)
  sd <- warren_root_well_function(td, par)
  s <- 0.868589*a*sd
  return(s)
}
#' @title
#' warren_root_WF_LT
#' @description
#' Calculates the well function of the Warrent and Root model
#' @param p Laplace tranform parameter
#' @param arg1 Value of parameter sigma
#' @param arg2 Value of parameter lambda
#' @return
#' A numeric vector with the dimensionless drawdown in the Laplace domain
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @export
#' @family warren_root functions
#' @references
#' Warren, J. & Root, P. The Behavior of Naturally Fractured Reservoirs. Society of
#' Petroleum Engineers Journal., 1963, 3.
#' @examples
#' p <- 0.5
#' sigma <- 1
#' lambda <- 1
#' s <- warren_root_WF_LT(p, sigma, lambda)
warren_root_WF_LT <- function(p, arg1, arg2){
  sigma <- arg1
  lambda <- arg2
  s <- (1./p)*besselK(sqrt(p+(lambda*sigma*p)/(sigma*p+lambda)), 0)
  return(s)
}
#' @title
#' warren_root_solution_dlogt
#' @description
#' Function to calculate the derivative of the drawdown with respect to the logarithm of
#' time for the Warren Root model
#' @param ptest A pumping_test object
#' @param a Slope of the Cooper-Jacob solution r
#' @param t0 Intercept with the time axis for the early time asymptote
#' @param t1 Intercept with the time axis for the late time asymptote
#' @param tm Time of the minimum derivative
#' @param t A numeric vector with the time values
#' @return
#' A numeric vector with the calculated drawdown
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @export
#' @family warren_root functions
#' @references
#' Warren, J. & Root, P. The Behavior of Naturally Fractured Reservoirs. Society of
#' Petroleum Engineers Journal., 1963, 3.
warren_root_solution_dlogt <- function(ptest, a, t0, t1, tm, t){
  if(class(ptest) != 'pumping_test'){
    stop('A pumping_test object is required as input')
  }
  sb <- warren_root_solution(ptest, a, t0, t1, tm, t)
  dl_sb <- log_derivative_central(t, sb)
  res <- dl_sb$y
  return(res)
}