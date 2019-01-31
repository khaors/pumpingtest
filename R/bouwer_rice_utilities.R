#' @section Bouwer Rice functions:
#'
#' These functions are used in the estimation of hydraulic parameters using a slug
#' test
#'
#' The functions included in this section are:
#'
#' Acoeff, Bcoeff, Ccoeff, bouwer_rice_well_function, bouwer_rice_solution_initial, bouwer_rice_calculate_parameters, bouwer_rice_solution, bouwer_rice_solution_dlogt
#'
#' @docType package
#' @name pumpingtest
NULL
#' @title
#' Acoeff
#' @description
#' Function to calculate the A coefficient of the Bouwer and Rice solution.
#' @param x A numeric value
#' @return
#' This function returns the value of the A coefficient
#' @export
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @family bouwer_rice functions
Acoeff <- function(x){
  A <- 0
  if(x < 2.554422663){
    A <- 1.638445671 + 0.166908063 * x + 0.000740459 *
      exp(6.17105281 * x - 1.054747686 * x * x)
  }
  else{
    A <- 11.00393028 - 170.7752217 * exp(-1.509639982 *x)
  }
  return(A)
}
#' @title
#' Bcoeff
#' @description
#' Function to calculate the A coefficient of the Bouwer and Rice solution.
#' @param x A numeric value
#' @return
#' This function returns the value of the B coefficient
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @export
#' @family bouwer_rice functions
Bcoeff <- function(x){
  B <- 0
  if(x < 2.596774459){
    B <- 0.174811819 + 0.060059188 * x + 0.007965502 *
      exp(2.053376868 * x - 0.007790328 * x * x)
  }
  else {
    B <- 4.133124586 - 93.06136936 * exp(-1.435370997 * x)
  }
  return(B)
}
#' @title
#' Ccoeff
#' @description
#' Function to calculate the C coefficient of the Bouwer and Rice solution.
#' @param x A numeric value
#' @return
#' This function returns the value of the C coefficient
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @export
#' @family bouwer_rice functions
Ccoeff <- function(x){
  C <- 0.0
  if(x < 2.200426117){
    C <- 0.074711376 + 1.083958569 * x + 0.00557352 *
      exp(2.929493814  * x - 0.001028433 * x * x)
  }
  else {
    C <- 15.66887372 - 178.4329289 * exp(-1.322779744 * x)
  }
  return(C)
}
#' @title
#' bouwer_rice_well_function
#' @description
#' Function to calculate the well function of the Bower-Rice solution
#' @param td A numeric vector with the dimensionless time
#' @param par A list with the value of the t0
#' @return
#' A numeric vector with the normalized drawdown
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @export
#' @family bouwer_rice functions
bouwer_rice_well_function <- function(td, par){
  t0 <- par$t0
  s <- exp(-t/t0)
  return(s)
}
#' @title
#' bouwer_rice_initial_solution
#' @description
#' Function to calculate the initial value of t0 that defines the Bower-Rice solution
#' @param ptest A pumping_test object
#' @return
#' A list with the value of t0
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @export
#' @family bouwer_rice functions
bouwer_rice_initial_solution <- function(ptest){
  if(class(ptest) != 'pumping_test'){
    stop('ERROR A pumping_test object is required as input')
  }
  t0 <- pracma::interp1(rev(ptest$s), rev(ptest$t), xi = 0.37, method = 'linear')
  res <- list(t0 = t0)
  return(res)
}
#' @title
#' bower_rice_solution
#' @description
#' Function to calculate the drawdown during a slug test using the Bower-Rice solution.
#' @param ptest A pumping_test object
#' @param t0 A numeric value with the t0
#' @param t A numeric vector with the time
#' @return
#' A numeric vector with the calculated drawdown
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @export
#' @family bouwer_rice functions
bouwer_rice_solution <- function(ptest, t0, t){
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