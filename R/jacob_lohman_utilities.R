#' @section jacob_lohman functions:
#' There are two types of functions: F and G. The Jacob_Lohman F functions are used to calculate
#' drawdown at specific distance from the well, while the Jacob_Lohman G Functions are used to
#' calculate the discharge at the well or a given distance from it.
#'
#' The F functions are:
#'
#' jacob_lohman_F_well_function, jacob_lohman_F_solution_initial, jacob_lohman_F_calculate_parameters, jacob_lohman_F_solution, jacob_lohman_F_solution_dlogt, jacob_lohman_F_WF_LT
#'
#' The G Functions are
#'
#' jacob_lohman_G_well_function, jacob_lohman_G_solution_initial, jacob_lohman_G_calculate_parameters, jacob_lohman_G_solution, jacob_lohman_G_solution_dlogt, jacob_lohman_G_WF_LT
#'
#' @docType package
#' @name pumpingtest
NULL
#' @title
#' jacob_lohman_F_well_function
#' @description
#' This calculates the well function of the Jacob-Lohman solution for a constant head aquifer
#' (artesian aquifer). This function calculates the drawdown.
#' @param td Numeric vector with the Dimensionless time
#' @param par A list with the parameters
#' @return
#' A numeric vector with the dimensionless drawdowns
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @export
#' @family jacob_lohman functions
#' @references
#' Jacob, C. E. & Lohman, S. W. Nonsteady Flow to a Well Of Constant Drawdown. American
#' Geophysical Union, 1952, 33, 10.
#' @examples
#' td <- logseq(2, 6, 100)
#' par <- list(rho = 5, coeffs = stehfest_coeffs(8))
#' sd <- jacob_lohman_F_well_function(td, par)
jacob_lohman_F_well_function <- function(td, par){
  coeffs <- par$coeffs
  arg1 <- par$rho
  sd <- stehfest_inversion(td, coeffs, jacob_lohman_WF_F_LT, arg1 = arg1)
  return(sd)
}
#' @title
#' jacob_lohman_F_solution_initial
#' @description
#' Function to calculate the initial parameters of aconstant head aquifer (artesian aquifer)
#' using the Jacob-Lohman solution
#' @param ptest A pumping_test object
#' @return
#' A list with
#' \itemize{
#'  \item a : Slope of the straight line fitted to the drawdown data using the
#' Cooper-Jacob approach
#'  \item t0 : Intercept of the straight line fitted to the drawdown data using the
#' Cooper-Jacob approach
#' }
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @export
#' @family jacob_lohman functions
#' @references
#' Jacob, C. E. & Lohman, S. W. Nonsteady Flow to a Well Of Constant Drawdown. American
#' Geophysical Union, 1952, 33, 10.
#' @examples
#' data(jacob_lohman)
#' s <- vector("numeric", length(jacob_lohman$t))
#' s[1:length(jacob_lohman$t)] <- 0
#' ptest <- pumping_test("Well1", Q = jacob_lohman$q, r = 0.84, t = jacob_lohman$t, s = s)
#' ptest$additional_parameters <- list(rw = 0.084, s0 = 28.142)
#' sol0 <- jacob_lohman_F_solution_initial(ptest)
#' print(sol0)
jacob_lohman_F_solution_initial <- function(ptest){
  if(class(ptest) != 'pumping_test'){
    stop('A pumping_test object is required as input')
  }
  t <- ptest$t
  n <- ceiling(length(t)/3);
  Q <- ptest$Q
  nq <- length(Q)
  if(nq != length(t)){
    stop('Discharge values for each time are required')
  }
  endp <- length(t)
  t <- t[n:endp]
  s <- 1./Q[n:endp]
  jl.df <- data.frame(t,s)
  jl.lm <- lm(s ~ log10(t), data = jl.df)
  t1 <- 10^(-jl.lm$coefficients[1]/jl.lm$coefficients[2])
  parameters <- list(a = unname(jl.lm$coefficients[2]), t0 = ptest$t[1])
  res <- parameters
  return(res)
}
#' @title
#' jacob_lohman_F_calculate_parameters
#' @description
#' Calculate hydraulic parameters of the Jacob-Lohman Solution
#' @param ptest A pumping_test object
#' @param par A list with the Slope (a) and the Intercept (t0) of the straight line fitted to
#' the drawdown data using the Cooper-Jacob approach
#' @param hydraulic Logical flag to indicate if hydraulic parameters are calculated. If False,
#' the the statistcal parameter (a and t0) are calculated.
#' @return
#' A list with
#' \itemize{
#' \item Tr: Transmisivity.
#' \item Ss: The Storage coefficient
#' \item radius_influence
#' }
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @family jacob_lohman functions
#' @export
#' @references
#' Jacob, C. E. & Lohman, S. W. Nonsteady Flow to a Well Of Constant Drawdown. American
#' Geophysical Union, 1952, 33, 10.
#' @examples
#' data(jacob_lohman)
#' s <- vector("numeric", length(jacob_lohman$t))
#' s[1:length(jacob_lohman$t)] <- 0
#' ptest <- pumping_test("Well1", Q = jacob_lohman$q, r = 0.84, t = jacob_lohman$t, s = s)
#' ptest$additional_parameters <- list(rw = 0.084, s0 = 28.142)
#' sol0 <- jacob_lohman_F_solution_initial(ptest)
#' jl.par <- jacob_lohman_F_calculate_parameters(ptest, sol0)
#' print(jl.par)
jacob_lohman_F_calculate_parameters <- function(ptest, par, hydraulic = TRUE){
  if(class(ptest) != 'pumping_test'){
    stop('A pumping_test object is required as input')
  }
  endp <- length(t)
  t <- ptest$t
  Q <- ptest$Q
  Q1 <- ptest$Q[1]
  r <- ptest$r
  if(hydraulic){
    a <- par$a
    t0 <- par$t0
    pos_valid_a <- abs(a) < 1.0e-15
    pos_valid_t0 <- abs(t0) < 1.0e-15
    if( sum(pos_valid_a) > 0 | sum(pos_valid_t0) > 0){
      stop('a or t0 or both are close to numerical precision')
    }
    #Compute the transmissivity, storativity and radius of influence
    Tr <- 0.1832339*Q1/a
    Ss <- 2.2458394*Tr*t0/r^2
    Ri <- 2.0*sqrt(Tr*t[endp]/Ss)
    res <- list(Tr = Tr, Ss = Ss, radius_influence = Ri)
  }
  else {
    Tr <- par$Tr
    Ss <- par$Ss
    a <- 0.1832339*Q1/Tr
    t0 <- (Ss*r^2)/(2.2458394*Tr)
    res <- list(a = a, t0 = t0)
  }
  return(res)
}
#' @title
#' jacob_lohman_F_solution
#' @description
#' Function to calculate the drawdown in a constant rate test using the Jacob-Lohman Solution
#' @param ptest A pumping_test object
#' @param a Slope of the straight line fitted to the drawdown data using the
#' Cooper-Jacob approach
#' @param t0 Intercept of the straight line fitted to the drawdown data using the
#' Cooper-Jacob approach
#' @param t Numeric vector with the times at which measurements were taken
#' @return
#' A numeric vector with the values of the drawdown
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @export
#' @family jacob_lohman functions
#' @references
#' Jacob, C. E. & Lohman, S. W. Nonsteady Flow to a Well Of Constant Drawdown. American
#' Geophysical Union, 1952, 33, 10.
#' @examples
#' data(jacob_lohman)
#' s <- vector("numeric", length(jacob_lohman$t))
#' s[1:length(jacob_lohman$t)] <- 0
#' ptest <- pumping_test("Well1", Q = jacob_lohman$q, r = 0.84, t = jacob_lohman$t, s = s)
#' ptest$additional_parameters <- list(rw = 0.084, s0 = 28.142)
#' sol0 <- jacob_lohman_F_solution_initial(ptest)
#' sol <- jacob_lohman_F_solution(ptest, sol0$a, sol0$t0, ptest$t)
#' print(sol)
jacob_lohman_F_solution <- function(ptest, a, t0 ,t){
  if(class(ptest) != 'pumping_test'){
    stop('A pumping_test object is required as input')
  }
  td <- t/4.0/exp(digamma(1))/t0
  r <- ptest$r
  rw <- ptest$additional_parameters$rw
  rho <- r/rw
  par <- list(rho = rho, coeffs = ptest$coeffs)
  W <- jacob_lohman_F_well_function(td, par)
  s <- ptest$additional_parameters$s0*W
  pos <- s < 1e-12
  s[pos] <- 1e-12
  return(s)
}
#' @title
#' jacob_lohman_WF_F_LT
#' @description
#' Function to calculate the well function of the Jacob-Lohman solution in the
#' Laplace domain
#' @param p Laplace parameter
#' @param arg1 Value of the dimensionless radius (rd)
#' @param ... Additional parameters
#' @return
#' Value of the dimensionless drawdown in the Laplace domain of the Jacob-Lohman solution
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @export
#' @family jacob_lohman functions
#' @references
#' Jacob, C. E. & Lohman, S. W. Nonsteady Flow to a Well Of Constant Drawdown. American
#' Geophysical Union, 1952, 33, 10.
#' @examples
#' p <- 0.1
#' rd <- 5
#' sp <- jacob_lohman_WF_F_LT(p, rd)
jacob_lohman_WF_F_LT <- function(p, arg1,...){
  num <- besselK(2.0*sqrt(p)*arg1, 0)
  den <- besselK(2.0*sqrt(p), 0)
  return( (1.0/p)*(num/den) )
}
#' @title
#' jacob_lohman_G_well_function
#' @description
#' This calculates the well function of the Jacob-Lohman solution for a constant head aquifer
#' (artesian aquifer). This function calculates the dimensionless discharge.
#' @param td Numeric vector with the Dimensionless time
#' @param par A list with the values of the Stehfest coefficients.
#' @return
#' A numeric vector with the dimensionless drawdowns
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @export
#' @family jacob_lohman functions
#' @references
#' Jacob, C. E. & Lohman, S. W. Nonsteady Flow to a Well Of Constant Drawdown. American
#' Geophysical Union, 1952, 33, 10.
#' @examples
#' td <- logseq(-2,4,50)
#' par <- list(coeffs = stehfest_coeffs(8))
#' sd <- jacob_lohman_G_well_function(td, par)
#' print(sd)
jacob_lohman_G_well_function <- function(td, par){
  coeffs <- par$coeffs
  qd <- stehfest_inversion(td, coeffs, jacob_lohman_G_WF_LT)
  return(qd)
}
#' @title
#' jacob_lohman_G_solution_initial
#' @description
#' Calculates the value of the initial parameters of a constant head test using the
#' Jacob-Lohman solution (Discharge)
#' @param ptest A pumping_test object
#' @return
#' A list with the values of
#' \itemize{
#' \item a Slope of the straight line fitted to the drawdown data using the
#' Cooper-Jacob approach
#' \item t0 Intercept of the straight line fitted to the drawdown data using the
#' Cooper-Jacob approach
#' }
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @export
#' @family jacob_lohman functions
#' @references
#' Jacob, C. E. & Lohman, S. W. Nonsteady Flow to a Well Of Constant Drawdown. American
#' Geophysical Union, 1952, 33, 10.
#' @examples
#' data(jacob_lohman)
#' s <- vector("numeric", length(jacob_lohman$t))
#' s[1:length(jacob_lohman$t)] <- 0
#' ptest <- pumping_test("Well1", Q = jacob_lohman$q, r = 0.0, t = jacob_lohman$t, s = s)
#' sol0 <- jacob_lohman_G_solution_initial(ptest)
#' print(sol0)
jacob_lohman_G_solution_initial <- function(ptest){
  if(class(ptest) != 'pumping_test'){
    stop('A pumping_test object is required as input')
  }
  t <- ptest$t
  n <- ceiling(length(t)/3);
  Q <- ptest$Q
  nq <- length(Q)
  if(nq != length(t)){
    stop('Discharge values for each time are required')
  }
  endp <- length(t)
  t <- t[n:endp]
  s <- 1./Q[n:endp]
  jl.df <- data.frame(t,s)
  jl.lm <- lm(s ~ log10(t), data = jl.df)
  t1 <- 10^(-jl.lm$coefficients[1]/jl.lm$coefficients[2])
  parameters <- list(a = unname(jl.lm$coefficients[2]), t0 = ptest$t[1])
  res <- parameters
  return(res)
}
#' @title
#' jacob_lohman_G_calculate_parameters
#' @description
#' Function to calculate the hydraulic parameters for the Jacob-Lohman solution (Discharge)
#' @param ptest A pumping_test object
#' @param par A list with the slope and the intercept (a and t0) of the straight line fitted to
#' the drawdown data using the Cooper-Jacob approach
#' @return
#' A list with
#' \itemize{
#' \item Tr Transmissivity
#' \item Ss Storage coefficient
#' \item radius_influence
#' }
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @family jacob_lohman functions
#' @export
#' @references
#' Jacob, C. E. & Lohman, S. W. Nonsteady Flow to a Well Of Constant Drawdown. American
#' Geophysical Union, 1952, 33, 10.
#' @examples
#' data(jacob_lohman)
#' s <- vector("numeric", length(jacob_lohman$t))
#' s[1:length(jacob_lohman$t)] <- 0
#' ptest <- pumping_test("Well1", Q = jacob_lohman$q, r = 0.0, t = jacob_lohman$t, s = s)
#' sol0 <- jacob_lohman_G_solution_initial(ptest)
#' jl.par <- jacob_lohman_G_calculate_parameters(ptest, sol0)
#' print(jl.par)
jacob_lohman_G_calculate_parameters <- function(ptest, par){
  if(class(ptest) != 'pumping_test'){
    stop('A pumping_test object is required as input')
  }
  if(class(par) != 'list'){
    stop('A list is required as input')
  }
  endp <- length(t)
  t <- ptest$t
  Q <- ptest$Q
  a <- par$a
  t0 <- par$t0
  Q1 <- ptest$Q[1]
  r <- ptest$r
  #Compute the transmissivity, storativity and radius of influence
  Tr <- 0.1832339*Q1/a
  Ss <- 2.2458394*Tr*t0/r^2
  Ri <- 2.0*sqrt(Tr*t[endp]/Ss)
  res <- list(Tr = Tr, Ss = Ss, radius_influence = Ri)
  return(res)
}
#' @title
#' jacob_lohmna_G_solution
#' @description
#' Function to calculate the discharge at the pumping well using the Jacob-Lohman solution for a
#' constant head test (artesian aquifer)
#' @param ptest A pumping_test object
#' @param a Slope of the straight line fitted tothe drawdown data using the Cooper-Jacob approach
#' @param t0 Intercept of the straight line fitted tothe drawdown data using the Cooper-Jacob approach
#' @param t Numeric vector with the time values
#' @return
#' A numeric vector with the calculate discharge using the Jacob-Lohman solution
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @export
#' @family jacob_lohman functions
#' @references
#' Jacob, C. E. & Lohman, S. W. Nonsteady Flow to a Well Of Constant Drawdown. American
#' Geophysical Union, 1952, 33, 10.
#' @examples
#' data(jacob_lohman)
#' s <- vector("numeric", length(jacob_lohman$t))
#' s[1:length(jacob_lohman$t)] <- 0
#' ptest <- pumping_test("Well1", Q = jacob_lohman$q, r = 0.0, t = jacob_lohman$t, s = s)
#' sol0 <- jacob_lohman_G_solution_initial(ptest)
#' sol <- jacob_lohman_G_solution(ptest, sol0$a, sol0$t0, ptest$t)
#' print(sol)
jacob_lohman_G_solution <- function(ptest, a, t0, t){
  if(class(ptest) != 'pumping_test'){
    stop('A pumping_test object is required as input')
  }
  td <- t/4.0/exp(digamma(1))/t0
  par <- list(coeffs = ptest$coeffs)
  W <- jacob_lohman_G_well_function(td, par)
  q <- (log(10)/2)*W/a
  return(q)
}
#' @title
#' jacob_lohman_G_WF_LT
#' @description
#' Function to calculate the well function of the Jacob-Lohman solution (Discharge) in the
#' Laplace domain
#' @param p Laplace parameter
#' @param ... Additional parameters
#' @return
#' Value of the dimensionless discharge in the Laplace domain of the Jacob-Lohman solution
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @export
#' @family jacob_lohman functions
#' @references
#' Jacob, C. E. & Lohman, S. W. Nonsteady Flow to a Well Of Constant Drawdown. American
#' Geophysical Union, 1952, 33, 10.
#' @examples
#' p <- 0.1
#' qp <- jacob_lohman_G_WF_LT(p)
#' print(qp)
jacob_lohman_G_WF_LT <- function(p,...){
  sp <- sqrt(p)
  q  <- besselK(sp, 1)/(sp*besselK(sp, 0))
  return(q)
  #Asymptotic
  #return( 0.5/p/besselK(2.0*sqrt(p),0))
}