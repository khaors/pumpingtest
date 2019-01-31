#' @section Papadopoulous-Cooper functions:
#' papadopulos_cooper_well_function, papadopulos_cooper_solution_initial, papadopulos_cooper_calculate_parameters, papadopulos_cooper_solution, papadopulos_cooper_solution_dlogt, papadopulos_cooper_WF_LT
#' @docType package
#' @name pumpingtest
NULL
#' @title
#' papadopulos_cooper_well_function
#' @description
#' Function to calculate the well function of the Papadopulos-Cooper solution that describe the drawdown
#' variations for large diameter wells
#' @param td Dimensionless time
#' @param par List with the  \eqn{cd} and \eqn{\rho} parameters.
#' The well storage (\eqn{cd} parameter) is defined as
#' \deqn{cd = \frac{r_{w}^{2} S}{r_{c}^{2}}}
#' where \eqn{r_{w}} is the radius of the well, \eqn{r_{c}} is the radius of the casing,
#' and \eqn{S} is the storage coefficient. The dimensionless radius (\eqn{\rho} parameter)
#'  is defined as:
#' \deqn{\rho=\frac{r}{r_{w}}}
#' where \eqn{r} is the distance to the observation well, and \eqn{r_{w}} is the radius
#' of the well.
#' @return
#' The value of the well function for the Papadopulos-Cooper model
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @family papadopulos_cooper functions
#' @export
#' @references
#' Papadopulos, I. S. & Cooper, H. H. Drawdown in a well of large diameter. Water
#' Resources Research, 1967, 3, 241-244.
#' @examples
#' td <- logseq(-2,4,30)
#' par <- list(cd = 0.1, rho = 0.1)
#' res <- papadopulos_cooper_well_function(td, par)
papadopulos_cooper_well_function <- function(td, par){
  cd <- par$cd
  rho <- par$rho
  #coeffs <- par$coeffs
  #res <- stehfest_inversion(td, coeffs, papadopulos_cooper_F_WF_LT, cd, rho)
  res <- papadopulos_cooper_well_function_cpp(td, cd, rho, 0.0)
  return(res)
}
#' @title
#' papadopulos_cooper_solution_initial
#' @description
#' Function to calculate the initial solution of the Papadopulos-Cooper model
#' @param ptest A pumping_test object
#' @return
#' A list with
#' \itemize{
#' \item a: Slope of the straight line fitted to the drawdown data using the
#' Cooper-Jacob approach
#' \item t0: Intercept of the straight line fitted to the drawdown data using the
#' Cooper-Jacob approach
#' }
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @family papadopulos_cooper functions
#' @export
#' @references
#' Papadopulos, I. S. & Cooper, H. H. Drawdown in a well of large diameter. Water
#' Resources Research, 1967, 3, 241-244.
#' @examples
#' data(papadopulos_cooper)
#' ptest <- pumping_test("Well1", Q=0.0050472, r = 3.048, t = as.numeric(papadopulos_cooper$t),
#'                       s = papadopulos_cooper$s)
#' sol0 <- papadopulos_cooper_solution_initial(ptest)
#' print(sol0)
papadopulos_cooper_solution_initial <- function(ptest){
  if(class(ptest) != 'pumping_test'){
    stop('A pumping_test object is required as input')
  }
  ptest_dlogt <- log_derivative_spline(ptest$t, ptest$s)
  endp <- length(ptest_dlogt$y)
  a <- log(10)*ptest_dlogt$y[endp]
  t0 <- ptest$t[endp]*exp(-ptest$s[endp]/ptest_dlogt$y[endp])
  res <- list(a = a, t0 = t0)
  return(res)
}
#' @title
#' papadopulos_cooper_calculate_parameters
#' @description
#' Function to calculate the hydraulic parameters of the Papadopulos-Cooper solution
#' @param ptest A pumping test object
#' @param par A list with the values of a (slope) and t0 (Intercept) of the straight line fitted
#'  to the drawdown data using the Cooper-Jacob approach.
#' @param hydraulic Logical flag to indicate if hydraulic parameters are calculated. If False,
#' the the statistcal parameters (a and t0) are calculated.
#' @return
#' A list with the hydraulic parameters:
#' \itemize{
#' \item Tr: Transmissivity
#' \item Ss: Storage coefficient
#' \item cd: Wellbore Storage
#' \item rho: Dimensionless radius
#' }
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @family papadopulos_cooper functions
#' @export
#' @references
#' Papadopulos, I. S. & Cooper, H. H. Drawdown in a well of large diameter. Water
#' Resources Research, 1967, 3, 241-244.
#' @examples
#' data(papadopulos_cooper)
#' ptest <- pumping_test("Well1", Q=0.0050472, r = 3.048, t = as.numeric(papadopulos_cooper$t),
#'                       s = papadopulos_cooper$s)
#' par <- list(rw = 0.6096, rc = 0.6096 )
#' ptest$additional_parameters <- par
#' sol0 <- papadopulos_cooper_solution_initial(ptest)
#' pc.par <- papadopulos_cooper_calculate_parameters(ptest, sol0)
#' print(pc.par)
papadopulos_cooper_calculate_parameters <- function(ptest, par, hydraulic = TRUE){
  if(class(ptest) != 'pumping_test'){
    stop('A pumping_test object is required as input')
  }
  if(class(par) != 'list'){
    stop('A list is required as input')
  }
  Q <- ptest$Q
  r <- ptest$r
  if(is.null(ptest$additional_parameters)){
    stop('The additional_parameters list is empty')
  }
  rw <- ptest$additional_parameters$rw
  rc <- ptest$additional_parameters$rc
  if(hydraulic){
    a <- par$a
    t0 <- par$t0
    pos_valid_a <- abs(a) < 1.0e-15
    pos_valid_t0 <- abs(t0) < 1.0e-15
    if( sum(pos_valid_a) > 0 | sum(pos_valid_t0) > 0){
      stop('a or t0 or both are close to numerical precision')
    }
    Tr <- 0.1832339*Q/a     # Compute the transmissivity
    Ss <- 2.25*Tr*t0/r^2    # Compute the storativity
    alpha <- (rc^2)/(2*rw^2*Ss)  # Compute the wellbore storage coefficient Cd
    rho <- r/rw
    res <- list(Tr = Tr, Ss = Ss, cd = alpha)
  }
  else {
    Tr <- par$Tr
    Ss <- par$Ss
    cd <- par$cd
    rho <- par$rho
    a <- 0.1832339*Q/Tr
    t0 <- (Ss*r^2)/(2.25*Tr)
    res <- list(a = a, t0 = t0)
  }
  return(res)
}
#' @title
#' papadopulos_cooper_solution
#' @description
#' Function to calculate the drawdown of the Papadopulos-Cooper solution
#' @param ptest A pumping_test object
#' @param a Slope of the straight line fitted to the drawdown data using the
#' Cooper-Jacob approach
#' @param t0 Intercept of the straight line fitted to the drawdown data using the
#' Cooper-Jacob approach
#' @param t Numeric vector with time
#' @return
#' A numeric vector with the calculated drawdown
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @export
#' @family papadopulos_cooper functions
#' @references
#' Papadopulos, I. S. & Cooper, H. H. Drawdown in a well of large diameter. Water
#' Resources Research, 1967, 3, 241-244.
#' @examples
#' data(papadopulos_cooper)
#' ptest <- pumping_test("Well1", Q=0.0050472, r = 3.048, t = as.numeric(papadopulos_cooper$t),
#'                       s = papadopulos_cooper$s)
#' ptest$additional_parameters <- list(rw = 0.6096, rc = 0.6096)
#' sol0 <- papadopulos_cooper_solution_initial(ptest)
#' sol <- papadopulos_cooper_solution(ptest, sol0$a, sol0$t0, ptest$t)
#' print(sol)
papadopulos_cooper_solution <- function(ptest, a, t0, t){
  if(class(ptest) != 'pumping_test'){
    stop('A pumping_test object is required as input')
  }
  par <- list(a = a, t0 =t0)
  hydr_par <- papadopulos_cooper_calculate_parameters(ptest, par)
  Tr <- hydr_par$Tr
  Ss <- hydr_par$Ss
  rw <- ptest$additional_parameters$rw
  rho <- ptest$r/rw
  alpha <- hydr_par$cd
  td <- 0.445268*(t/t0)*(rho^2)
  par1 <- list(cd = alpha, rho = rho,
               coeffs = ptest$coeffs)
  sd <- papadopulos_cooper_well_function(td, par1)
  s <- 0.868589*a*sd
  return(s)
}
#' @title
#' papadopulos_cooper_solution_dlogt
#' @description
#' This function calcultes the derivative of the drawdown with respect to the log of time for
#' Papadopulos-Cooper solution
#' @param ptest A pumping_test object
#' @param a Slope of the straight line fitted to the drawdown data using the
#' Cooper-Jacob approach
#' @param t0 Intercept of the straight line fitted to the drawdown data using the
#' Cooper-Jacob approach
#' @param t Numeric vector with time
#' @return
#' A list with the x and y coordinates of the derivative
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @family papadopulos_cooper functions
#' @export
#' @references
#' Papadopulos, I. S. & Cooper, H. H. Drawdown in a well of large diameter. Water
#' Resources Research, 1967, 3, 241-244.
#' @examples
#' data(papadopulos_cooper)
#' ptest <- pumping_test("Well1", Q=0.0050472, r = 3.048, t = as.numeric(papadopulos_cooper$t),
#'                       s = papadopulos_cooper$s)
#' ptest$additional_parameters <- list(rw = 0.6096, rc = 0.6096)
#' sol0 <- papadopulos_cooper_solution_initial(ptest)
#' sol_dlogt <- papadopulos_cooper_solution_dlogt(ptest, sol0$a, sol0$t0, ptest$t)
#' print(sol_dlogt)
papadopulos_cooper_solution_dlogt <- function(ptest, a, t0, t){
  sb <- papadopulos_cooper_solution(ptest, a, t0, t)
  dl_sb <- log_derivative_central(t, sb)
  res <- dl_sb$y
  return(res)
}
#' @title
#' papadopulos_cooper_F_WF_LT
#' @description
#' Function to calculate the drawdown of the Papadopulos-Cooper solution in the Laplace domain
#' @param p Laplace parameter
#' @param arg1 Dimensionless wellbore storage coefficient (cd)
#' @param arg2 Dimensionless radius (rho)
#' @return
#' A numeric vector with the drawdown in the Laplace domain
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @family papadopulos_cooper functions
#' @export
#' @references
#' Papadopulos, I. S. & Cooper, H. H. Drawdown in a well of large diameter. Water
#' Resources Research, 1967, 3, 241-244.
#' @examples
#' p <- 0.5
#' cd <- 1e-2
#' rho <- 0.5
#' res <- papadopulos_cooper_F_WF_LT(p, cd, rho)
#' print(res)
papadopulos_cooper_F_WF_LT <- function(p, arg1, arg2) {
  num <- besselK(2.0*sqrt(p),0)
  var1 <- (sqrt(p)/arg2)
  den1 <- besselK(var1,1)
  var2 <- p/(arg1*arg2^2)
  den2 <- besselK(var1,0)
  den <- p*(var1*den1+var2*den2)
  return( num/den )
}
#' @title
#' papadopulos_cooper_F_WF_LT_dlogt
#' @description
#' Calculates the derivative of the drawdown with respect to the logarithm of time for the
#' papadopulos_cooper solution in the Laplace domain
#' @param p Laplace parameter
#' @param arg1 Value of dimensionless wellbore storage coefficient (cd)
#' @param arg2 Value of dimensionless radius (rho)
#' @return
#' A numeric vector with the derivative of the drawdown with respect to the logarithm
#' of time for thepapadopulos_cooper solution in the Laplace domain
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @family papadopulos_cooper functions
#' @export
#' @references
#' Papadopulos, I. S. & Cooper, H. H. Drawdown in a well of large diameter. Water
#' Resources Research, 1967, 3, 241-244.
#' @examples
#' data(papadopulos_cooper)
papadopulos_cooper_F_WF_LT_dlogt <- function(p, arg1, arg2){
  u <- papadopulos_cooper_F_WF_LT(p, arg1, arg2)
  cd <- arg1
  rho <- arg2
  alpha <- 1.0
  sp <- sqrt(p)
  k1 <- besselK(2.0*sp, 1)
  k0 <- besselK(2.0*sp/rho, 0)
  k2 <- besselK(2.0*sp/rho, 2)
  term1_den <- (p*k0/rho^2) + (p^(3/2))*k1/rho
  term1 <- k1/(sp*term1_den)
  term2 <- k0/(rho^2) - sp*k1/(rho^3) + 3*sp*k1/(2.0*rho) +
    p*(-k0-k2)/(2.0*rho^2)
  term2_den <- ((p*k0 /(alpha*rho^2) + (p^(3/2)*k1/rho)))^2
  term2 <- (besselK(2.0*sp, 0)*term2)/term2_den
  dudp <- term1 + term2
  res <- -u+p*dudp
  return(res)
}
