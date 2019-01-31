#' @section Agarwal functions (Pumping tests with skin effect):
#'
#' The functions included in this section are:
#'
#' agarwal_skin_well_function, agarwal_skin_solution_initial, agarwal_skin_calculate_parameters, agarwal_skin_solution, agarwal_skin_WF_LT, agarwal_skin_solution_dlogt
#'
#' @docType package
#' @name pumpingtest
NULL
#' @title
#' agarwal_skin_well_function
#' @description
#' Function to calculate the well function using the Agarwal solution (skin and wellbore storage)
#' @param td Numeric vector with time
#' @param par A list with the values of wellbore storage(cd), dimensionless radius (rd)
#'  and skin factor (sigma)
#' @return
#' Numeric vector with the values of dimensionless drawdown
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @export
#' @family agarwal_skin functions
#' @references
#' Agarwal, R.; Al-Hussainy, R. & Ramey, H. An Investigation of Wellbore Storage and
#' Skin Effect in Unsteady Liquid Flow: I. Analytical Treatment. Society of Petroleum
#' Engineers Journal., 1970
#' @examples
#' data(papadopulos_cooper)
agarwal_skin_well_function <- function(td, par){
  coeffs <- par$coeffs
  arg1 <- par$cd
  arg2 <- par$rd
  arg3 <- par$sigma
  sd <- stehfest_inversion(td, coeffs, agarwal_skin_WF_LT, arg1 = arg1,
                           arg2 = arg2, arg3 = arg3)
  return(sd)
}
#' @title
#' agarwal_skin_solution_initial
#' @description
#' Function to calculate the initial solution of the Agarwal solution (skin and wellbore storage)
#' @param ptest A pumping_test object
#' @return
#' A list with :
#' \itemize{
#' \item a: Slope of the straight line fitted to the drawdown data using the
#' Cooper-Jacob approach
#' \item t0: Intercept of the straight line fitted to the drawdown data using the
#' Cooper-Jacob approach
#' \item sigma: skin factor
#' }
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @export
#' @family agarwal_skin functions
#' @references
#' Agarwal, R.; Al-Hussainy, R. & Ramey, H. An Investigation of Wellbore Storage and
#' Skin Effect in Unsteady Liquid Flow: I. Analytical Treatment. Society of Petroleum
#' Engineers Journal., 1970
#' @examples
#' data(papadopulos_cooper)
agarwal_skin_solution_initial <- function(ptest){
  if(class(ptest) != 'pumping_test'){
    stop('A pumping_test object is required as input')
  }
  #initial <- cooper_jacob_solution_initial(ptest)
  #return(initial)
  rc <- ptest$additional_parameters$rc
  rw <- ptest$additional_parameters$rw
  id <- floor(2.0*length(ptest$t)/3.0)
  #print(id)
  endp <- length(ptest$t)
  ptest1 <- pumping_test(id=ptest$id, Q = ptest$Q, r = ptest$r,
                         t = ptest$t[id:endp], s = ptest$s[id:endp])
  sol0 <- theis_solution_initial(ptest1)
  #print(sol0)
  dptest_logt <- NULL
  dptest_logt <- log_derivative(ptest$t, ptest$s, method = "spline", d = 10)
  endp <- length(dptest_logt$y)
  endp1 <- length(ptest$t)
  a <- log(10)*dptest_logt$y[endp]
  #print( c(ptest$t[endp1],-ptest$s[endp1],dptest_logt$y[endp]))
  t0 <- ptest$t[endp1]*exp(-ptest$s[endp1]/dptest_logt$y[endp])
  res <- list(a = a, t0 =t0, sigma = 1.0)
  return(res)
}
#' @title
#' agarwal_skin_calculate_parameters
#' @description
#' Function to calculate the hydraulic parameters for Agarwal solution (skin and wellbore storage)
#' @param ptest A pumping_test object
#' @param par A list with the values of slope and intercept (a and t0) of the straight
#' line fitted to the drawdown data using the Cooper-Jacob approach, and the value of
#' skin factor(sigma)
#' @param hydraulic Logical flag to indicate if hydraulic parameters are calculated. If False,
#' the the statistcal parameter (a and t0) are calculated.
#' @return
#' A list with
#' \itemize{
#' \item Tr: Transmissivity
#' \item Ss: Storage coefficient
#' \item cd: Wellbore storage
#' \item rd: Dimensionless radius
#' }
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@mgmail.com}
#' @export
#' @family agarwal_skin functions
#' @references
#' Agarwal, R.; Al-Hussainy, R. & Ramey, H. An Investigation of Wellbore Storage and
#' Skin Effect in Unsteady Liquid Flow: I. Analytical Treatment. Society of Petroleum
#' Engineers Journal., 1970
#' @examples
#' data(papadopulos_cooper)
agarwal_skin_calculate_parameters <- function(ptest, par, hydraulic = TRUE){
  if(class(ptest) != 'pumping_test'){
    stop('A pumping_test object is required as input')
  }
  Q <- ptest$Q
  r <- ptest$r
  rw <- ptest$additional_parameters$rw
  rc <- ptest$additional_parameters$rc
  if(hydraulic){
    a <- par$a
    t0 <- par$t0
    sigma <- par$sigma
    pos_valid_a <- abs(a) < 1.0e-15
    pos_valid_t0 <- abs(t0) < 1.0e-15
    pos_valid_sigma <- abs(sigma) < 1.0e-15
    if( sum(pos_valid_a) > 0 | sum(pos_valid_t0) > 0 | sum(pos_valid_sigma) > 0){
      stop('a or t0 or sigma or all of them are close to numerical precision')
    }
    # Transmissivity
    Tr <- 0.1832339*Q/a
    #Storativity
    Ss <- 2.25*Tr*t0/(r^2)
    #Wellbore coefficient
    cd <- (rc^2)/2/rw^2/Ss
    #Dimensionless radius rd
    rd <- r/rw
    #
    res <- list(Tr = Tr, Ss = Ss, cd = cd, rd = rd)
  }
  else {
    rw <- ptest$additional_parameters$rw
    Tr <- par$Tr
    Ss <- par$Ss
    cd <- par$cd
    rd <- par$rd
    a <- 0.1832339*Q/Tr
    t0 <- (Ss*r^2)/(2.25*Tr)
    sigma <- 0.5*log((Ss*rw^2)/(4*exp(-digamma(1))*Tr*t0))
    res <- list(a = a, t0 = t0, sigma = sigma)
  }
  return(res)
}
#' @title
#' agarwal_skin_solution
#' @description
#' Function to calculate the agarwal solution (skin and wellbore storage)
#' @param ptest A pumping_test object
#' @param a Slope of the straight line fitted to the drawdown data using the
#' Cooper-Jacob approach
#' @param t0 Intercept of the straight line fitted to the drawdown data using the
#' Cooper-Jacob approach
#' @param sigma The value of the skin factor
#' @param t Numeric vector with time
#' @return
#' A numeric vector with calculated drawdown
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @export
#' @family agarwal_skin functions
#' @references
#' Agarwal, R.; Al-Hussainy, R. & Ramey, H. An Investigation of Wellbore Storage and
#' Skin Effect in Unsteady Liquid Flow: I. Analytical Treatment. Society of Petroleum
#' Engineers Journal., 1970
#' @examples
#' data(papadopulos_cooper)
agarwal_skin_solution <- function(ptest, a, t0, sigma, t){
  if(class(ptest) != 'pumping_test'){
    stop('A pumping_test object is required as input')
  }
  Q <- ptest$Q
  r <- ptest$r
  rc <- ptest$additional_parameters$rc
  rw <- ptest$additional_parameters$rw
  par <- list(a = a, t0 = t0, sigma = sigma)
  hydr.par <- agarwal_skin_calculate_parameters(ptest, par)
  Tr <- hydr.par$Tr
  Ss <- hydr.par$Ss
  cd <- (rc^2)/(2*(rw^2)*Ss)
  rd <- r/rw
  td <- (0.445268*t)/(t0*(rd^2))
  par <- list(cd = cd, rd = rd, coeffs = ptest$coeffs, sigma = sigma)
  W <- agarwal_skin_well_function(td, par)
  s <- (2.0/log(10))*a*W
  pos <- s < 0
  s[pos] <- 0
  return(s)
}
#' @title
#' agarwal_skin_WF_LT
#' @description
#' Function to calculate the well function of the Agarwal solution (skin and wellbore storage) in the Laplace domain
#' @param p Laplace Transform parameter
#' @param arg1 value of wellbore coefficient (cd)
#' @param arg2 value of dimensionless radius (rd)
#' @param arg3 value of skin factor (sigma)
#' @return
#' A numeric vector with the calculated drawdown in the Laplace domain
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @export
#' @family agarwal_skin functions
#' @references
#' Agarwal, R.; Al-Hussainy, R. & Ramey, H. An Investigation of Wellbore Storage and
#' Skin Effect in Unsteady Liquid Flow: I. Analytical Treatment. Society of Petroleum
#' Engineers Journal., 1970
#' @examples
#' data(papadopulos_cooper)
agarwal_skin_WF_LT <- function(p, arg1, arg2, arg3){
  cd <- arg1
  rd <- arg2
  sigma <- arg3
  sp <- sqrt(p)
  k0 <- besselK(sp, 0)
  k1 <- besselK(sp, 1)
  s <- besselK(rd*sp,0)/(p*(((1.0+p*cd*sigma)*sp*k1)+(cd*p*k0)))
  return(s)
}
#' @title
#' agarwal_skin_solution_dlogt
#' @description
#' Function to calculate the derivative of drawdown with respect to log of time of the Agarwal solution (skin and wellbore storage)
#' @param ptest A pumping_test object
#' @param a Slope of the straight line fitted to the drawdown data using the
#' Cooper-Jacob approach
#' @param t0 Intercept of the straight line fitted to the drawdown data using the
#' Cooper-Jacob approach
#' @param sigma Value of the skin factor
#' @param t Numeric vector with time
#' @return
#' A list with
#' \itemize{
#' \item x: time coordinates
#' \item y: Derivative value
#' }
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @export
#' @family agarwal_skin functions
#' @references
#' Agarwal, R.; Al-Hussainy, R. & Ramey, H. An Investigation of Wellbore Storage and
#' Skin Effect in Unsteady Liquid Flow: I. Analytical Treatment. Society of Petroleum
#' Engineers Journal., 1970
#' @examples
#' data(papadopulos_cooper)
agarwal_skin_solution_dlogt <- function(ptest, a, t0, sigma, t){
  if(class(ptest) != 'pumping_test'){
    stop('A pumping_test object is required as input')
  }
  sb <- agarwal_skin_solution(ptest, a, t0, sigma, t)
  dl_sb <- log_derivative_central(t, sb)
  res <- dl_sb$y
  return(res)
}