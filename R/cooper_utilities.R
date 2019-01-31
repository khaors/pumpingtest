#' @section Cooper functions:
#' cooper_well_function, cooper_solution_initial, cooper_calculate_parameters, cooper_solution, cooper_solution_dlogt, cooper_WF_LT
#' @docType package
#' @name pumpingtest
NULL
#' @title
#' cooper_well_function
#' @description
#' Function to calculate the well function of the cooper solution for slug test
#' @param td Numeric vector with dimensionless time
#' @param par List with the stehfest coefficients and cd (dimensionless storage coefficient)
#' @return
#' A numeric vector with the value of the Cooper well function, that is, the dimensionless drawdown
#' @family cooper functions
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @export
#' @references
#' Cooper, H. H.; Bredehoeft, J. D. & Papadopulos, I. S. Response of a finite-diameter
#' well to an instantaneous charge of water Water Resources Research, 1967, 3, 263-269
#' @examples
#' td <- logseq(-2, 4, 50)
#' par <- list(cd = 1, rw = 1)
#' res <- cooper_well_function(td, par)
cooper_well_function <- function(td, par){
  cd <- par$cd
  #coeffs <- par$coeffs
  #res <- stehfest_inversion(td, coeffs, cooper_WF_LT, cd)
  res <- cooper_well_function_cpp(td, cd, 0.0, 0.0)
  return(res)
}
#' @title
#' cooper_solution_initial
#' @description
#' Function to calculate the initial values of the Cooper solution
#' @param ptest A pumping_test object
#' @return
#' A list with
#' \itemize{
#'   \item cd: Initial value of the wellbore storage coefficient
#'   \item t0: Intercept of the straight line fitted to the drawdown data using the
#' Cooper-Jacob approach
#' }
#' @family cooper functions
#' @importFrom pracma interp1
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @export
#' @references
#' Cooper, H. H.; Bredehoeft, J. D. & Papadopulos, I. S. Response of a finite-diameter
#' well to an instantaneous charge of water Water Resources Research, 1967, 3, 263-269
#' @examples
#' data(cooper_slug)
#' ptest <- pumping_test("Well1", Q = 0.0, r = 0.0, t = cooper_slug$t, s = cooper_slug$s)
#' sol0 <- cooper_solution_initial(ptest)
#' print(sol0)
cooper_solution_initial <- function(ptest){
  if(class(ptest) != 'pumping_test'){
    stop('A pumping_test object is required as input')
  }
  t <- ptest$t
  s <- ptest$s
  dsdlnt <- log_derivative_central(t, -s)
  dmax <- max(dsdlnt$y)
  pdmax <- which.max(dsdlnt$y)
  #print(c(dmax,pdmax))
  #load('R/sysdata.rda')
  cd0 <- exp(interp1(cooper$dsmax, log(cooper$cdmax), dmax))
  res.sort <- sort(s, index.return = TRUE)
  s.sort <- res.sort$x
  s.ix <- res.sort$ix
  t.sort <- t[s.ix]
  #pos_valid <- !duplicated(t.sort)
  #t.sort <- t.sort[pos_valid]
  pos_valid <- !duplicated(s.sort)
  s.sort <- s.sort[pos_valid]
  t.sort <- t.sort[pos_valid]
  t50 <- interp1(s.sort, t.sort, 0.5)
  t0 <- t50/exp(pracma::interp1(log(cooper$cdmax), log(td05$td05), log(cd0)))
  res <- list(cd = cd0, t0 = t0)
  return(res)
}
#' @title
#' cooper_calculate_parameters
#' @description
#' Function to calculate the hydraulic parameters of a cooper solution for a slug test
#' @param ptest A pumping_test object
#' @param par A list with the dimensionless wellbore storage coefficient (cd) and the
#' Intercept of the straight line fitted to the drawdown data using the
#' Cooper-Jacob approach (t0)
#' @param hydraulic Logical flag to indicate if hydraulic parameters are calculated. If
#' False, the the statistcal parameter are calculated.
#' @return
#' A list with
#' \itemize{
#'   \item Tr: Transmisivity
#'   \item Ss: Storage coefficient
#'   \item radius_influence: radius of influence
#' }
#' @family cooper functions
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @export
#' @references
#' Cooper, H. H.; Bredehoeft, J. D. & Papadopulos, I. S. Response of a finite-diameter
#' well to an instantaneous charge of water Water Resources Research, 1967, 3, 263-269
#' @examples
#' data(cooper_slug)
#' ptest <- pumping_test("Well1", Q = 0, r = 0, t = cooper_slug$t, s = cooper_slug$s)
#' par <- list(rw = 0.071, rc = 0.025)
#' ptest$additional_parameters <- par
#' sol0 <- cooper_solution_initial(ptest)
#' hydr.par <- cooper_calculate_parameters(ptest,sol0)
#' print(hydr.par)
cooper_calculate_parameters <- function(ptest, par, hydraulic = TRUE){
  if(class(ptest) != 'pumping_test'){
    stop('A pumping_test object is required as input')
  }
  t <- ptest$t
  s <- ptest$s
  if(is.null(ptest$additional_parameters)){
    stop('The additional_parameters list is empty')
  }
  rc <- ptest$additional_parameters$rc
  rw <- ptest$additional_parameters$rw
  if(hydraulic){
    endpos <- length(t)
    cd <- par$cd
    t0 <- par$t0
    pos_valid_cd <- abs(cd) < 1.0e-15
    pos_valid_t0 <- abs(t0) < 1.0e-15
    if( sum(pos_valid_cd) > 0 | sum(pos_valid_t0) > 0){
      stop('cd or t0 or both are close to numerical precision')
    }
    Tr <- 0.5*rc^2/t0;
    Ss <- 0.5/cd*(rc/rw)^2
    #
    n <- 0.462
    m <- 0.495
    xlim <- 2.5
    tdl <- cd*t[endpos]/t0
    x <- tdl^n/cd^m
    if(x < xlim){
      radius_influence <- rw*3.54*tdl^n
    }
    else {
      radius_influence <- rw*8.37*cd^m
    }
    res <- list(Tr = Tr, Ss = Ss, radius_influence = radius_influence)
  }
  else {
    Tr <- par$Tr
    Ss <- par$Ss
    t0 <- (0.5*rc^2)/Tr
    cd <- 0.5/Ss*(rc/rw)^2
    res <- list(cd = cd, t0 = t0)
  }
  return(res)
}
#' @title
#' cooper_solution
#' @description
#' Function to calculate the drawdown using the Cooper solution for a slug test
#' @param ptest A pumping_test object
#' @param cd Value of the dimensionless wellbore storage coefficient
#' @param t0 Intercept of the straight line fitted to the drawdown data using the
#' Cooper-Jacob approach
#' @param t Numeric vector with the value of time
#' @return
#' A numeric vector with the calculated drawdown
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @family cooper functions
#' @export
#' @references
#' Cooper, H. H.; Bredehoeft, J. D. & Papadopulos, I. S. Response of a finite-diameter
#' well to an instantaneous charge of water Water Resources Research, 1967, 3, 263-269
#' @examples
#' data(cooper_slug)
#' ptest <- pumping_test("Well1", Q = 0, r = 0, t = cooper_slug$t, s = cooper_slug$s)
#' sol0 <- cooper_solution_initial(ptest)
#' sol <- cooper_solution(ptest, sol0$cd, sol0$t0, ptest$t)
#' print(sol)
cooper_solution <- function(ptest, cd, t0 , t){
  if(class(ptest) != 'pumping_test'){
    stop('A pumping_test object is required as input')
  }
  #par <- list(cd = cd, t0 = t0)
  #hydr_par <- papadopulos_cooper_calculate_parameters(ptest, par)
  #Tr <- hydr_par$Tr
  #Ss <- hydr_par$Ss
  td <- t/t0
  par <- list(cd = cd, coeffs = ptest$coeffs)
  s <- cooper_well_function(td, par)
  return(s)
}
#' @title
#' cooper_WF_LT
#' @description
#' Cooper's well function for slug tests in the Laplace domain
#' @param p Laplace parameter
#' @param arg1 Value of the dimensionless wellbore storage coefficient (cd)
#' @param arg2 Value or the dimensionless radius (rd)
#' @return
#' The value of the Cooper's solution for slug test in the Laplace domain
#' @family cooper functions
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @export
#' @references
#' Cooper, H. H.; Bredehoeft, J. D. & Papadopulos, I. S. Response of a finite-diameter
#' well to an instantaneous charge of water Water Resources Research, 1967, 3, 263-269
#' @examples
#' p <- seq(1,10,0.1)
#' res <- cooper_WF_LT(p, 0.1, 0.1)
#' print(res)
cooper_WF_LT <- function(p, arg1, arg2){
  cd <- arg1
  rd <- arg2
  if(p < 0){
    p <-1e-3
  }
  sp <- sqrt(p)
  if(is.na(sp)){
    print(p)
  }
  res_num <- cd * besselK(sp, 0)
  res_den <- (cd*p*besselK(sp, 0) + sp*besselK(sp,1))
  return(res_num/res_den)
}
#' @title
#'cooper_solution_dlogt
#' @description
#' Function to calculate the derivative of the drawdown with respect to the log of time for the
#' Cooper solution for slug tests
#' @param ptest A pumping_test object
#' @param cd Value of the dimensionless storage coefficient (cd)
#' @param t0 Intercept of the straight line fitted to the drawdown data using the
#' Cooper-Jacob approach
#' @param t A numeric vector with the time values
#' @return
#' A numeric vector with the values of the drawdown
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @family cooper functions
#' @export
#' @references
#' Cooper, H. H.; Bredehoeft, J. D. & Papadopulos, I. S. Response of a finite-diameter
#' well to an instantaneous charge of water Water Resources Research, 1967, 3, 263-269
#' @examples
#' data(cooper_slug)
#' ptest <- pumping_test("Well1", Q = 0, r = 0, t = cooper_slug$t, s = cooper_slug$s)
#' sol0 <- cooper_solution_initial(ptest)
#' sol_dlogt <- cooper_solution_dlogt(ptest, sol0$cd, sol0$t0, ptest$t)
#' print(sol_dlogt)
cooper_solution_dlogt <- function(ptest, cd, t0, t){
  sb <- cooper_solution(ptest, cd, t0, t)
  dl_sb <- log_derivative_central(t, -sb)
  res <- dl_sb$y
  return(res)
}