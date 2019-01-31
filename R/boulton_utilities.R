#' @section Boulton functions:
#' boulton_well_function, boulton_solution_initial, boulton_calculate_parameters, boulton_solution, boulton_solution_dlogt, boulton_WF_LT
#' @docType package
#' @name pumpingtest
NULL
#' @title
#' boulton_well_function
#' @description
#' Well function of the Boulton model for the phreatic aquifer
#' @param td Numeric dimensionless time
#' @param par list with the delay parameter (phi), the sigma parameter and the coeffs vector
#' @return
#' Value of the well function of the Boulton model
#' @family boulton functions
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @export
#' @references
#' Boulton, N. The drawdown of the water-table under non-steady conditions near a pumped
#' well in an unconfined formation. Proceedings of the Institution of Civil Engineers,
#' 1954, 3, 564-579
#' @examples
#' td <- logseq(-2, 4, 30)
#' par <- list(sigma =0.1, phi = 0.1)
#' W.boulton <- boulton_well_function(td, par)
#' plot(1/td, W.boulton, type = "p", log = "xy")
boulton_well_function <- function(td, par){
  sigma <- par$sigma
  phi1 <- par$phi
  res <- boulton_well_function_cpp(td, sigma, phi1)
  return(res)
}
#' @title
#' boulton_solution_initial
#' @description
#' Function to calculate the initial values of the hydraulic parameters of the Boulton model
#' @param ptest A pumping_test object
#' @return
#' A list with:
#' \itemize{
#' \item a: Slope of the straight line fitted to the drawdown data using the
#' Cooper-Jacob approach
#' \item t0: Intercept of the straight line fitted to the drawdown data using the
#' Cooper-Jacob approach (early time)
#' \item t1: Intercept of the straight line fitted to the drawdown data using the
#' Cooper-Jacob approach (late time)
#' \item phi: Delay parameter. Dimensionless parameter defined as
#' \deqn{\phi = \frac{\alpha_{1} r^{2} S}{T}}
#' where \eqn{\alpha_{1}} is a fitting parameter without a physical estimation, \eqn{r} is
#' the distance between the pumping and observation well, \eqn{S} is the storage coefficient and
#' \eqn{T} is the transmissivity.
#' }
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @export
#' @references
#' Boulton, N. The drawdown of the water-table under non-steady conditions near a pumped
#' well in an unconfined formation. Proceedings of the Institution of Civil Engineers,
#' 1954, 3, 564-579
#' @importFrom stats lm
#' @examples
#' data(boulton)
#' ptest <- pumping_test("Well1", Q = 0.03, r = 20, t = boulton$t, s = boulton$s)
#' boulton.sol <- boulton_solution_initial(ptest)
boulton_solution_initial <- function(ptest){
  if(class(ptest) != 'pumping_test'){
    stop('A pumping_test object is required as input')
  }
  n <- floor(length(ptest$t)/3)
  endp <- length(ptest$t)
  ptest1 <- as.data.frame(cbind(ptest$t[n:endp],ptest$s[n:endp]))
  names(ptest1) <- c('t', 's')
  initial <- lm(s ~ log10(t), data = ptest1)
  t1 <- 10^(-initial$coefficients[1]/initial$coefficients[2]);
  res <- list(a = unname(initial$coefficients[2]), t0 = ptest$t[1],
              t1 = unname(t1), phi = 1.0e-4)
  return(res)
}
#' @title
#' boulton_calculate_parameters
#' @description
#' Function to calculate the hydraulic parameters of the Boulton model
#' @param ptest A pumping_test object
#' @param par List with the parameters a (slope), t0 and t1 are the intercepts of the
#' straight line fitted to the drawdown using the Cooper-Jacob method for early and late
#' times respectively.
#' @param hydraulic Logical flag to indicate if hydraulic parameters are calculated. If False,
#' the the statistcal parameter are calculated.
#' @return
#' A list with the values of transmisivity (Tr), storage coefficient(Ss), Drainage Porosity
#' (omegad) and radius of influence (radius_influence)
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @family boulton functions
#' @export
#' @references
#' Boulton, N. The drawdown of the water-table under non-steady conditions near a pumped
#' well in an unconfined formation. Proceedings of the Institution of Civil Engineers,
#' 1954, 3, 564-579
#' @examples
#' data(boulton)
#' ptest <- pumping_test("Well1", Q = 0.03, r = 20, t = boulton$t, s = boulton$s)
#' boulton_sol0 <- boulton_solution_initial(ptest)
#' boulton.par <- boulton_calculate_parameters(ptest, boulton_sol0)
boulton_calculate_parameters <- function(ptest, par, hydraulic = TRUE){
  if(class(ptest) != 'pumping_test'){
    stop('A pumping_test object is required as input')
  }
  if(class(par) != 'list'){
    stop('A list is required as input')
  }
  Q <- ptest$Q
  r <- ptest$r
  t <- ptest$t
  if(hydraulic){
    a <- par$a
    t0 <- par$t0
    t1 <- par$t1
    pos_valid_a <- abs(a) < 1.0e-15
    pos_valid_t0 <- abs(t0) < 1.0e-15
    pos_valid_t1 <- abs(t1) < 1.0e-15
    if( sum(pos_valid_a) > 0 | sum(pos_valid_t0) > 0 | sum(pos_valid_t1) > 0){
      stop('a or t0 or t1 or all of them are close to numerical precision')
    }
    phi <- par$phi
    Tr <- 0.1832339*Q/a
    Ss <- 2.2458394*Tr*t0/r^2
    omegad <- 2.2458394*Tr*t1/r^2-Ss
    endp <- length(t)
    Ri <- 2*sqrt(Tr*t[endp]/omegad)
    res <- list(Tr = Tr, Ss = Ss, omegad = omegad,
                radius_influence = Ri)
  }
  else {
    Tr <- par$Tr
    Ss <- par$Ss
    omegad <- par$omegad
    a <- 0.1832339*Q/Tr
    t0 <- (Ss*r^2)/(2.2458394*Tr)
    t1 <- (r^2)*(omegad+Ss)/(2.2458394*Tr)
    res <- list(a = a, t0 = t0, t1 = t1, phi = par$phi)
  }
  return(res)
}
#' @title
#' boulton_solution
#' @description
#' Function to calculate the drawdown of a Boulton model
#' @param ptest A pumping_test object
#' @param a Slope of the straight line fitted to the drawdown data using the
#' Cooper-Jacob approach
#' @param t0 Intercept of the straight line fitted to the drawdown data using the
#' Cooper-Jacob approach (early time)
#' @param t1 Intercept of the straight line fitted to the drawdown data using the
#' Cooper-Jacob approach (late time)
#' @param phi Delay parameter. Dimensionless parameter defined as
#' \deqn{\phi = \frac{\alpha_{1} r^{2} S}{T}}
#' where \eqn{\alpha_{1}} is a fitting parameter without a physical estimation, \eqn{r} is
#' the distance between the pumping and observation well, \eqn{S} is the storage coefficient and
#' \eqn{T} is the transmissivity.
#' @param t Numeric vector with the time values
#' @return
#' A vector with the calculated drawdown
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @family boulton functions
#' @export
#' @references
#' Boulton, N. The drawdown of the water-table under non-steady conditions near a pumped
#' well in an unconfined formation. Proceedings of the Institution of Civil Engineers,
#' 1954, 3, 564-579
#' @examples
#' data(boulton)
#' ptest <- pumping_test("Well1", Q = 0.03, r = 20, t = boulton$t, s = boulton$s)
#' boulton_sol0 <- boulton_solution_initial(ptest)
#' boulton_sol1 <- boulton_solution(ptest, boulton_sol0$a, boulton_sol0$t0,
#'                 boulton_sol0$t1, boulton_sol0$phi, boulton$t)
boulton_solution <- function(ptest, a, t0, t1, phi, t){
  # Tr, Ss, Sy, a1
  # IMPORTANT:
  # Here 8 terms in the inversion are used, otherwise
  # it gives wrong results
  Q <- ptest$Q
  r <- ptest$r
  td <- 0.445268*t/t0;
  sigma <- t0/(t0+t1)
  phi1 <- 2*phi*t0
  sd1 <- 0.868589*a*stehfest_inversion(td, ptest$coeffs, boulton_WF_LT, sigma, phi1)
  pos <- is.na(sd1)
  if(sum(pos) > 0){
    sd1[pos] <- 0.0
  }
  return(sd1)
}
#' @title
#' boulton_WF_LT
#' @description
#' Function to calculate the Boulton solution in the Laplace domain.
#' @param p Laplace parameter
#' @param arg1 Value of sigma parameter defined as
#' \deqn{\sigma = \frac{S}{S_{y}} }
#' @param arg2 Value of the delay parameter(phi):
#' \deqn{\phi = \frac{\alpha_{1} r^{2} S}{T}}
#' @return
#' The value of the drawdown in the Laplace domain
#' @family boulton functions
#' @export
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @references
#' Boulton, N. The drawdown of the water-table under non-steady conditions near a pumped
#' well in an unconfined formation. Proceedings of the Institution of Civil Engineers,
#' 1954, 3, 564-579
#' @examples
#' sigma <- 0.1
#' phi <- 0.1
#' p <- 0.5
#' res <- boulton_WF_LT(p, sigma, phi)
boulton_WF_LT <- function(p, arg1, arg2){
  sigma <- abs(arg1)
  phi <- abs(arg2)
  p <- abs(p)
  sp <- sqrt( p + (phi*p)/(sigma*(p+phi)) )
  res <- besselK(sp, 0)/p
  return(res)
}
#' @title
#' boulton_solution_dlogt
#' @description
#' Function to calculate the derivative of the drawdown with respect to the
#' log of time using the Boulton solution
#' @param ptest value
#' @param a Slope of the straight line fitted to the drawdown data using the
#' Cooper-Jacob approach
#' @param t0 Intercept of the straight line fitted to the drawdown data using the
#' Cooper-Jacob approach (early time)
#' @param t1 Intercept of the straight line fitted to the drawdown data using the
#' Cooper-Jacob approach (late time)
#' @param phi Delay parameter. Dimensionless parameter defined as
#' \deqn{\phi = \frac{\alpha_{1} r^{2} S}{T}}
#' where \eqn{\alpha_{1}} is a fitting parameter without a physical estimation, \eqn{r} is
#' the distance between the pumping and observation well, \eqn{S} is the storage coefficient and
#' \eqn{T} is the transmissivity.
#' @param t value
#' @return
#' This function returns
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @export
#' @family boulton functions
#' @references
#' Boulton, N. The drawdown of the water-table under non-steady conditions near a pumped
#' well in an unconfined formation. Proceedings of the Institution of Civil Engineers,
#' 1954, 3, 564-579
#' @examples
#' data(boulton)
#' ptest <- pumping_test('Well1', Q = 0.03, r = 20, t = boulton$t, s = boulton$s)
#' boulton.sol <- boulton_solution_initial(ptest)
#' boulton.dsol <- boulton_solution_dlogt(ptest, boulton.sol$a, boulton.sol$t0,
#'                 boulton.sol$t1, boulton.sol$phi, boulton$t)
boulton_solution_dlogt <- function(ptest, a, t0, t1, phi, t){
  sb <- boulton_solution(ptest, a, t0, t1, phi, t)
  dl_sb <- log_derivative_central(t, sb)
  res <- dl_sb$y
  return(res)
}