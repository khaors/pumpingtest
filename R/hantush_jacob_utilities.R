#' @section Hantush-Jacob functions:
#' hantush_jacob_well_function, hantush_jacob_solution_initial, hantush_jacob_calculate_parameters, hantush_jacob_solution, hantush_jacob_solution_dlogt, hantush_jacob_WF_LT
#' @docType package
#' @name pumpingtest
NULL
#' @title
#' hantush_jacob_well_function
#' @description
#' Well function of the Hantush-Jacob model
#' @param td Dimensionless time
#' @param par List with the beta parameter and the coefficients used in the Stehfest
#' inversion approach. The beta parameter is defined as
#' \deqn{\beta = \frac{r}{4} \sqrt{\frac{k_{h}^{'}S_{s}^{'}}{TS}}}
#' where \eqn{T} and \eqn{S} are the transmissivity and storage coefficient of the
#' aquifer, and \eqn{k_{h}^{'}} and \eqn{S_{s}^{'}} are the vertical hydraulic
#' conductivity and the Storage coefficient of the aquitard respectively.
#' @return
#' The value of the dimensionless drawdown of the Hantush-Jacob model
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @export
#' @family hantush_jacob functions
#' @references
#' Hantush, M. S. Flow to wells in aquifers separated by a semipervious layer. Journal
#' of Geophysical Research, 1967, 72, 1709.
#' @examples
#' td <- logseq(-2, 4, 50)
#' par <- list( beta = 0.5)
#' res <- hantush_jacob_well_function(td, par)
#' print(res)
hantush_jacob_well_function <- function(td, par){
  beta <- par$beta
  res <- hantush_jacob_well_function_cpp(td, beta, 0.0, 0.0)
  return(res)
}
#' @title
#' hantush_jacob_solutionb_initial
#' @description
#' Function to calculate the initial parameters of the Hantush-Jacob model
#' @param ptest A pumping_test object
#' @return
#' A list with
#' \itemize{
#' \item a: Slope of the straight line fitted to the drawdown data using the
#' Cooper-Jacob approach
#' \item t0: Intercept of the straight line fitted to the drawdown data using the
#' Cooper-Jacob approach
#' \item lambda: Leakage factor defined as
#' \deqn{\lambda = \sqrt{\frac{b^{'}T}{k_{h}^{'}}} }
#' }
#' where \eqn{b^{'}} and \eqn{k_{h}^{'}} are the thickness and the vertical hydraulic
#' conductivity of the aquitard, and \eqn{T} is the transmissivity of the underlying
#' aquifer.
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @export
#' @family hantush_jacob functions
#' @references
#' Hantush, M. S. Flow to wells in aquifers separated by a semipervious layer. Journal
#' of Geophysical Research, 1967, 72, 1709.
#' @examples
#' data(hantush_jacob)
#' ptest <- pumping_test('Well1', Q = 6.309e-3, r = 3.048, t=hantush_jacob$t, s = hantush_jacob$s)
#' sol0 <- hantush_jacob_solution_initial(ptest)
#' print(sol0)
hantush_jacob_solution_initial <- function(ptest){
  if(class(ptest) != 'pumping_test'){
    stop('A pumping_test object is required as input')
  }
  n <- floor(length(ptest$t)/3)
  ptest1 <- ptest
  ptest1$t <- ptest$t[1:n]
  ptest1$s <- ptest$s[1:n]
  initial <- fit(ptest1,'theis')
  a <- initial$parameters$a
  t0 <- initial$parameters$t0
  #print(t0)
  # n <- floor(length(ptest$t)/3)
  endp <- length(ptest$t)
  # ptest1 <- as.data.frame(cbind(ptest$t[n:endp],ptest$s[n:endp]))
  # names(ptest1) <- c('t', 's')
  # initial <- lm(s ~ log10(t), data = ptest1)
  # t1 <- 10^(-initial$coefficients[1]/initial$coefficients[2])
  sm <- ptest$s[endp]
  rt <- exp(-sm/a*2.3/2+0.1)
  if(rt > 1){
    rt <- -log(sm/a*2.3/2)
  }
  res <- list(a = a, t0 = t0,
              lambda = unname(1/rt))
  return(res)
}
#' @title
#' hantush_jacob_calculate_parameters
#' @description
#' Function to calculate the hydraulic parameters of the Hantush-Jacob model
#' @param ptest A pumping_test object
#' @param par A list with the values of a (slope) and t0 (Intercept) of the straight line fitted
#'  to the drawdown data using the Cooper-Jacob approach,  and lambda is the leakage factor.
#' @param hydraulic Logical flag to indicate if hydraulic parameters are calculated. If False,
#' the the statistcal parameter are calculated.
#' @return
#' A list with
#' \itemize{
#' \item Tr = Transmissivity
#' \item Ss = Storage coefficient
#' \item Ka = Aquitard Hydraulic Conductivity
#' \item radius_influence
#' }
#' @family hantush_jacob functions
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @export
#' @references
#' Hantush, M. S. Flow to wells in aquifers separated by a semipervious layer. Journal
#' of Geophysical Research, 1967, 72, 1709.
#' @examples
#' data(hantush_jacob)
#' ptest <- pumping_test('Well1', Q = 6.309e-3, r = 3.048, t=hantush_jacob$t, s = hantush_jacob$s)
#' ptest$additional_parameters$B <- 6.096
#' sol0 <- hantush_jacob_solution_initial(ptest)
#' hantush.par <- hantush_jacob_calculate_parameters(ptest,sol0)
#' print(hantush.par)
hantush_jacob_calculate_parameters <- function(ptest, par, hydraulic = TRUE){
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
  thickness_aquitard <- ptest$additional_parameters$B
  if(hydraulic){
    a <- par$a
    t0 <- par$t0
    pos_valid_a <- abs(a) < 1.0e-15
    pos_valid_t0 <- abs(t0) < 1.0e-15
    if( sum(pos_valid_a) > 0 | sum(pos_valid_t0) > 0){
      stop('a or t0 or both are close to numerical precision')
    }
    lambda <- par$lambda
    Tr <- 0.1832339*Q/a
    Ss <- 2.2458394*Tr*t0/r^2
    beta <-  r*lambda
    Ka <- Tr*thickness_aquitard/beta^2
    res <- list(Tr = Tr, Ss = Ss, Ka = Ka)
  }
  else {
    Tr <- par$Tr
    Ss <- par$Ss
    Ka <- par$Ka
    a <- 0.1832339*Q/Tr
    t0 <- (Ss*r^2)/(2.2458394*Tr)
    beta <- sqrt(Tr*thickness_aquitard/Ka)
    lambda <- beta/r
    res <- list(a = a, t0 = t0, lambda = lambda)
  }
  return(res)
}
#' @title
#' hantush_jacob_solution
#' @description
#' Function to calculate the drawdown of a Hantush-Jacob model
#' @param ptest A pumping_test object
#' @param a Slope of the straight line fitted to the drawdown data using the
#' Cooper-Jacob approach
#' @param t0 Intercept of the straight line fitted to the drawdown data using the
#' Cooper-Jacob approach
#' @param lambda Leakage factor defined as
#' \deqn{\lambda = \sqrt{\frac{b^{'}T}{k_{h}^{'}}}}
#' where \eqn{b^{'}} and \eqn{k_{h}^{'}} are the thickness and the vertical hydraulic
#' conductivity of the aquitard, and \eqn{T} is the transmissivity of the underlying
#' aquifer.
#' @param t Numeric vector with the time values
#' @return
#' A numeric vector with the calculated drawdown
#' @family hantush_jacob functions
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @export
#' @references
#' Hantush, M. S. Flow to wells in aquifers separated by a semipervious layer. Journal
#' of Geophysical Research, 1967, 72, 1709.
#' @examples
#' data(hantush_jacob)
#' ptest <- pumping_test('Well1', Q = 6.309e-3, r = 3.048, t=hantush_jacob$t, s = hantush_jacob$s)
#' ptest$additional_parameters$B <- 6.096
#' sol0 <- hantush_jacob_solution_initial(ptest)
#' sol <- hantush_jacob_solution(ptest,sol0$a, sol0$t0, sol0$lambda, ptest$t)
#' print(sol)
hantush_jacob_solution <- function(ptest, a, t0, lambda, t){
  if(class(ptest) != 'pumping_test'){
    stop('A pumping_test object is required as input')
  }
  Q <- ptest$Q
  r <- ptest$r
  par <- list(a = a, t0 = t0, lambda = lambda)
  hydr_par <- hantush_jacob_calculate_parameters(ptest, par)
  Tr <- hydr_par$Tr
  Ss <- hydr_par$Ss
  beta <- r/lambda
  #print(c(Tr,Ss,beta))
  u <- (Ss*r^2)/(4*Tr*t)
  td <- 1/u
  par <- list(beta = beta, coeffs = ptest$coeffs)
  #print(beta)
  W_hantush_jacob <- stehfest_inversion(td, ptest$coeffs, hantush_jacob_WF_LT,beta)
  #hantush_jacob_well_function(td, par)
  s <- (Q/(4*pi*Tr))*W_hantush_jacob
  return(s)
}
#' @title
#' hantush_jacob_solution_logdt
#' @description
#' Function to calculate the derivative of drawdown with respect to the log of time
#' @param ptest A pumping_test object
#' @param a Slope of the straight line fitted to the drawdown data using the
#' Cooper-Jacob approach
#' @param t0 Intercept of the straight line fitted to the drawdown data using the
#' Cooper-Jacob approach
#' @param lambda Leakage factor defined as
#' \deqn{\lambda = \sqrt{\frac{b^{'}T}{k_{h}^{'}}}}
#' where \eqn{b^{'}} and \eqn{k_{h}^{'}} are the thickness and the vertical hydraulic
#' conductivity of the aquitard, and \eqn{T} is the transmissivity of the underlying
#' aquifer.
#' @param t Numeric vector with the time values
#' @return
#' This function returns a list with the x and y values of the derivative of drawdown with
#' respect to log of time
#' @family hantush_jacob functions
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @export
#' @references
#' Hantush, M. S. Flow to wells in aquifers separated by a semipervious layer. Journal
#' of Geophysical Research, 1967, 72, 1709.
#' @examples
#' data(hantush_jacob)
#' ptest <- pumping_test('Well1', Q = 6.309e-3, r = 3.048, t=hantush_jacob$t, s = hantush_jacob$s)
#' ptest$additional_parameters$B <- 6.096
#' sol0 <- hantush_jacob_solution_initial(ptest)
#' sol_dlogt <- hantush_jacob_solution_dlogt(ptest, sol0$a, sol0$t0, sol0$lambda, ptest$t)
#' print(sol_dlogt)
hantush_jacob_solution_dlogt <- function(ptest, a, t0, lambda, t){
  sb <- hantush_jacob_solution(ptest, a, t0, lambda, t)
  dl_sb <- log_derivative_central(t, sb)
  res <- dl_sb$y
  return(res)
}
#' @title
#' hantush_jacob_WF_LT
#' @description
#' Calculates the Hantush-Jacob solution (drawdown) in the Laplace domain
#' @param p Laplace parameter
#' @param arg1 Beta parameter. The beta parameter is defined as
#' \deqn{\beta = \frac{r}{4} \sqrt{\frac{k_{h}^{'}S_{s}^{'}}{TS}}}
#' where \eqn{T} and \eqn{S} are the transmissivity and storage coefficient of the
#' aquifer, and \eqn{k_{h}^{'}} and \eqn{S_{s}^{'}} are the vertical hydraulic
#' conductivity and the Storage coefficient of the aquitard respectively.
#' @param ... Additional parameters
#' @return
#' A numeric vector with the calculated drawdown in the Laplace domain
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @family hantush_jacob functions
#' @export
#' @references
#' Hantush, M. S. Flow to wells in aquifers separated by a semipervious layer. Journal
#' of Geophysical Research, 1967, 72, 1709.
#' @examples
#' p <- 0.5
#' beta <- 0.5
#' res <- hantush_jacob_WF_LT(p, beta)
#' print(res)
hantush_jacob_WF_LT <- function(p, arg1, ...){
  beta <- arg1
  pm <- abs(4*p+beta^2)
  sp <- sqrt(pm)
  return( (2.0/p)*besselK(sp,0))
}
## -u-pdu/dp
#' @title
#' hantush_jacob_WF_LT_dlogt
#' @description
#' Calculates the values of the derivative of drawdown with respect to logarithm of
#' time for the Hantush-Jacob solution in the Laplace domain
#' @param p Laplace parameter
#' @param arg1 Value of beta
#' @param ... Additional parameters
#' @return
#' Numeric vector with the values of the derivative of drawdown with respect to logarithm of
#' time for the Hantush-Jacob solution in the Laplace domain
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @family hantush_jacob functions
#' @export
#' @references
#' Hantush, M. S. Flow to wells in aquifers separated by a semipervious layer. Journal
#' of Geophysical Research, 1967, 72, 1709.
#' @examples
#' p <- 0.5
#' beta <- 0.5
#' res <- hantush_jacob_WF_LT_dlogt(p, beta)
#' print(res)
hantush_jacob_WF_LT_dlogt <- function(p, arg1, ...){
  u <- hantush_jacob_WF_LT(p, arg1, ...)
  beta <- arg1
  k0 <- besselK(sqrt(4.0*p+beta^2), 0)
  k1 <- besselK(sqrt(4.0*p+beta^2), 1)
  term1 <- 2.0*k0/p^2
  term2 <- 4.0*k1/(p*sqrt(4.0*p+beta^2))
  dudp <- term1 + term2
  res <- -u + p*dudp
  return(res)
}
