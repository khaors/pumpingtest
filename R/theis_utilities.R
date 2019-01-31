#' @section Theis functions:
#'
#' The Theis solution describes the drawdown variations in a confined, homogeneous, isotropic,
#' aquifer of infinite horizontal extent due to the water extraction from a well that fully
#' penetrates the aquifer.
#'
#' The Functions included here are:
#'
#' theis_well_function, theis_solution_initial, theis_calculate_parameters, theis_solution, theis_solution_dlogt, rtheis_solution
#'
#' @docType package
#' @name pumpingtest
NULL
#' @title
#' theis_well_function_small
#' @description
#' Function to calculate the theis well function for values of \eqn{u < 1}
#' @param u Dimensionless \eqn{u} variable defined as:
#' \deqn{ u=\frac{r^{2}S}{4Tt} }
#' where \eqn{r} is the distance between the pumping and observation well in m, \eqn{S} is the
#' storage coefficient (dimensionless), \eqn{T} is the Transmisivity in m2/s and \eqn{t} is the time in s
#' @return
#' This function returns the value of the dimensionless drawdown
#' @family theis functions
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @export
#' @references
#' Theis, C. The relation between the lowering of the piezometric surface and the rate
#' and duration of discharge of a well using groundwater storage. Transactions of the
#' American Geophysical Union, 1935, 2, 519-524
#' @examples
#' u <- seq(0.01, 1, 0.01)
#' Ws <- theis_well_function_small(u)
#' plot(1/u, Ws, type = "l", log = "xy")
theis_well_function_small <- function(u){
  a <- c(-0.57721566, 0.99999193, -0.24991055, 0.05519968, -0.00976004,
         0.00107857)
  W<- -log(u) + a[1] + a[2]*u + a[3]*u**2 + a[4]*u**3 + a[5]*u**4 + a[6]*u**5
  return(W)
}
#' @title
#' theis_well_function_large
#' @description
#' Function to calculate the theis well function for values of \eqn{u > 1}
#' @param u Dimensionless \eqn{u} variable defined as:
#' \deqn{ u=\frac{r^{2}S}{4Tt} }
#' where \eqn{r} is the distance between the pumping and observation well in m, \eqn{S} is the
#' storage coefficient (dimensionless), \eqn{T} is the Transmisivity in m2/s and \eqn{t} is the time in s
#' @return
#' This function returns the value of the dimensionless drawdown
#' @family theis functions
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @export
#' @references
#' Theis, C. The relation between the lowering of the piezometric surface and the rate
#' and duration of discharge of a well using groundwater storage. Transactions of the
#' American Geophysical Union, 1935, 2, 519-524
#' @examples
#' u <- seq(1.0, 20.0, 0.5)
#' res <- theis_well_function_large(u)
theis_well_function_large <- function(u){
  a <- c(2.334733, 0.250621)
  b <- c(3.330657, 1.681534)
  num <- exp(-u)*(u^2+a[1]*u+a[2])
  den <- u*(u^2+b[1]*u+b[2])
  W <- num/den
  return(W)
}
#' @title
#' exponential_integral
#' @description
#' Function to calculate the exponential integral used in the Theis well function
#' @param x Numeric value
#' @return
#' This function returns the value of the exponential integral
#' @family theis functions
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @export
#' @examples
#' x <- seq(0.01,10,0.05)
#' Ei <- exponential_integral(x)
#' plot(x, Ei, type = "p", log = 'y')
exponential_integral <- function(x){
  # NIST Handbook Oliver pag167
  res <- log(x) + 0.57721566
  for(i in 1:10){
    num <- x^i
    den <- i*factorial(i)
    res <- res + num/den
  }
  return(res)
}
#' @title
#' theis_well_function
#' @description
#' Function to evaluate the Theis well function:
#'\deqn{ W(u) = \int_{u}^{\infty} \frac{\exp{-u}}{u} du}
#' @param td Dimensionless time \eqn{t_n} defined as
#'\deqn{ t_d=\frac{4Tt}{r^{2}S} }
#' where \eqn{r} is the distance between the pumping and observation well in m, \eqn{S} is the
#' storage coefficient (dimensionless), \eqn{T} is the Transmisivity in m2/s and \eqn{t} is the time in s
#' @param par additional parameters
#' @return
#' This function returns the value of the Theis well function
#' @family theis functions
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @export
#' @references
#' Theis, C. The relation between the lowering of the piezometric surface and the rate
#' and duration of discharge of a well using groundwater storage. Transactions of the
#' American Geophysical Union, 1935, 2, 519-524
#' @examples
#' u <- seq(0.05,30,0.05)
#' td <- 1.0/u
#' W <- theis_well_function(td, par)
#' plot(td, W, type = "l", log = "xy")
theis_well_function <- function(td, par = NULL){
  res <- theis_well_function_cpp(td, 0.0, 0.0, 0.0)
  return(res)
}
#' @title
#' theis_solution_initial
#' @description
#' Function to estimate the initial parameters of the Theis solution
#' @param ptest Object of a pumping_test class
#' @return
#' A list with the values of \eqn{a} and \eqn{t0}
#' @family theis functions
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @export
#' @references
#' Theis, C. The relation between the lowering of the piezometric surface and the rate
#' and duration of discharge of a well using groundwater storage. Transactions of the
#' American Geophysical Union, 1935, 2, 519-524
#' @examples
#' # Data from a confined aquifer
#' data(theis)
#' ptest <- pumping_test('Well1', Q = 1.3e-3, r = 200, t = theis$t, s = theis$s)
#' res <- theis_solution_initial(ptest)
#' print(res)
theis_solution_initial <- function(ptest){
  if(class(ptest) != 'pumping_test'){
    stop('A pumping_test object is required as input')
  }
  res <- cooper_jacob_solution_initial(ptest)
  return(res)
}
#' @title
#' theis_calculate_parameters
#' @description
#' Function to calculate the Transmisivity, Storage Coefficient and Radius of influence
#' for the Theis model
#' @param ptest An object of the pumping_test class
#' @param par List with the values of a and t0 obtained using the theis_solution_initial
#' @param hydraulic Logical flag to indicate if hydraulic parameters are calculated. If False, the the statistcal parameter are calculated.
#' function
#' @return
#' A list with
#' \itemize{
#'   \item Tr: Transmisivity
#'   \item Ss: Storage coefficient
#'   \item radius_influence
#' }
#' @family theis functions
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @export
#' @references
#' Theis, C. The relation between the lowering of the piezometric surface and the rate
#' and duration of discharge of a well using groundwater storage. Transactions of the
#' American Geophysical Union, 1935, 2, 519-524
#' @examples
#' # Data from a confined aquifer
#' data(theis)
#' ptest <- pumping_test('Well1', Q = 1.3e-3, r = 200, t = theis$t, s = theis$s)
#' res <- theis_solution_initial(ptest)
#' # Calculate hydraulic parameters from fit parameters
#' theis_hydrpar <- theis_calculate_parameters(ptest, res)
#' print(theis_hydrpar)
#' # Calculate fit parameters from hydraulic parameters
#' res <- list(Tr = 1.5e-3, Ss = 2e-5)
#' theis_fitpar <- theis_calculate_parameters(ptest, res, FALSE)
theis_calculate_parameters <- function(ptest, par, hydraulic = TRUE){
  if(class(ptest) != 'pumping_test') {
    stop('A pumping_test object is expected as input')
  }
  if(class(par) != 'list'){
    stop('A list is expected as input')
  }
  Q <- ptest$Q
  r <- ptest$r
  if(hydraulic){
    t <- ptest$t
    a <- par$a
    t0 <- par$t0
    pos_valid_a <- abs(a) < 1.0e-15
    pos_valid_t0 <- abs(t0) < 1.0e-15
    if( sum(pos_valid_a) > 0 | sum(pos_valid_t0) > 0){
      warning('a or t0 or both are close to numerical precision')
      #t1 <- 10^(-initial$coefficients[1]/initial$coefficients[2])
      #parameters <- list(a = unname(initial$coefficients[2]), t0 = unname(t1))
      # Work with the original slope and intercept of the linear regression line
      A <- -log10(t0)*a # Intercept
      C <- a  #Slope
      Tr <- 0.1832339*Q/C
      Dd <- 10^(A/C + log10((ptest$r^2)/2.25))
      Ss <- Tr/Dd
      endp <- length(t)
      Ri <- 2.0*sqrt(Tr*t[endp]/Ss)
    }
    else{
      Tr <- 0.1832339*Q/a
      Ss <- 2.2458394*Tr*t0/r^2
      endp <- length(t)
      Ri <- 2.0*sqrt(Tr*t[endp]/Ss)
    }
    res <- list(Tr = Tr, Ss = Ss, radius_influence = Ri)
  }
  else {
    Tr <-par$Tr
    Ss <- par$Ss
    a <- 0.1832339*Q/Tr
    t0 <- (Ss*r^2)/(2.2458394*Tr)
    res <- list(a = a ,t0 = t0)
  }
  return(res)
}
#' @title
#' theis_solution
#' @description
#' Function to calculate the drawdown using the Theis solution
#' @param ptest A pumpint test object
#' @param a slope of the straight line fitted to the drawdown data using the
#' Cooper-Jacob approach
#' @param t0 intercept of the straight line fitted to the drawdown data using the
#' Cooper-Jacob approach
#' @param t Numeric vector with the time values
#' @return
#' A numeric vector with the drawdown
#' @family theis functions
#' @export
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @references
#' Theis, C. The relation between the lowering of the piezometric surface and the rate
#' and duration of discharge of a well using groundwater storage. Transactions of the
#' American Geophysical Union, 1935, 2, 519-524
#' @examples
#' # Data from a confined aquifer
#' data(theis)
#' ptest <- pumping_test('Well1', Q = 1.3e-3, r = 200, t = theis$t, s = theis$s)
#' res <- theis_solution_initial(ptest)
#' s_model <- theis_solution(ptest,res$a, res$t0, theis$t)
#' #Compare drawdown from model with measurements
#' mse <- sqrt(mean((s_model-theis$s)^2))
#' print(mse)
theis_solution <- function(ptest, a, t0, t){
  if(ptest$Q < 1.0e-10){
    stop("The discharge rate is not defined")
  }
  if(ptest$r < 1.0e-10){
    stop("The distance to the observation well is not defined")
  }
  Q <- ptest$Q
  r <- ptest$r
  invalid_a <- abs(a) < 1.0e-15
  invalid_t0 <- abs(t0) < 1.0e-15
  if(sum(invalid_a) > 0 | sum(invalid_t0) > 0){
    warning('a or t0 or both are close to numerical precision')
  }
  input_par <- list(a = a, t0 = t0)
  par <- theis_calculate_parameters(ptest, input_par)
  par$coeffs <- ptest$coeffs
  Tr <- par$Tr
  Ss <- par$Ss
  u <- (Ss*r**2)/(4*Tr*t)
  td <- 1.0/u
  W <- theis_well_function(td, par)
  s <- (Q/(4*pi*Tr))*W
  return(s)
}
#' @title
#' theis_WF_LT
#' @description
#' Function to calculate the drawdown in the Laplace domain
#' @param p Laplace parameter
#' @param ... Additional parameters
#' @return
#' A numeric vector with the calculated drawdown in Laplace domain
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail}
#' @export
#' @family theis functions
#' @references
#' Theis, C. The relation between the lowering of the piezometric surface and the rate
#' and duration of discharge of a well using groundwater storage. Transactions of the
#' American Geophysical Union, 1935, 2, 519-524
#' @examples
#' coeffs <- stehfest_coeffs(8)
#' td <- logseq(-2,4,50)
#' W <- stehfest_inversion(td, coeffs, theis_WF_LT)
#' plot(td, W, log = "xy")
theis_WF_LT <- function(p, ...){
  sp <- sqrt(p)
  sd <- (2.0/p)*besselK(2.0*sp, 0)
  return(sd)
}
#' @title
#' theis_WF_LT_dlogt
#' @description
#' Function to calculate the derivative of drawdown wrt log of time in the Laplace domain
#' @param p Laplace parameter
#' @param ... Additional parameters
#' @return
#' A numeric vector with the derivative of drawdown wrt log of time
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @family theis functions
#' @export
#' @references
#' Theis, C. The relation between the lowering of the piezometric surface and the rate
#' and duration of discharge of a well using groundwater storage. Transactions of the
#' American Geophysical Union, 1935, 2, 519-524
#' @examples
#' coeffs <- stehfest_coeffs(8)
#' td <- logseq(-0.5,4,50)
#' Wdlogt <- stehfest_inversion(td, coeffs, theis_WF_LT_dlogt)
#' plot(td, Wdlogt, log = "xy")
theis_WF_LT_dlogt <- function(p, ...){
  u <- theis_WF_LT(p, ...)
  sp <- sqrt(p)
  dudp <- 2.*besselK(2.0*sp, 0)/(p^2) + 2.0*besselK(2.0*sp, 1)/(p^(3/2))
  res <- -u+p*dudp
  return(res)
}
#' @title
#' theis_solution_dlogt
#' @description
#' Function to calculate the derivative of the drawdown with respect to the log of time
#' @param ptest A pumping_test object
#' @param a Slope of the straight line fitted to the drawdown data using the
#' Cooper-Jacob approach
#' @param t0 Intercept of the straight line fitted to the drawdown data using the
#' Cooper-Jacob approach
#' @param t value
#' @return
#' This function returns the derivative of the drawdown with respect to the
#' derivative of log time
#' @family theis functions
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @export
#' @references
#' Theis, C. The relation between the lowering of the piezometric surface and the rate
#' and duration of discharge of a well using groundwater storage. Transactions of the
#' American Geophysical Union, 1935, 2, 519-524
#' @examples
#' data(theis)
#' ptest <- pumping_test("Well1", Q = 1.3e-3, r = 200, t = theis$t, s = theis$s)
#' sol0 <- theis_solution_initial(ptest)
#' resd <- theis_solution_dlogt(ptest, sol0$a, sol0$t0, theis$t)
#' print(resd)
theis_solution_dlogt <- function(ptest, a, t0, t){
  Q <- ptest$Q
  r <- ptest$r
  par <- list(a = a, t0 = t0)
  res <- theis_calculate_parameters(ptest, par)
  Tr <- res$Tr
  Ss <- res$Ss
  u <- (Ss*r**2)/(4.0*Tr*t)
  #res <- (Q/(4*pi*Tr))*((1-u)/t)
  #res <- (Q/(4.0*pi*Tr))*exp(-u)
  td <- 1/u
  Wdlogt <- stehfest_inversion(td, ptest$coeffs, theis_WF_LT_dlogt)
  res <- (Q/(4.0*pi*Tr))*Wdlogt
  return(res)
}
