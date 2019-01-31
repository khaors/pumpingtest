#' @section Neuzil functions:
#'
#' These functions are used in the estimation of pulse tests.
#'
#' The functions included in this section are:
#' neuzil_well_function, neuzil_solution_initial, neuzil_calculate_parameters, neuzil_solution, neuzil_WF_LT, neuzil_solution_dlogt
#'
#' @docType package
#' @name pumpingtest
NULL
#' @title
#' neuzil_well_function
#' @description
#' Calculates the well function of the Neuzil solution
#' @param td Numeric vector with the dimensionless time
#' @param par A list with the wellbore storage coefficient (cd) and stehfest coeffs (coeffs)
#' @return
#' A numeric vector with the dimensionless drawdown calculated using the Neuzil solution
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @export
#' @family neuzil functions
#' @references
#' Neuzil, C. E. On conducting the modified 'Slug' test in tight formations. Water
#' Resources Research, 1982, 18, 2, 439-441
#' @examples
#' td <- logseq(-1, 3, 50)
#' par <- list(cd = 0.01, coeffs = stehfest_coeffs(8) )
#' W <- neuzil_well_function(td, par)
#' plot(td, W, type = "l", log = "xy", main = "NEUZIL SOLUTION: Well function")
neuzil_well_function <- function(td, par){
  arg1 <- par$cd
  td <- td*arg1
  coeffs <- par$coeffs
  sd <- stehfest_inversion(td, coeffs, neuzil_WF_LT, arg1 = arg1)
  return(sd)
}
#' @title
#' neuzil_solution_initial
#' @description
#' Calculates the initial values of parameters of Neuzil solution
#' @param ptest A pumping_test object
#' @return
#' A list with
#' \itemize{
#' \item cd: Wellbore storage coefficient
#' \item t0: Intercept
#' }
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @export
#' @family neuzil functions
#' @references
#' Neuzil, C. E. On conducting the modified 'Slug' test in tight formations. Water
#' Resources Research, 1982, 18, 2, 439-441
#' @examples
#' data(neuzil)
#' ptest.neuzil <- pumping_test("Well1", Q = 0, r = 0.0, t = neuzil$t, s = neuzil$s)
#' par <- list(rw = 0.067, ceff = 2.723e-9, vs = 1.59)
#' ptest.neuzil$additional_parameters <- par
#' sol0 <- neuzil_solution_initial(ptest.neuzil)
#' print(sol0)
neuzil_solution_initial <- function(ptest){
  if(class(ptest) != 'pumping_test'){
    stop('A pumping_test object is required as input')
  }
  res <- cooper_solution_initial(ptest)
  if(res$cd > 1e6){
    res$cd <- 50
  }
  return(res)
}
#' @title
#' neuzil_calculate_parameters
#' @description
#' Function to calculate the hydraulic parameters for the Neuzil solution
#' @param ptest A pumping_test object
#' @param par A list with the wellbore storage coefficient (cd) and t0 parameters
#' @param hydraulic Logical flag to indicate if hydraulic parameters are calculated.
#' If False, the the statistcal parameter (a and t0) are calculated.
#' @return
#' A list with
#' \itemize{
#' \item Tr: Transmissivity
#' \item Ss: Storage coefficient
#' \item radius_influence
#' }
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @export
#' @family neuzil functions
#' @references
#' Neuzil, C. E. On conducting the modified 'Slug' test in tight formations. Water
#' Resources Research, 1982, 18, 2, 439-441
#' @examples
#' data(neuzil)
#' ptest.neuzil <- pumping_test("Well1", Q = 0, r = 0.0, t = neuzil$t, s = neuzil$s)
#' par <- list(rw = 0.067, ceff = 2.723e-9, vs = 1.59)
#' ptest.neuzil$additional_parameters <- par
#' sol0 <- neuzil_solution_initial(ptest.neuzil)
#' neuzil.par <- neuzil_calculate_parameters(ptest.neuzil, sol0)
#' print(neuzil.par)
neuzil_calculate_parameters <- function(ptest, par, hydraulic = TRUE){
  if(class(ptest) != 'pumping_test'){
    stop('A pumping_test object is required as input')
  }
  if(class(par) != 'list'){
    stop('a list is required as input')
  }
  if(is.null(ptest$additional_parameters)){
    stop('The additional_parameters list is empty')
  }
  ceff <- ptest$additional_parameters$ceff
  vs <- ptest$additional_parameters$vs
  rw <- ptest$additional_parameters$rw
  rho <- 1000
  g <- 9.80
  sw <- vs*ceff*rho*g
  if(hydraulic){
    cd <- par$cd
    t0 <- par$t0
    pos_valid_cd <- abs(cd) < 1.0e-15
    pos_valid_t0 <- abs(t0) < 1.0e-15
    if( sum(pos_valid_cd) > 0 | sum(pos_valid_t0) > 0){
      stop('cd or t0 or both are close to numerical precision')
    }
    #Compute the transmissivity and storativity
    Tr <- 0.5*sw/(t0*pi)
    Ss <- 0.5*sw/(cd*pi*rw^2)
    # Compute the radius of investigation according to Guyonnet et al. (1993)
    n <- 0.462
    m <- 0.495
    xlim <- 2.5
    endp <- length(ptest$t)
    tdl <- cd*t(endp)/t0
    x <- (tdl^n)/(cd^m);
    if(x < xlim){
      radius_influence <- rw*3.54*(tdl^n)
    }
    else{
      radius_influence <- rw*8.37*(cd^m)
    }
    res <- list(Tr = Tr, Ss = Ss, radius_influence = radius_influence)
  }
  else {
    Tr <- par$Tr
    Ss <- par$Ss
    t0 <- 0.5*sw/(Tr*pi)
    cd <- 0.5*sw/(Ss*pi*rw^2)
    res <- list(t0 = t0, cd = cd)
  }
  return(res)
}
#' @title
#' neuzil_solution
#' @description
#' Calculates the drawdown using Neuzil solution
#' @param ptest A pumping_test object
#' @param cd Value of the wellbore storage coefficient
#' @param t0 Intercept
#' @param t Numeric vector with time
#' @return
#' A numeric vector with the calculated drawdown
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @export
#' @family neuzil functions
#' @references
#' Neuzil, C. E. On conducting the modified 'Slug' test in tight formations. Water
#' Resources Research, 1982, 18, 2, 439-441
#' @examples
#' data(neuzil)
#' ptest.neuzil <- pumping_test("Well1", Q = 0, r = 0.0, t = neuzil$t, s = neuzil$s)
#' par <- list(rw = 0.067, ceff = 2.723e-9, vs = 1.59)
#' ptest.neuzil$additional_parameters <- par
#' sol0 <- neuzil_solution_initial(ptest.neuzil)
#' sol <- neuzil_solution(ptest.neuzil, sol0$cd, sol0$t0, ptest.neuzil$t)
#' print(sol)
neuzil_solution <- function(ptest, cd, t0, t){
  if(class(ptest) != 'pumping_test'){
    stop('A pumping_test object is required as input')
  }
  td <- t/t0
  par <- list(cd = cd, coeffs = ptest$coeffs)
  s <- neuzil_well_function(td, par)
  return(s)
}
#' @title
#' neuzil_WF_LT
#' @description
#' Calculates the drawdown using Neuzil function in the Laplace domain
#' @param p Laplace parameter
#' @param arg1 Value of the wellbore storage coefficient (cd)
#' @param ... Additional parameters
#' @return
#' A numeric vector with the calculated drawdown in the Laplace domain
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @export
#' @family neuzil functions
#' @references
#' Neuzil, C. E. On conducting the modified 'Slug' test in tight formations. Water
#' Resources Research, 1982, 18, 2, 439-441
#' @examples
#' p <- 0.5
#' cd <- 0.05
#' s <- neuzil_WF_LT(p, cd)
#' print(s)
neuzil_WF_LT <- function(p, arg1, ...){
  cd <- arg1
  sp <- sqrt(abs(p))
  k0 <- besselK(sp, 0)
  k1 <- besselK(sp, 1)
  s <- (cd*k0)/(cd*p*k0+sp*k1)
  return(s)
}
#' @title
#' neuzil_solution_dlogt
#' @description
#' Function to calculates the derivative of neuzil solution
#' @param ptest A pumping_test object
#' @param cd Value of the wellbore storage coefficient (cd)
#' @param t0 Intercept
#' @param t A Numeric vector with time
#' @return
#' A numeric vector with the derivative of drawdown
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @export
#' @family neuzil functions
#' @references
#' Neuzil, C. E. On conducting the modified 'Slug' test in tight formations. Water
#' Resources Research, 1982, 18, 2, 439-441
#' @examples
#' data(neuzil)
#' ptest.neuzil <- pumping_test("Well1", Q = 0, r = 0.0, t = neuzil$t, s = neuzil$s)
#' par <- list(rw = 0.067, ceff = 2.723e-9, vs = 1.59)
#' ptest.neuzil$additional_parameters <- par
#' sol0 <- neuzil_solution_initial(ptest.neuzil)
#' sol_dlogt <- neuzil_solution_dlogt(ptest.neuzil, sol0$cd, sol0$t0, ptest.neuzil$t)
#' print(sol_dlogt)
neuzil_solution_dlogt <- function(ptest, cd, t0, t){
  sb <- neuzil_solution(ptest, cd, t0, t)
  dl_sb <- log_derivative_central(t, -sb)
  res <- dl_sb$y
  return(res)
}