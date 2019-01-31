#' @section Cooper-Jacob functions:
#' cooper_jacob_well_function, cooper_jacob_solution_initial, cooper_jacob_calculate_parameters, cooper_jacob_solution
#' @docType package
#' @name pumpingtest
NULL
#' @title
#' cooper_jacob_well_function
#' @description
#' Function to calculate the well function of the Cooper-Jacob solution
#' @param u dimensionless variable u
#' @return
#' This function returns the dimensionless drawdown using the Cooper-Jacob solution
#' @family cooper_jacob functions
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @export
#' @references
#' Cooper, H. & Jacob, C. A generalized graphical method for evaluating formation
#' constants and summarizing well field history. Transactions of the American
#' Geophysical Union, 1946, 27, 526-534
#' @examples
#' td <- logseq(-2, 4, 50)
#' u <- 1.0/td
#' res <- cooper_jacob_well_function(u)
#' plot(u,res, type="p", log="xy")
cooper_jacob_well_function <- function(u){
  pos_valid <- u > 0.3
  #if(sum(pos_valid) > 0){
  #  warning("The value of u > 0.3. The Cooper-Jacob solution is not accurate")
  #}
  return(-0.5772157-log(u))
}
#' @title
#' cooper_jacob_solution_initial
#' @description
#' Function to calculate the initial values of the Cooper-Jacob Solution
#' @param ptest A pumping_test object
#' @return
#' A list with the values of a (slope) and t0 (intercept) of the straight line fitted to the
#' drawdown data
#' @family cooper_jacob functions
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @export
#' @references
#' Cooper, H. & Jacob, C. A generalized graphical method for evaluating formation
#' constants and summarizing well field history. Transactions of the American
#' Geophysical Union, 1946, 27, 526-534
#' @importFrom stats lm
#' @examples
#' data(theis)
#' ptest <- pumping_test("Well1", Q= 1.3e-3, r = 200, t = theis$t, s = theis$s)
#' cooper_jacob_sol0 <- cooper_jacob_solution_initial(ptest)
#' print(cooper_jacob_sol0)
cooper_jacob_solution_initial <- function(ptest){
  n <- 1 #floor(length(ptest$t)/3)
  endp <- length(ptest$t)
  ptest1 <- as.data.frame(cbind(ptest$t[n:endp],ptest$s[n:endp]))
  names(ptest1) <- c('t', 's')
  initial <- lm(s ~ log10(t), data = ptest1)
  t1 <- 10^(-initial$coefficients[1]/initial$coefficients[2])
  parameters <- list(a = unname(initial$coefficients[2]), t0 = unname(t1))
  res <- parameters
  return(res)
}
#' @title
#' cooper_jacob_calculate_parameters
#' @description
#' Function to calculate the transmissivity, storage coefficient and radius of influence
#' using the Cooper-Jacob solution
#' @param ptest A pumping_test object
#' @param par list with the value of a (slope) and t0 (Intercept) of the straight line
#' fitted to the drawdown data
#' @param hydraulic Logical flag to indicate if hydraulic parameters are calculated. If
#' False, the the statistcal parameter are calculated.
#' @return
#' A list with the transmissivity, storage coefficient and radious of influence
#' @family cooper_jacob functions
#' @export
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @references
#' Cooper, H. & Jacob, C. A generalized graphical method for evaluating formation
#' constants and summarizing well field history. Transactions of the American
#' Geophysical Union, 1946, 27, 526-534
#' @examples
#' data(theis)
#' ptest <- pumping_test("Well1", Q= 1.3e-3, r = 200, t = theis$t, s = theis$s)
#' cooper_jacob_sol0 <- cooper_jacob_solution_initial(ptest)
#' par <- cooper_jacob_calculate_parameters(ptest, cooper_jacob_sol0)
#' print(par)
cooper_jacob_calculate_parameters <- function(ptest, par, hydraulic = TRUE){
  if(class(ptest) != 'pumping_test') {
    stop('A pumping_test object is expected as input')
  }
  if(class(par) != 'list'){
    stop('A list is expected as input')
  }
  Q <- ptest$Q
  r <- ptest$r
  t <- ptest$t
  if(hydraulic){
    a <- par$a
    t0 <- par$t0
    invalid_a <- abs(a) < 1.0e-15
    invalid_t0 <- abs(t0) < 1.0e-15
    if(sum(invalid_a) > 0 | sum(invalid_t0) > 0){
      stop('a or t0 or both are close to numerical precision')
    }
    endp <- length(t)
    Tr <- 0.1832339*Q/a
    Ss <- 2.2458394*Tr*t0/r^2
    #print(t[endp])
    Ri <- 2*sqrt(Tr*t[endp]/Ss)
    res <- list(Tr = Tr, Ss = Ss, radius_influence = Ri)
  }
  else {
    Tr <- par$Tr
    Ss <- par$Ss
    a <- 0.1832339*Q/Tr
    t0 <- (Ss*r^2)/(2.2458394*Tr)
    res <- list(a = a, t0 = t0)
  }
  return(res)
}
#' @title
#' cooper_jacob_solution
#' @description
#' Function to calculate the drawdown using the Cooper-Jacob solution for a confined aquifer
#' @param ptest A pumping_test object
#' @param a Slope of the straight line fitted to the drawdown data using the
#' Cooper-Jacob approach
#' @param t0 Intercept of the straight line fitted to the drawdown data using the
#' Cooper-Jacob approach
#' @param t Numeric vector with the values of time
#' @return
#' This function returns a numeric vector with the calculated drawdown
#' @family cooper_jacob functions
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @export
#' @references
#' Cooper, H. & Jacob, C. A generalized graphical method for evaluating formation
#' constants and summarizing well field history. Transactions of the American
#' Geophysical Union, 1946, 27, 526-534
#' @examples
#' data(theis)
#' ptest <- pumping_test("Well1", Q= 1.3e-3, r = 200, t = theis$t, s = theis$s)
#' sol0 <- cooper_jacob_solution_initial(ptest)
#' sol <- cooper_jacob_solution(ptest, sol0$a, sol0$t0, ptest$t)
#' print(sol)
cooper_jacob_solution <- function(ptest, a, t0, t){
  Q <- ptest$Q
  r <- ptest$r
  if(abs(a) < 1.0e-15 | abs(t0) < 1.0e-15){
    stop('a or t0 or both are close to numerical precision')
  }
  par <- list(a = a, t0 = t0)
  parameters <- cooper_jacob_calculate_parameters(ptest, par)
  Tr <- parameters$Tr
  Ss <- parameters$Ss
  u <- (Ss*r**2)/(4*Tr*t)
  pos_warning <- u > 0.3
  #if(sum(pos_warning) > 0){
  #  warning("The value of u > 0.3. The Cooper-Jacob solution is not accurate")
  #}
  W <- cooper_jacob_well_function(u)
  s <- (Q/(4.0*pi*Tr))*W;
  return(s)
}
#' @title
#' cooper_jacob_solution_dlogt
#' @description
#' Function to calculate the derivative of the drawdown with respect to log of time of the
#' Cooper-Jacob solution
#' @param ptest A pumping_test object
#' @param a Slope of the straight line fitted to the drawdown data using the
#' Cooper-Jacob approach
#' @param t0 Intercept of the straight line fitted to the drawdown data using the
#' Cooper-Jacob approach
#' @param t Numeric vector with the time at which measurements were taken
#' @return
#' A numeric vector with the values of the derivative of drawdown with respect to log of time
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @family cooper_jacob functions
#' @export
#' @references
#' Cooper, H. & Jacob, C. A generalized graphical method for evaluating formation
#' constants and summarizing well field history. Transactions of the American
#' Geophysical Union, 1946, 27, 526-534
#' @examples
#' data(theis)
#' ptest <- pumping_test("Well1", Q= 1.3e-3, r = 200, t = theis$t, s = theis$s)
#' sol0 <- cooper_jacob_solution_initial(ptest)
#' sol_dlogt <- cooper_jacob_solution_dlogt(ptest, sol0$a, sol0$t0, ptest$t)
#' print(sol_dlogt)
cooper_jacob_solution_dlogt <- function(ptest, a, t0, t){
  if(abs(a) < 1.0e-15 | abs(t0) < 1.0e-15){
    stop('a or t0 or both are close to numerical precision')
  }
  n <- length(t)
  res <- vector("numeric", length = n)
  hydr_par <- cooper_jacob_calculate_parameters(ptest,list(a= a, t0 = t0))
  Tr <- hydr_par$Tr
  res[1:n] <- 1.0*(ptest$Q/(4*pi*Tr))
  return(res)
}