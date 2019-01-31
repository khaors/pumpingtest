#' @section Flow Dimension functions:
#' flow_dimension
#' @docType package
#' @name pumpingtest
NULL
#' @title
#' flow_dimension
#' @description
#' Function to calculate the flow dimension using different approaches
#' @param t Numeric vector with time values
#' @param s Numeric vector with drawdown values
#' @param d Derivative parameter. if method equals to bourdet then d is equal to
#' the number of adjacent values used in the derivative calculation. If method is
#' equal to spline then d is equal to the number of knots used in the interpolation
#' of the drawdown data. In this case a value of d=20 to d=30 is recommended. If
#' method is equal to spane then d is equal to the number of points used in the linear
#' regression approach.
#' @param method Method to calculate the derivative. See log_derivative
#' @return
#' This function returns a list with components named as x and y that contains the
#' flow dimension y evaluated at specific points x.
#' @author
#' Oscar Garcia-Cabrejo, \email{khaors@gmail.com}
#' @family flow_dimension functions
#' @export
flow_dimension <- function(t, s, d = 2, method = "central"){
  flow_dim <- NULL
  if(method == 'central'){
    log_d <- log_derivative_central(t, s)
    y <- log(log_d$y)
    flow_dim <- log_derivative_central(t, y)
    flow_dim$n <- 2-2*flow_dim$y
  }
  else if(method == 'horner'){
    log_d <- log_derivative_horner(t, s)
    y <- log(log_d$y)
    flow_dim <- log_derivative_horner(t, y)
    flow_dim$n <- 2-2*flow_dim$y
  }
  else if(method == 'bourdet'){
    log_d <- log_derivative_bourdet(t, s, d)
    y <- log(log_d$y)
    flow_dim <- log_derivative_bourdet(t, y, d)
    flow_dim$n <- 2-2*flow_dim$y
  }
  else if(method =='spline'){
    log_d <- log_derivative_spline(t, s, n = d)
    y <- log(log_d$y)
    flow_dim <- log_derivative_spline(t, y, n = d)
    flow_dim$n <- 2-2*flow_dim$y
  }
  else if(method == 'spane'){
    log_d <- log_derivative_spane(t, s, n = d)
    y <- log(log_d$y)
    flow_dim <- log_derivative_spane(t, y, n = d)
    flow_dim$n <- 2-2*flow_dim$y
  }
  else if(method == 'smoothspline'){
    log_d <- log_derivative_smoothspline(t, s)
    y <- log(log_d$y)
    flow_dim <- log_derivative_smoothspline(t, y)
    flow_dim$n <- 2-2*flow_dim$y
  }
  else if(method == 'kernelreg'){
    log_d <- log_derivative_kernelreg(t, s)
    y <- log(log_d$y)
    flow_dim <- log_derivative_kernelreg(t, y)
    flow_dim$n <- 2-2*flow_dim$y
  }
  else if(method == 'lokern'){
    log_d <- log_derivative_lokern(t, s)
    y <- log(log_d$y)
    flow_dim <- log_derivative_lokern(t, y)
    flow_dim$n <- 2-2*flow_dim$y
  }
  else if(method == 'locpol'){
    log_d <- log_derivative_lokern(t, s)
    y <- log(log_d$y)
    flow_dim <- log_derivative_locpol(t, y)
    flow_dim$n <- 2-2*flow_dim$y
  }
  else if(method == "lpridge"){
    log_d <- log_derivative_lpridge(t, s)
    y <- log(log_d$y)
    flow_dim <- log_derivative_lpridge(t, y)
    flow_dim$n <- 2-2*flow_dim$y
  }
  else {
    stop("ERROR: Unknown derivative type. Please check and try again")
  }
  return(flow_dim)
  
  
}