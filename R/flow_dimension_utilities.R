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
#' This function returns a list with components:
#' \itemize{
#' \item x: Numeric vector with times at which the second derivative is evaluated
#' \item y: Numeric vector with the values of the second derivative of drawdown
#' \item n: Numeric vector with the flow dimension values evaluated at each time
#' }
#' @author
#' Oscar Garcia-Cabrejo, \email{khaors@gmail.com}
#' @family flow_dimension functions
#' @export
#' @examples
#' data(boulton)
#' t <- boulton$t
#' s <- boulton$s
#' boulton.fdim <- flow_dimension(t,s, method = "smoothspline")
#' plot(boulton.fdim$x, boulton.fdim$n, type = "p", log = "x", 
#' ylim = c(0, 10))
flow_dimension <- function(t, s, d = 2, method = "central"){
  flow_dim <- NULL
  if(method == 'central'){
    res1a <- log_derivative_central(t, s, return.pos = F, log = F)
    pos <- res1a$y > 0
    res1b <- log_derivative_central(res1a$x[pos],
                                    log(res1a$y[pos]*log(10)),
                                     return.pos = F)
    flow_dim$x <- res1b$x
    flow_dim$y <- res1b$y
    flow_dim$n <- 2-2*res1b$y
  }
  else if(method == 'horner'){
    res1a <- log_derivative_horner(t, s, return.pos = F, log = F)
    pos <- res1a$y > 0
    res1b <- log_derivative_horner(res1a$x[pos],
                                   log(res1a$y[pos]*log(10)),
                                    return.pos = F)
    flow_dim$x <- res1b$x
    flow_dim$y <- res1b$y
    flow_dim$n <- 2-2*res1b$y
  }
  else if(method == 'bourdet'){
    res1a <- log_derivative_bourdet(t, s, return.pos = F, 
                                    log = F, d = d)
    pos <- res1a$y > 0.0
    res1b <- log_derivative_bourdet(res1a$x[pos], 
                                    log(res1a$y[pos]*log(10)),
                                     return.pos = F, d = d)
    flow_dim$x <- res1b$x
    flow_dim$y <- res1b$y
    flow_dim$n <- 2-2*res1b$y
  }
  # else if(method =='spline'){
  #   log_d <- log_derivative_spline(t, s, n = d)
  #   y <- log(log_d$y)
  #   flow_dim <- log_derivative_spline(t, y, n = d)
  #   flow_dim$n <- 2-2*flow_dim$y
  # }
  else if(method == 'spane'){
    log_d <- log_derivative_spane(t, s, n = d, return.pos = F)
    y <- log(log_d$y)
    flow_dim <- log_derivative_horner(t, y, return.pos = F)
    flow_dim$n <- 2-2*flow_dim$y
  }
  else if(method == 'smoothspline'){
    res1a <- log_derivative_smoothspline(t, s, return.pos = F, 
                                         log = F)
    pos <- res1a$y > 0 
    res1b <- log_derivative_smoothspline(res1a$x[pos],
                                         log(res1a$y[pos]*log(10)),
                                         return.pos = F)
    flow_dim$x <- res1b$x
    flow_dim$y <- res1b$y
    flow_dim$n <- 2-2*res1b$y
  }
  else if(method == 'kernelreg'){
    log_d <- log_derivative_kernelreg(t, s)
    y <- log(log_d$y)
    flow_dim <- log_derivative_kernelreg(t, y)
    flow_dim$n <- 2-2*flow_dim$y
  }
  else if(method == 'lokern'){
    res1a <- log_derivative_lokern(t, s, return.pos = F, log = F)
    pos <- res1a$y > 0
    res1b <- log_derivative_lokern(res1a$x[pos], 
                                   log(res1a$y[pos]*log(10)),
                                    return.pos = F)
    flow_dim$x <- res1b$x
    flow_dim$y <- res1b$y
    flow_dim$n <- 2-2*res1b$y
  }
  else if(method == 'locpol'){
    res1a <- log_derivative_locpol(t, s, return.pos = F, log = F)
    pos <- res1a$y > 0
    res1b <- log_derivative_locpol(res1a$x[pos], 
                                   log(res1a$y[pos]*log(10)),
                                    return.pos = F)
    flow_dim$x <- res1b$x
    flow_dim$y <- res1b$y
    flow_dim$n <- 2-2*res1b$y
  }
  else if(method == "lpridge"){
    res1a <- log_derivative_lpridge(t, s, return.pos = F, log = F)
    pos <- res1a$y > 0
    res1b <- log_derivative_lpridge(res1a$x[pos], 
                                    log(res1a$y[pos]*log(10)),
                                   return.pos = F)
    flow_dim$x <- res1b$x
    flow_dim$y <- res1b$y
    flow_dim$n <- 2-2*res1b$y   
  }
  else {
    stop("ERROR: Unknown derivative type. Please check and try again")
  }
  return(flow_dim)
  
  
}