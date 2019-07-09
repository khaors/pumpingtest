#' @section Utility functions:
#' logseq, log_derivative, log_derivative_central, log_derivative_bourdet, log_derivative_spline, log_derivative_spane,
#' log_derivative_smoothspline, log_derivative_kernel_reg, log_derivative_locpol, log_derivative_lokern,
#' log_derivative_lpridge
#' @docType package
#' @name pumpingtest
NULL
#' @title
#' logseq
#' @description
#' Function to generate a sequence with log increments
#' @param from,to log10 of initial and final points in the sequence
#' @param n step
#' @return
#' Sequence with log increments
#' @author 
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @family utility functions
#' @export
#' @examples
#' s <- logseq(0, 2, 10)
logseq <- function( from = 1, to = 1, n) {
  #lfrom <- log10(from)
  #lto <- log10(to)
  exp(log(10)*seq(from, to, length.out = n))
}
#' @title
#' log_derivative
#' @description
#' Function to calculate the derivative of drawdown wrt the log of time
#' @param t vector with time values
#' @param s vector with drawdown values
#' @param d Derivative parameter. if method equals to bourdet then d is equal to the number of
#' adjacent values used in the derivative calculation. If method is equal to spline then d is equal
#' to the number of knots used in the interpolation of the drawdown data. In this case a value of d=20 to
#' d=30 is recommended. If method is equal to spane then d is equal to the number of points used in the
#' linear regression approach.
#' @param method Method to calculate the derivative (central, horner, bourdet and
#' spline)
#' @param return.pos Logical flag to return only the positive values of the log-derivative (default = TRUE)
#' @param log Logical flag to indicate that natural logarithm (a log to the base e) is used in the derivative 
#' calculation (default = TRUE). Logarithm to the base 10 is used if FALSE.
#' @return
#' This function returns a list with components named as x and y that contains the
#' log_derivative y evaluated at specific points x.
#' @family log_derivative functions
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @export
#' @examples
#' # Load test data
#' data(theis)
#' # calculate derrivative of s wrt to logt using central differences
#' logd.central <- log_derivative(theis$t, theis$s)
#' # calculate derrivative of s wrt to logt using Horner method
#' logd.horner <- log_derivative(theis$t, theis$s, method = "horner")
#' # calculate derrivative of s wrt to logt using Bourdet method
#' logd.bourdet <- log_derivative(theis$t, theis$s, method = "bourdet")
#' # plot results
#' plot(logd.central$x, logd.central$y, type = "p")
#' points(logd.horner$x, logd.horner$y, col = "red")
#' points(logd.bourdet$x, logd.bourdet$y, col = "blue")
log_derivative <- function(t, s, d = 2, method = 'central', return.pos = T, 
                           log = T){
  if(method == 'central'){
    log_d <- log_derivative_central(t, s, return.pos, log)
  }
  else if(method == 'horner'){
    log_d <- log_derivative_horner(t, s, return.pos, log)
  }
  else if(method == 'bourdet'){
    log_d <- log_derivative_bourdet(t, s, d, return.pos, log)
  }
  else if(method =='spline'){
    log_d <- log_derivative_spline(t, s, n = d, log)
  }
  else if(method == 'spane'){
    log_d <- log_derivative_spane(t, s, n = d, return.pos, log)
  }
  else if(method == 'smoothspline'){
    log_d <- log_derivative_smoothspline(t, s, return.pos, log)
  }
  else if(method == 'kernelreg'){
    log_d <- log_derivative_kernelreg(t, s, bw = d, return.pos, log)
  }
  else if(method == 'lokern'){
    log_d <- log_derivative_lokern(t, s, return.pos, log)
  }
  else if(method == 'locpol'){
    log_d <- log_derivative_locpol(t, s, return.pos, log)
  }
  else if(method == "lpridge"){
    log_d <- log_derivative_lpridge(t, s, return.pos, log)
  }
  #  else if(method == "wavelet"){
  #    log_d <- log_derivative_wavelet(t, s)
  #  }
  else {
    stop("ERROR: Unknown derivative type. Please check and try again")
  }
  return(log_d)
}
#' @title
#' log_derivative_bourdet
#' @description
#' Function to calculate the derivative of the drawdown with respect to the derivative
#' of log time using the approach proposed by Bourdet
#' @param t Numeric vector with the time
#' @param s Numeric vector with the drawdown
#' @param d Numeric value
#' @param return.pos Logical flag to return only the positive values of the log-derivative (default = TRUE)
#' @param log Logical flag to indicate that natural logarithm (a log to the base e) is used in the derivative 
#' calculation (default = TRUE). Logarithm to the base 10 is used if FALSE.
#' @return
#' A list with
#' \itemize{
#'   \item x: Numeric vector with the x coordinates where the log-derivative is evaluated
#'   \item y: Numeric vector with the value of the log-derivative
#' }
#' @family log_derivative functions
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @export
#' @references
#'  Bourdet, D.; Whittle, T.; Douglas, A. & Pirard, Y. A new set of type curves
#'  simplifies well test analysis World Oil, 1983.
#' @examples
#' data(boulton)
#' t <- boulton$t
#' s <- boulton$s
#' ptest <- pumping_test('Well1', Q = 0.03, r = 20, t = t, s = s)
#' #
#' dptest.bd <- log_derivative_bourdet(ptest$t, ptest$s, d = 10)
#' plot(t, s, type= "p", log = "xy", ylim = c(1e-3,2))
#' points(dptest.bd$x, dptest.bd$y, col = "red")
log_derivative_bourdet <- function(t, s, d = 2, return.pos = T, log = T){
  t1 <- NULL
  if(log){
    t1 <- log(t)
  } else {
    t1 <- log10(t)
  }
  dx <- t1[2:length(t)]-t1[1:(length(t)-1)]
  dy <- s[2:length(t)]-s[1:(length(t)-1)]
  dx1 <- dx[1:(length(t)-2*d+1)]
  dx2 <- dx[(2*d):length(t)]
  dy1 <- dy[1:(length(s)-2*d+1)]
  dy2 <- dy[(2*d):length(s)]
  xd <- t[d:(length(t)-d)]
  yd<-(((dx2*dy1)/dx1)+((dx1*dy2)/dx2))/(dx1+dx2)
  if(return.pos){
    pos_valid <- !is.na(yd) & yd > 1.0e-10
    xd <- xd[pos_valid]
    yd <- yd[pos_valid]
  }
  pos <- !is.na(yd)
  results <- list(x = xd[pos], y = yd[pos])
  return(results)
}
#' @title
#' log_derivative_horner
#' @description
#' Function to calculate the derivative of the drawdown with respect to the derivative
#' of log time using the approach proposed by Horner
#' @param t Numeric vector with the time
#' @param s Numeric vector with the drawdown
#' @param return.pos Logical flag to return only the positive values of the log-derivative (default = TRUE)
#' @param log Logical flag to indicate that natural logarithm (a log to the base e) is used in the derivative 
#' calculation (default = TRUE). Logarithm to the base 10 is used if FALSE.
#' @return
#' A list with
#' \itemize{
#'   \item x: Numeric vector with the x coordinates where the log-derivative is evaluated
#'   \item y: Numeric vector with the value of the log-derivative
#' }
#' @family log_derivative functions
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @export
#' @references
#' Horne, R. N. Modern Well Test Analysis: A Computer-Aided Approach Petroway,
#' Incorporated, 1990
#' @examples
#' data(boulton)
#' t <- boulton$t
#' s <- boulton$s
#' ptest <- pumping_test('Well1', Q = 0.03, r = 20, t = t, s = s)
#' #
#' dptest.hn <- log_derivative_horner(ptest$t, ptest$s)
#' plot(t, s, type= "p", log = "xy", ylim = c(1e-3,2))
#' points(dptest.hn$x, dptest.hn$y, col = "red")
log_derivative_horner <- function(t, s, return.pos = T, log = T){
  end <- length(t)
  t1 <- t[1:(end-2)]
  t2 <- t[2:(end-1)]
  t3 <- t[3:end];
  s1 <- s[1:(end-2)];
  s2 <- s[2:(end-1)]
  s3 <- s[3:end];
  d1 <- NULL
  d2 <- NULL
  d3 <- NULL
  if(log){
    d1 <- (log(t2/t1)*s3)/(log(t3/t2)*log(t3/t1))
    d2 <- (log(t3*t1/t2^2)*s2)/(log(t3/t2)*log(t2/t1))
    d3 <- (log(t3/t2)*s1)/(log(t2/t1)*log(t3/t1));
  } else {
    d1 <- (log10(t2/t1)*s3)/(log10(t3/t2)*log10(t3/t1))
    d2 <- (log10(t3*t1/t2^2)*s2)/(log10(t3/t2)*log10(t2/t1))
    d3 <- (log10(t3/t2)*s1)/(log10(t2/t1)*log10(t3/t1));
  }
  yd <- d1+d2-d3;
  xd <- t2;
  if(return.pos){
    pos_valid <- !is.na(yd) & yd > 1.0e-10
    xd <- xd[pos_valid]
    yd <- yd[pos_valid]
  }
  pos <- !is.na(yd)
  results <- list(x = xd[pos], y = yd[pos])
  return(results)
}
#' @title
#' log_derivative_central
#' @description
#' Function to calculate the derivative of the drawdown with respect to the derivative
#' of log time using the central finite differences
#' @param t Numeric vector with the time
#' @param s Numeric vector with the drawdown
#' @param return.pos Logical flag to return only the positive values of the log-derivative (default = TRUE)
#' @param log Logical flag to indicate that natural logarithm (a log to the base e) is used in the derivative 
#' calculation (default = TRUE). Logarithm to the base 10 is used if FALSE.
#' @return
#' A list with
#' \itemize{
#'   \item x: Numeric vector with the x coordinates where the log-derivative is evaluated
#'   \item y: Numeric vector with the value of the log-derivative
#' }
#' @family log_derivative functions
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @export
#' @examples
#' data(boulton)
#' t <- boulton$t
#' s <- boulton$s
#' ptest <- pumping_test('Well1', Q = 0.03, r = 20, t = t, s = s)
#' #
#' dptest.cn <- log_derivative_central(ptest$t, ptest$s)
#' plot(t, s, type = "p", log = "xy", ylim = c(1e-3,2))
#' points(dptest.cn$x, dptest.cn$y, col = "red")
log_derivative_central <- function(t, s, return.pos = T, log = T){
  pos <- t > 0.0
  t <- t[pos]
  s <- s[pos]
  dx <- t[2:length(t)]-t[1:length(t)-1]
  dy <- s[2:length(t)]-s[1:length(t)-1]
  xd <- sqrt(t[1:(length(t)-1)]*t[2:length(t)])
  yd <- xd*dy/dx
  pos.valid <- rep(TRUE, length(t))
  if(return.pos){
    pos_valid <- !is.na(yd) & yd > 1e-10
    xd <- xd[pos_valid]
    yd <- yd[pos_valid]
  }
  results <- list(x = xd, y = yd)
  return(results)
}
#' @title
#' log_derivative_spline
#' @description
#' Function to calculate the derivative of the drawdown with respect to the derivative
#' of log time using the interpolation of the measured points with spline
#' @param t Numeric vector with the time
#' @param s Numeric vector with the drawdown
#' @param n Number of points where the derivative is calculated
#' @param return.pos Logical flag to return only the positive values of the log-derivative (default = TRUE)
#' @param log Logical flag to indicate that natural logarithm (a log to the base e) is used in the derivative 
#' calculation (default = TRUE). Logarithm to the base 10 is used if FALSE.
#' @return
#' A list with
#' \itemize{
#'   \item x: Numeric vector with the x coordinates where the log-derivative is evaluated
#'   \item y: Numeric vector with the value of the log-derivative
#' }
#' @family log_derivative functions
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @importFrom stats spline
#' @export
#' @references
#' Renard, P.; Glenz, D. , Mejias, M. Understanding diagnostic plots for well-test estimation, Hydrogeology Journal, 2008, 17, 589-600
#' @examples
#' # Load test data
#' data(boulton)
#' t <- boulton$t
#' s <- boulton$s
#' ptest <- pumping_test('Well1', Q = 0.03, r = 20, t = t, s = s)
#' #
#' dptest.sp <- log_derivative_spline(ptest$t, ptest$s, n = 30)
#' plot(t, s, type="p", log="xy", ylim = c(1e-3,2))
#' points(dptest.sp$x, dptest.sp$y, col = "red")
log_derivative_spline <- function(t, s, n = 20, return.pos = T, log = T){
  #print(n)
  pos_valid <- t  > 0
  t <- t[pos_valid]
  s <- s[pos_valid]
  min_t <- min(t)
  max_t <- max(t)
  t1 <- logseq(log10(min_t), log10(max_t), n)
  s1 <- spline(t, s, xout = t1)
  end_p <- length(s1$x)
  x <- s1$x[2:(end_p-1)];
  y <- x*(s1$y[3:end_p]-s1$y[1:(end_p-2)])/(s1$x[3:end_p]-s1$x[1:(end_p-2)])
  if(return.pos){
    pos_valid <- !is.na(y) & y > 1.0e-10
    x <- x[pos_valid]
    y <- y[pos_valid]
  }
  res <- list(x = x, y =y)
  return(res)
}
#' @title
#' log_derivative_spane
#' @description
#' Function to calculate the derivative of the drawdown with respect to the derivative
#' of log time using the approach proposed by Spane(1993) based on linear regression
#' @param t Numeric vector with the time
#' @param s Numeric vector with the drawdown
#' @param n Number of points where the derivative is calculated
#' @param return.pos Logical flag to return only the positive values of the log-derivative (default = TRUE)
#' @param log Logical flag to indicate that natural logarithm (a log to the base e) is used in the derivative 
#' calculation (default = TRUE). Logarithm to the base 10 is used if FALSE.
#' @return
#' A list with
#' \itemize{
#'   \item x: Numeric vector with the x coordinates where the log-derivative is evaluated
#'   \item y: Numeric vector with the value of the log-derivative
#' }
#' @family log_derivative functions
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @importFrom stats lm
#' @export
#' @references
#' Spane, F. & Wurstner, S. DERIV: a computer program for calculating pressure derivatives for use in hydraulic test analysis Ground Water, 1993, 31, 814-824
#' @examples
#' data(boulton)
#' t <- boulton$t
#' s <- boulton$s
#' ptest <- pumping_test('Well1', Q = 0.03, r = 20, t = t, s = s)
#' #
#' dptest.sp <- log_derivative_spane(ptest$t, ptest$s, n = 6)
#' plot(t, s, type = "p", log = "xy", ylim = c(1e-3, 2))
#' points(dptest.sp$x, dptest.sp$y, col="red")
log_derivative_spane <- function(t, s, n = 2, return.pos = T, log = T){
  pos <- t > 0.0
  t <- t[pos]
  s <- s[pos]
  x <- t
  y <- s
  nd <- length(x)
  xd <- vector("numeric", length = (nd-n))
  yd <- vector('numeric', length = (nd-n))
  pos <- 1
  for(i in (n+1):(nd-n-1)){
    beginp <- i
    if(log){
      xc1 <- log(x[(beginp-n):beginp])
      xc2 <- log(x[beginp:(beginp+n)])
    } else {
      xc1 <- log10(x[(beginp-n):beginp])
      xc2 <- log10(x[beginp:(beginp+n)])
    }
    yc1 <- y[(beginp-n):beginp]
    yc2 <- y[beginp:(beginp+n)]
    current_data1 <- data.frame(cbind(xc1, yc1))
    current_data2 <- data.frame(cbind(xc2, yc2))
    current.lm1 <- lm(yc1 ~ xc1, data = current_data1)
    current.lm2 <- lm(yc2 ~ xc2, data = current_data2)
    xd[pos] <- x[beginp]
    m1 <- current.lm1$coefficients[2]
    m2 <- current.lm2$coefficients[2]
    dx1 <- max(xc1)-min(xc1)
    dx2 <- max(xc2)-min(xc2)
    yd[pos] <-  (m1*dx2+m2*dx1)/(2*(dx1+dx2))
    pos <- pos + 1
  }
  if(return.pos){
    pos_valid <- !is.na(yd) & yd > 1e-10
    xd <- xd[pos_valid]
    yd <- yd[pos_valid]
  }
  pos <- !is.na(yd) & yd > 1e-10
  results <- list(x = xd[pos], y = yd[pos])
  return(results)
}
#' @title
#' log_derivative_smoothspline
#' @description
#' Function to calculate the derivative of the drawdown with respect to the derivative
#' of log time fitting a smoothing spline to the measured data (Generalized Cross
#' Validation),
#' @param t Numeric vector with the time
#' @param s Numeric vector with the drawdown
#' @param return.pos Logical flag to return only the positive values of the log-derivative (default = TRUE)
#' @param log Logical flag to indicate that natural logarithm (a log to the base e) is used in the derivative 
#' calculation (default = TRUE). Logarithm to the base 10 is used if FALSE.
#' @return
#' A list with
#' \itemize{
#'   \item x: Numeric vector with the x coordinates where the log-derivative is evaluated
#'   \item y: Numeric vector with the value of the log-derivative
#' }
#' @family log_derivative functions
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @importFrom stats smooth.spline
#' @importFrom stats predict
#' @export
#' @examples
#' data(boulton)
#' t <- boulton$t
#' s <- boulton$s
#' ptest <- pumping_test('Well1', Q = 0.03, r = 20, t = t, s = s)
#' #
#' dptest.sm <- log_derivative_smoothspline(ptest$t, ptest$s)
#' plot(t,s,type="p",log="xy",ylim=c(1e-3,2))
#' lines(dptest.sm$x,dptest.sm$y,col="red")
log_derivative_smoothspline <- function(t, s, return.pos = T, log = T){
  pos <- t > 0.0
  t <- t[pos]
  s <- s[pos]
  x <- t
  y <- s
  logx <- log(x)
  if(!log){
    logx <- log10(x)
  }
  res <- smooth.spline(x = logx, y = y, cv = FALSE)
  t1 <- log10(t)
  t1a <- logseq(from = 1.01*min(t1), to = 0.99*max(t1), n = length(t1))
  logt1a <- log(t1a)
  if(!log){
    logt1a <- log10(t1a)
  }
  res1 <- predict(res, logt1a, deriv = 1)
  y <- res1$y
  if(return.pos){
    if(log){
      x <- exp(res1$x)
      y <- abs(res1$y)
    } else{
      x <- 10**(res1$x)
      y <- abs(res1$y)/log(10)
    }
  }
  results <- list(x = x, y = y, res = res)
  return(results)
}
#' @title
#' log_derivative_kernelreg
#' @description
#' Function to calculate the derivative of the drawdown with respect to the derivative
#' of log time using kernel regression of the measured data.
#' @param t Numeric vector with the time
#' @param s Numeric vector with the drawdown
#' @param bw bandwidth
#' @param return.pos Logical flag to return only the positive values of the log-derivative (default = TRUE)
#' @param log Logical flag to indicate that natural logarithm (a log to the base e) is used in the derivative 
#' calculation (default = TRUE). Logarithm to the base 10 is used if FALSE.
#' @return
#' A list with
#' \itemize{
#'   \item x: Numeric vector with the x coordinates where the log-derivative is evaluated
#'   \item y: Numeric vector with the value of the log-derivative
#'   \item res: Kernel regression object
#' }
#' @family log_derivative functions
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @importFrom KernSmooth locpoly
#' @importFrom KernSmooth dpill
#' @export
#' @examples
#' data(boulton)
#' t <- boulton$t
#' s <- boulton$s
#' ptest <- pumping_test('Well1', Q = 0.03, r = 20, t = t, s = s)
#' #
#' dptest.ks <- log_derivative_kernelreg(ptest$t, ptest$s)
#' plot(t,s,type="p",log="xy",ylim=c(1e-3,2))
#' lines(dptest.ks$x,dptest.ks$y,col="red")
log_derivative_kernelreg <- function(t, s, bw = NULL, return.pos = T, log = T){
  pos <- t > 0.0
  t <- t[pos]
  s <- s[pos]
  ntvalues <- length(t)
  incr <- 0
  if(is.null(bw)){
    bw <- NaN
    while(is.nan(bw)){
      if(log){
        current.t <- log(t+incr)
      } else {
        current.t <- log10(t+incr)
      }
      bw <- dpill(x = current.t, y = s, blockmax = 5, divisor = 20,
                  trim = 0.01,
                  proptrun = 0.05,
                  gridsize = ntvalues)
      incr <- incr + 1
    }
  }
  #
  logt <- log(t)
  if(!log){
    logt <- log10(t)
  }
  res <- locpoly(x = logt, y = s, drv = 1L, degree = 2, kernel = "normal",
                 bandwidth = bw,
                 gridsize = ntvalues)
  pos.valid <- seq(1, ntvalues, by = 1)
  if(return.pos){
    pos.valid <- res$y > 0 & !is.na(res$y)
    if(log){
      x <- exp(res$x[pos.valid])
      y <- res$y[pos.valid]/log(10)
    } else {
      x <- 10**(res$x[pos.valid])
      y <- res$y[pos.valid]
    }
  }
  result <- list(x = x, y = y, res = res)
  return(result)
}
#' @title
#' log_derivative_lokern
#' @description
#' Function to calculate the derivative of the drawdown with respect to the derivative
#' of log time using local kernel regression of the measured data.
#' @param t Numeric vector with the time
#' @param s Numeric vector with the drawdown
#' @param return.pos Logical flag to return only the positive values of the log-derivative (default = TRUE)
#' @param log Logical flag to indicate that natural logarithm (a log to the base e) is used in the derivative 
#' calculation (default = TRUE). Logarithm to the base 10 is used if FALSE.
#' @return
#' A list with
#' \itemize{
#'   \item x: Numeric vector with the x coordinates where the log-derivative is evaluated
#'   \item y: Numeric vector with the value of the log-derivative
#'   \item res: Local Kernel Regression object
#' }
#' @family log_derivative functions
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @importFrom lokern lokerns
#' @export
#' @examples
#' data(boulton)
#' t <- boulton$t
#' s <- boulton$s
#' ptest <- pumping_test('Well1', Q = 0.03, r = 20, t = t, s = s)
#' #
#' dptest.lp <- log_derivative_lokern(ptest$t, ptest$s)
#' plot(t, s, type = "p", log = "xy", ylim = c(1e-3,2))
#' points(dptest.lp$x, dptest.lp$y, col = "red")
log_derivative_lokern <- function(t, s, return.pos = T, log = T){
  pos <- t > 0.0
  t <- t[pos]
  s <- s[pos]
  ntdat <- length(t)
  current.t <- log(t)
  if(!log){
    current.t <- log10(t)
  }
  res <- lokerns(current.t, s, deriv = 1, n.out = length(t), 
                 x.out = current.t)
  pos <- rep(T, ntdat)
  if(return.pos){
    pos <- res$est > 1.0e-5
  } 
  if(log){
    x <- exp(res$x.out[pos])
    y = res$est[pos] 
  }
  else {
    x <- 10**(res$x.out[pos])
    y = res$est[pos] 
  }
  results <- list(x = x, y = y, res = res)
  return(results)
}
#' @title
#' log_derivative_locpol
#' @description
#' Function to calculate the derivative of the drawdown with respect to the derivative
#' of log time using local polynomial regression of the measured data.
#' @param t Numeric vector with the time
#' @param s Numeric vector with the drawdown
#' @param return.pos Logical flag to return only the positive values of the log-derivative (default = TRUE)
#' @param log Logical flag to indicate that natural logarithm (a log to the base e) is used in the derivative 
#' calculation (default = TRUE). Logarithm to the base 10 is used if FALSE.
#' @return
#' A list with
#' \itemize{
#'   \item x: Numeric vector with the x coordinates where the log-derivative is evaluated
#'   \item y: Numeric vector with the value of the log-derivative
#'   \item res: Local Polinomial Regression object
#' }
#' @family log_derivative functions
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @importFrom locpol regCVBwSelC
#' @importFrom locpol locPolSmootherC
#' @importFrom locpol gaussK
#' @export 
#' @examples
#' data(boulton)
#' t <- boulton$t
#' s <- boulton$s
#' ptest <- pumping_test('Well1', Q = 0.03, r = 20, t = t, s = s)
#' #
#' dptest.lp <- log_derivative_locpol(ptest$t, ptest$s)
#' plot(t, s, type = "p", log = "xy", ylim = c(1e-3,2))
#' points(dptest.lp$x, dptest.lp$y, col = "red")
log_derivative_locpol <- function(t, s, return.pos = T, log = T){
  pos <- t > 0.0
  t <- t[pos]
  s <- s[pos]
  teval <- t
  logt <- log(t)
  logteval <- log(t)
  if(!log){
    logt <- log10(t)
    logteval <- log10(t)
  }
  bw <- regCVBwSelC(logt, s, deg = 2, kernel = gaussK,
                    weig = rep(1, length(t)))
  res <- locPolSmootherC(logt, s, logteval, bw = bw,
                         deg = 2, kernel = gaussK)
  if(log){
    x <- exp(res$x)
  } else {
    x <- 10**(res$x)
  }
  results <- list(x = x, y = res$beta1, res = res)
  return(results)
}
#' @title
#' log_derivative_lpridge
#' @description
#' Function to calculate the derivative of the drawdown with respect to the derivative
#' of log time using local ridge regression of the measured data.
#' @param t Numeric vector with the time
#' @param s Numeric vector with the drawdown
#' @param return.pos Logical flag to return only the positive values of the log-derivative (default = TRUE)
#' @param log Logical flag to indicate that natural logarithm (a log to the base e) is used in the derivative 
#' calculation (default = TRUE). Logarithm to the base 10 is used if FALSE.
#' @return
#' A list with
#' \itemize{
#'   \item x: Numeric vector with the x coordinates where the log-derivative is evaluated
#'   \item y: Numeric vector with the value of the log-derivative
#'   \item res: A LPRIDGE Regression object
#' }
#' @family log_derivative functions
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @importFrom lpridge lpridge
#' @export
#' @examples
#' data(boulton)
#' t <- boulton$t
#' s <- boulton$s
#' ptest <- pumping_test('Well1', Q = 0.03, r = 20, t = t, s = s)
#' #
#' dptest.lp <- log_derivative_lpridge(ptest$t, ptest$s)
#' plot(t, s, type = "p", log = "xy", ylim = c(1e-3,2))
#' points(dptest.lp$x, dptest.lp$y, col = "red")
log_derivative_lpridge <- function(t, s, return.pos = T, log = T){
  pos <- t > 0.0
  t <- t[pos]
  s <- s[pos]
  #
  logt <- log(t)
  if(!log){
    logt <- log10(t)
  }
  bw <- bandwidth_lpridge(logt, s)
  p.min <- which.min(bw$cverror)
  h <- bw$h[p.min]
  res <- lpridge(logt, s, bandwidth = h, deriv = 1, x.out = logt)
  x <- exp(res$x)
  if(!log){
    x <- 10**(res$x)
  }
  results <- list(x = x, y = res$est, res = res)
}
