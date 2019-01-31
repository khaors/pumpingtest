#' @section Smoothing functions:
#' bandwidth_lpridge, smoothing.pumping_test
#' @docType package
#' @name pumpingtest
NULL
#' @title
#' bandwidth_lpridge
#' @description
#' Function to calculate the ridge regression bandwidth
#' @param t Numeric vector
#' @param s Numeric vector
#' @param h.range values
#' @param seed Random seed
#' @return
#' This function returns a list with the following entries:
#' \itemize{
#' \item h: bandwidth values
#' \item cverror: Crossvalidation error
#' }
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @family smoothing functions
#' @importFrom lpridge lpridge
#' @export
bandwidth_lpridge <- function(t, s, h.range = c(1e-2,1), seed = 12345){
  set.seed(seed)
  ntvalues <- length(t)
  ntvalues1 <- round(0.75*ntvalues)
  tmp <- sample(ntvalues, ntvalues)
  est_samples <- tmp[1:ntvalues1]
  test_samples <- tmp[(ntvalues1+1):ntvalues]
  t_est <- t[est_samples]
  s_est <- s[est_samples]
  t_test <- t[test_samples]
  s_test <- s[test_samples]
  h <- seq(h.range[1], h.range[2], by = h.range[1])
  nh <- length(h)
  cverror <- vector('numeric', length = nh)
  for(ih in 1:nh){
    res <- lpridge(log(t_est), s_est, bandwidth = h[ih], x.out = log(t_test))
    current_var <- mean((res$est-s_test)**2)
    cverror[ih] <- current_var
  }
  results <- list(h = h, cverror = cverror)
}
#' @title
#' smoothing.pumping_test
#' @description
#' Function to smooth or filter the drawdown data using different methodologies.
#' @param ptest A pumping_test object
#' @param ntvalues Integer, number of values
#' @param method Character string specifying the method to be used in the smoothing.
#' Current options include:
#' \itemize{
#' \item spline: smooth.spline
#' \item kernsmooth: KernSmooth
#' \item lokern: local kernel regression
#' \item locpoly: local polynomial regression
#' \item lpridge: lpridge
#' \item wavelet: Wavelet thresholding
#' }
#' @importFrom stats smooth.spline
#' @importFrom KernSmooth locpoly
#' @importFrom KernSmooth dpill
#' @importFrom lokern lokerns
#' @importFrom locpol regCVBwSelC
#' @importFrom locpol locPolSmootherC
#' @importFrom locpol gaussK
#' @importFrom lpridge lpridge
#' @importFrom stats approx
#' @importFrom pracma ceil
#' @export
smoothing.pumping_test <- function(ptest, ntvalues = NULL,
                                   method = c('spline', 'kernsmooth', 'lokern',
                                              'locpol', 'lpridge')){
  if(class(ptest) != 'pumping_test'){
    stop('ERROR: a pumping_test object is required as input')
  }
  if(is.null(ntvalues)){
    ntvalues <- length(ptest$t)
  }
  #
  res <- NULL
  results <- list()
  #
  if(method == 'spline'){
    res <- smooth.spline(x = log(ptest$t), y = ptest$s)
    results$t <- exp(res$x)
    results$s <- res$y
    results$residuals <- ptest$s - res$y
    results$method <- method
    results$smooth <- res
  }
  else if(method == 'kernsmooth'){
    bw <- NaN
    incr <- 0
    while(is.nan(bw)){
      bw <- dpill(x = log(ptest$t+incr), y = ptest$s, blockmax = 5, divisor = 20,
                  trim = 0.01,
                  proptrun = 0.05,
                  gridsize = ntvalues)
      incr <- incr + 1
    }
    mnt <- min(log(ptest$t))
    mxt <- max(log(ptest$t))
    res <- locpoly(x = log(ptest$t), y = ptest$s, degree = 2,
                   kernel = "normal",
                   bandwidth = bw,
                   gridsize = ntvalues)
    results$t <- exp(res$x)
    results$s <- res$y
    results$residuals <- NULL #ptest$s - res$y
    results$method <- method
    results$bw <- bw
    results$smooth <- res
  }
  else if(method == 'lokern'){
    res <- lokerns(log(ptest$t), ptest$s, deriv = 0, n.out = ntvalues,
                   x.out = log(ptest$t))
    results$t <- exp(res$x.out)
    results$s <- res$est
    results$method <- method
    results$smooth <- res
    results$residuals <- res$est - ptest$s
  }
  else if(method == 'locpol'){
    bw <- regCVBwSelC(log(ptest$t), ptest$s, deg = 0, kernel = gaussK,
                      weig = rep(1, length(ptest$t)))
    mnt <- min(log10(ptest$t))
    mxt <- max(log10(ptest$t))
    teval <- ptest$t #logseq(mnt, mxt, length(ptest$t))
    res <- locPolSmootherC(log(ptest$t), ptest$s, log(teval), bw = bw,
                           deg = 2, kernel = gaussK)
    results$t <- teval
    results$s <- res$beta0
    results$residuals <- results$s-ptest$s
    results$smooth <- res
  }
  else if(method == 'lpridge'){
    bw <- bandwidth_lpridge()
    res <- lpridge(log(ptest$t), ptest$s, bandwidth = bw,  deriv = 0, x.out = log(ptest$t))
    results$t <- exp(res$x.out)
    results$s <- res$est
    results$method <- method
    results$smooth <- res
  }
  # else if(method == 'wavelet'){
  #   nt <- length(ptest$t)
  #   n2b <- ceil(log(nt)/log(2))
  #   t.out <- logseq(min(log10(ptest$t)), max(log10(0.99*ptest$t)), 2^n2b)
  #   s.approx <- approx(log10(ptest$t), ptest$s,xout = log10(t.out), method='linear')
  #   waveletwmap <- wd(s.approx$y, family="DaubLeAsymm", filter.number=10)
  #   softthreshwmap <- threshold(waveletwmap, type="soft", policy="universal")
  #   hardthreshwmap <- threshold(waveletwmap, type="hard", policy="universal")
  #   s.soft <- wr(softthreshwmap)
  #   s.hard <- wr(hardthreshwmap)
  #   results$t <- t.out
  #   results$s <- s.soft
  #   results$method <- method
  #   results$smooth  <- s.soft
  # }
  else{
    stop('ERROR: the requested smoothing method is not available')
  }
  return(results)
}
