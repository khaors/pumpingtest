#' @section twalk Auxiliary functions:
#' psi1, traverse_kernel, walk_kernel, hop_kernel, blow_kernel, log_ratio_density_blow, log_ratio_density_hop, run_twalk_MCMC
#' @docType package
#' @name pumpingtest
NULL
#' @title
#' psi1
#' @description
#' return a random number with Psi distribution
#' @return
#' A numeric value
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @export
#' @references
#' Christen, J. A. & Fox, C. A general purpose sampling algorithm for continuous distributions (the t-walk).
#' Bayesian Analysis, 2010, 5, 2, 263-281.
#' @family twalk_auxiliary_function functions
psi1 <- function(){
  at = 6; # PARAMETER_at in distribution
  rand <- runif(1)
  if (rand < (at-1)/(2*at)){
    beta <- rand^(1/(1+at))
  }
  else{
    beta <- rand^(1/(1-at));
  }
  return(beta)
}
#' @title
#' traverse_kernel
#' @description
#' PDF of the traverse move
#' @param x Value
#' @param xp Value
#' @param phi Logical vector with the parameters to change.
#' @return
#' A numeric value
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @export
#' @references
#' Christen, J. A. & Fox, C. A general purpose sampling algorithm for continuous distributions (the t-walk).
#' Bayesian Analysis, 2010, 5, 2, 263-281.
#' @family twalk_auxiliary_function functions
traverse_kernel <- function (x, xp, phi){
  # simulate traverse kernel
  beta <- psi1() # simulate psi1 to give beta
  h <- x;
  h[phi] <- xp[phi] + beta*(xp[phi] - x[phi])
  nI <- sum(phi);
  lgr <- (nI-2)*log(beta); # since psi1(beta) = psi1(1/beta)
  res <- list(h = h, lgr = lgr)
}
#' @title
#' walk_kernel
#' @description
#' Walk step
#' @param x Numeric value
#' @param xp Numeric value
#' @param phi Logical vector with the parameters to change.
#' @return
#' Numeric value
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @export
#' @references
#' Christen, J. A. & Fox, C. A general purpose sampling algorithm for continuous distributions (the t-walk).
#' Bayesian Analysis, 2010, 5, 2, 263-281.
#' @family twalk_auxiliary_function functions
walk_kernel <- function(x, xp, phi){
  #[h,lgr] = Simh2(x,xp)
  #walk kernel
  aw <- 1.5; # PARAMETER_aw in distribution
  u <- runif(length(x))
  z <- (aw/(1+aw))*(-1 + 2*u + aw*u^2) # z ~ 1/sqrt(1+z) in [-a/(1+a),a]
  h <-  x + phi*(x - xp)*z;
  lgr <- 0; # As z has density 1/sqrt(1+z)
  res <- list(h = h, lgr = lgr)
  return(res)
}
#' @title
#' hop_kernel
#' @description
#' Hop kernel
#' @param x Numeric value
#' @param xp Numeric value
#' @param phi Logical vector with the parameters to change.
#' @return
#' Numeric value
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @export
#' @references
#' Christen, J. A. & Fox, C. A general purpose sampling algorithm for continuous distributions (the t-walk).
#' Bayesian Analysis, 2010, 5, 2, 263-281.
#' @family twalk_auxiliary_function functions
hop_kernel <- function(x, xp, phi){
  # hop kernel
  if (sum(phi)>0){
    sigma <- max(abs(xp[phi] - x[phi]))
  }
  else{
    sigma <- 0
  }
  h <- x + phi*sigma/3.*rnorm(length(x)); # iid standard normal
  #evaluate log of ratio of densities g(x|h,x')/g(h|x,x')
  lgr <- log_ratio_density_hop(x,h,xp,phi) - log_ratio_density_hop(h,x,xp,phi)
  return(lgr)
}
#' @title
#' log_ratio_density_hop
#' @description
#' Function to evaluate the proposal distribution of the hop step
#' @param h Numeric value
#' @param x Numeric value
#' @param xp Numeric value
#' @param phi Logical vector with the parameters to change.
#' @return
#' Numeric value
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @export
#' @references
#' Christen, J. A. & Fox, C. A general purpose sampling algorithm for continuous distributions (the t-walk).
#' Bayesian Analysis, 2010, 5, 2, 263-281.
#' @family twalk_auxiliary_function functions
log_ratio_density_hop <- function(h, x, xp, phi){
  # log density for hop move, up to constant not affecting ratio
  if(sum(phi)>0){
    sigma <- max(abs(xp[phi] - x[phi]))
  }
  else{
    sigma <- 0
  }
  nI <- sum(phi);
  lg <- -nI*log(sigma) - sum((h-x)^2)/(2*sigma^2)
  return(lg)
}
#' @title
#' blow_kernel
#' @description
#' Function to calculate blow step
#' @param x Numeric value
#' @param xp Numeric value
#' @param phi Logical vector with the parameters to change.
#' @return
#' Numeric value
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @export
#' @references
#' Christen, J. A. & Fox, C. A general purpose sampling algorithm for continuous distributions (the t-walk).
#' Bayesian Analysis, 2010, 5, 2, 263-281.
#' @family twalk_auxiliary_function functions
blow_kernel <- function(x, xp, phi){
  # [h,lgr] = Simh4(x,xp)
  # blow kernel
  if (sum(phi)>0){
    sigma <- max(abs(xp[phi] - x[phi]))
  }
  else{
    sigma <- 0
  }
  h <- x
  h[phi] <- xp[phi] + sigma*rnorm(length(x[phi])); # iid standard normal
  
  # evaluate log of ratio of densities g(x|h,x')/g(h|x,x')
  lgr <- log_ratio_density_blow(x,h,xp,phi) - log_ratio_density_blow(h,x,xp,phi);
  return(lgr)
}
#' @title
#' log_ratio_density_blow
#' @description
#' Function to evaluate the proposal distribution of the hop step
#' @param h Numeric value
#' @param x Numeric value
#' @param xp Numeric value
#' @param phi Logical vector with the parameters to change.
#' @return
#' Numeric value
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @export
#' @references
#' Christen, J. A. & Fox, C. A general purpose sampling algorithm for continuous distributions (the t-walk).
#' Bayesian Analysis, 2010, 5, 2, 263-281.
#' @family twalk_auxiliary_function functions
log_ratio_density_blow <- function(h, x, xp, phi){
  # log density for blow move, up to constant not affecting ratio
  if (sum(phi)>0){
    sigma <- max(abs(xp[phi] - x[phi]))
  }
  else{
    sigma <- 0
  }
  nI = sum(phi);
  lg = -nI*log(sigma) - sum((h-xp)^2)/(2*sigma^2)
  return(lg)
}
#' @title
#' run_twalk_MCMC
#' @description
#' Function that implements the transverse walk sampling strategy to obtain realizations of the
#' posterior distribution of the pumping test parameters.
#' @param startvalue1 Numeric vector with the values of the hydraulic parameters used to initialize the chain.
#' @param startvalue2 Numeric vector with the values of the hydraulic parameters used to initialize the chain.
#' It must be different than the vector startvalue1.
#' @param ptest A pumping test object
#' @param model A character string specifying the model that describe the flow during the pumping test
#' @param prior.pdf A character vector with the names of the PDF of the hydraulic parameters.
#' Currently only uniform ('unif') and normal ('norm') distributions are supported.
#' @param prior.parameters A numeric vector with the parameters of the PDF of the hydraulic parameters.
#' The min and max values of the hydraulic parameters are specified in the case of a uniform distribution.
#' The mean and sigma values of the hydraulic parameters are specified in the case of a normal distribution.
#' @param iterations Number of iterations to run the chain.
#' @return
#' This function returns a list with the following entries:
#' \itemize{
#' \item xxp = xxp
#' \item logdensity: matrix with the values of the logdensity of the accepted parameters.
#' \item acceptance: the acceptance rate
#' }
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @importFrom pracma mod eps
#' @export
#' @references
#' Christen, J. A. & Fox, C. A general purpose sampling algorithm for continuous distributions (the t-walk).
#' Bayesian Analysis, 2010, 5, 2, 263-281.
#' @family twalk_auxiliary_function functions
run_twalk_MCMC <- function(startvalue1, startvalue2, ptest, model, prior.pdf,
                           prior.parameters, iterations = 10000){
  # Initial values
  x <- startvalue1
  xp <- startvalue2
  #
  #cat('x= ', x, '\n')
  #cat('xp= ', xp, '\n')
  #print(prior.parameters)
  #
  ltx <- posterior(startvalue1, ptest, model, prior.pdf, prior.parameters)
  ltxp <- posterior(startvalue2, ptest, model, prior.pdf, prior.parameters)
  if(is.infinite(ltx) | is.infinite(ltxp)){
    stop('ERROR: the posterior density at starting values is infinite')
  }
  #cat('ltx= ', ltx, '\n')
  #cat('ltxp= ', ltxp, '\n')
  #
  MoveRatio <- c(0, 60, 60, 1, 1) # first ratio for move 0 -- i.e. don't move.
  endp <- length(MoveRatio)
  MoveProb <- cumsum(MoveRatio[1:endp-1]/sum(MoveRatio))
  n1 <- 4 # expected number of parameters to be moved
  pphi <- min(length(x),n1)/length(x)
  #
  prop <- vector('numeric', length(MoveRatio))       #collect number of proposals
  acc <- vector('numeric', length(MoveRatio))        #collect acceptance ratios
  xxp <- matrix(0.0, (iterations+1), 2*length(x))    #to store output
  lt <- matrix(0.0, (iterations+1), 2)               #store values of log target density
  
  xxp[1,] <- cbind(x, xp)
  lt[1,] <- cbind(ltx, ltxp)
  lalpha <- 0
  y <- 0
  yp <- 0
  lty <- 0
  ltyp <- 0
  #
  for(iter in 1:iterations){
    #if(mod(iter, 10000) == 0){
    #  cat("iteration = ", iter , "\n")
    #}
    rand <- runif(1)
    dir <- rand < 0.5                       # choose a direction
    #cat("rand =", rand, '\n')
    #cat('dir = ', dir, '\n')
    kernum <- sum(MoveProb <= rand);        # choose a kernel (0 to #Moves-1)
    #cat('kernum= ', kernum, '\n')
    phi <- runif(length(x)) < pphi;         # set parameters to move
    #cat('phi= ', phi, '\n')
    prop[kernum+1] = prop[kernum+1] + 1;    #
    if(kernum == 0){
      # do nothing, ensures aperiodicity ... not needed (see paper)
      if(dir) yp <- xp
      else y <- x
      lgr <- 0 # always accepted
    }
    if(kernum == 1){
      # traverse move
      if (dir) {
        res <- traverse_kernel(xp, x, phi)
        yp <- res$h
        lgr <- res$lgr
      }
      else {
        #[y,lgr] =
        res <- traverse_kernel(x, xp, phi)
        y <- res$h
        lgr <- res$lgr
      }
    }
    if(kernum == 2){
      # walk move
      if (dir){
        #[yp,lgr] =
        res <- walk_kernel(xp, x, phi)
        yp <- res$h
        lgr <- res$lgr
      }
      else {
        #[y,lgr]
        res <- walk_kernel(x, xp, phi)
        y <- res$h
        lgr <- res$lgr
      }
    }
    if(kernum == 3){
      # hop move
      if (dir) {
        res <- hop_kernel(xp, x, phi)
      }
      else {
        #[y,lgr]
        res <- hop_kernel(x, xp, phi)
      }
    }
    if(kernum == 4){
      # blow move
      if (dir) {
        #[yp,lgr] =
        res <- blow_kernel(xp,x,phi)
      }
      else {
        #[y,lgr]
        res <- blow_kernel(x,xp,phi)
      }
    }
    #print(res)
    #cat('res-lgr', res$lgr, '\n')
    if (dir){
      y <- x
      lty <- ltx
      par <- yp
      args <- list(par = par, ptest = ptest, model = model, prior.pdf = prior.pdf,
                   prior.parameters = prior.parameters)
      ltyp <- do.call(posterior, args = args)
      lalpha <- ltyp - ltxp + lgr
    }
    else {
      yp <- xp
      ltyp <- ltxp
      par <- y
      args <- list(par = par, ptest = ptest, model = model, prior.pdf = prior.pdf,
                   prior.parameters = prior.parameters)
      lty <- do.call(posterior, args = args)
      #cat('yp= ', yp, '\n')
      #cat('ltyp= ', ltyp, '\n')
      #cat('lty= ', lty, '\n')
      lalpha <- lty - ltx + lgr
      if(is.na(lalpha)) lalpha <- 0
    }
    #
    #cat('lalpha=', lalpha, '\n')
    if (rand < exp(lalpha)){
      #accepted
      acc[kernum+1] <-  acc[kernum+1] + sum(phi)/length(x)
      x <- y
      ltx <- lty
      xp <- yp
      ltxp <- ltyp
    }
    xxp[(iter+1),] <- cbind(x, xp)
    lt[(iter+1),] = cbind(ltx, ltxp)
  }#iter
  #calculate acceptance rates and return results
  eps <- eps(1)
  acc <- c(acc/(prop+eps), sum(acc)/iterations)
  res <- list(xxp = xxp, logdensity = lt, acceptance = acc)
  return(res)
}