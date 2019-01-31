#' @section Adaptative MCMC Auxiliary functions:
#' prior, posterior, proposalfunction_cov, proposalfunction, run_adap_metropolis_MCMC
#' @docType package
#' @name pumpingtest
NULL
#' @title
#' prior
#' @description
#' Function to calculate the prior distribution of the parameters used in the drawdown calculations
#' with analytical models.
#' @param par A numeric vector with the values of the drawdown model parameters
#' @param prior.pdf A character vector with the names of the probability density functions for
#' each parameters. Currently only 'unif' and 'norm' are supported.
#' @param prior.parameters A matrix with the parameters of the probability density functions. If
#' the pdf is 'unif' then the min and max must be specified, whereas in the case of a normal
#' distribution the mean and standard deviation.
#' @return
#' A numeric value with the sum of the logarithms of the prior distributions
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @family amcmc_auxiliary_function functions
#' @export
prior <- function(par, prior.pdf, prior.parameters){
  npar <- length(par)
  par.unif <- list(x = -1, min=-1,max=-1, log = TRUE)
  par.norm <- list(x = -1, mean = 0, sd = 1.0, log = TRUE)
  prior <- 0.0
  current_prior <- 0.0
  for(ipar in 1:npar){
    current_pdf <- paste('d', prior.pdf[ipar], sep = '')
    if(prior.pdf[ipar] == 'unif'){
      par.unif$x <- par[ipar]
      par.unif$min <- prior.parameters[ipar,1]
      par.unif$max <- prior.parameters[ipar,2]
      current_prior <- do.call(current_pdf, args = par.unif)
    }
    else {
      par.norm$x <- par[ipar]
      par.norm$mean <- prior.parameters[ipar,1]
      par.norm$sd <- prior.parameters[ipar,2]
      current_prior <- do.call(current_pdf, args = par.norm)
    }
    prior <- prior + current_prior
  }
  return(prior)
}
#' @title
#' posterior
#' @description
#' Function to calculate the posterior distribution of the parameters used in the drawdown calculations
#' with analytical models.
#' @param par A numeric vector with the values of the drawdown model parameters.
#' @param ptest A pumping_test object.
#' @param model A character string specifying the model used in the parameter estimation.
#' @param prior.pdf A character vector with the names of the probability density functions for
#' each parameters. Currently only 'unif' and 'norm' are supported.
#' @param prior.parameters A matrix with the parameters of the probability density functions. If
#' the pdf is 'unif' then the min and max must be specified, whereas in the case of a normal
#' distribution the mean and standard deviation.
#' @return
#' A numeric value with the sum of the logarithms of the posterior distributions
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @family amcmc_auxiliary_function functions
#' @export
posterior <- function(par, ptest, model, prior.pdf, prior.parameters){
  res <- loglikelihood(par, ptest, model) + prior(par, prior.pdf, prior.parameters)
  return(res)
}
#' @title
#' proposalfunction_cov
#' @description
#' Function to propose new parameter values from a multivariate normal distribution with a
#' given mean and covariance matrix. This function is used to improve the convergence of a
#' MCMC.
#' @param mn A numeric vector with the mean value of the parameters
#' @param cv The covariance matrix
#' @return
#' A vector with realizations of a multivariate normal distribution with a given mean and
#' covariance
#' @importFrom  mvtnorm rmvnorm
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @family amcmc_auxiliary_function functions
#' @export
proposalfunction_cov <-  function(mn,cv){
  d <- length(mn)
  cov_par <- ((2.308^2)/d)*cv
  res <- rmvnorm(d, mean = mn, sigma = cov_par)
  return(diag(res))
}
#' @title
#' proposalfunction
#' @description
#' Function to propose new parameter values from a multivariate normal distribution with a
#' given mean and variance.
#' @param param A numeric vector with the mean value of the parameters.
#' @param proposal.sigma A numeric vector with the variance of the parameters.
#' @return
#' A vector with realizations of a multivariate normal distribution with a given mean and
#' varianced
#' @importFrom stats rnorm
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @family amcmc_auxiliary_function functions
#' @export
proposalfunction <- function(param, proposal.sigma){
  #print(param)
  #c(.05,2,0.005)
  nd <- length(param)
  return(rnorm(nd, mean = param, sd = proposal.sigma))
}
#' @title
#' run_adap_metropolis_MCMC
#' @description
#' Function to sample for a complex probabibility density function using MCMC with the
#' adaptative Metropolis algorithm proposed by Roberts and Rosenthal(2009).
#' @param startvalue A numeric vector with the fit parameters of the pumping test.
#' @param iterations An integer with the number of iterations to run the chain.
#' @param iter_update_par An integet specifying the number of iterations to update the
#' covariance matrix.
#' @param ptest A pumping_test object.
#' @param model A character string with the name of the model used in the parameter estimation.
#' @param prior.pdf A character vector with the distributions of the fit parameters (
#' 'unif' and 'norm' are currently supported).
#' @param prior.parameters A matrix with the parameters of the distributions (min and max for
#' uniform distributions, mean and sd for normal distributions)
#' @param proposal.sigma A numeric vector with the standard deviations of the proposal distribution.
#' @param cov.corr A logical flag indicating if the covariance matrix must be corrected for positive
#' definiteness.
#' @return
#' A matrix with the sampled values of the fit parameters.
#' @details
#' This function implements the adaptative MCMC proposed by Roberts and Rosenthal (2009), in
#' which the proposal distribution \eqn{Q_{n}(x, \cdot)} is given by:
#' \deqn{Q_{n}(x, \cdot) = \left\{
#' \begin{aligned}
#' &(1-\theta)N(x, (2.38)^{2}\Sigma_{n}/d) + \theta N(x,(0.1)^{2}I_{d}/d), &\Sigma_{n}\text{ is positive definite} \\
#' &N(x,(0.1)^{2}I_{d}), &\Sigma_{n}\text{ is not positive definitive}\\
#' \end{aligned}
#' \right.}
#' where
#' \itemize{
#' \item \eqn{\theta \in (0,1)}: control parameters
#' \item \eqn{N()}: Normal distribution
#' \item \eqn{\Sigma_{n}}: empirical covariance matrix
#' \item \eqn{d}: number of parameters
#' \item \eqn{I_{d}}: identity matrix of size \eqn{d}.
#' }
#' This proposal function is implemented in the function proposalfunction_cov.
#' @importFrom pracma mod
#' @importFrom Matrix nearPD
#' @importFrom matrixcalc is.positive.definite
#' @importFrom stats cov runif
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @references
#'  Roberts, G. O. & Rosenthal, J. S. Examples of adaptive MCMC Journal of Computational
#'  and Graphical Statistics, 2009, 18, 349-367.
#' @family amcmc_auxiliary_function functions
#' @export
run_adap_metropolis_MCMC <- function(startvalue, iterations = 10000,
                                     iter_update_par = 100,  ptest, model,
                                     prior.pdf, prior.parameters, proposal.sigma,
                                     cov.corr = FALSE){
  nd <- length(startvalue)
  mn <- vector('numeric', length=nd)
  mn[1:nd] <- 0.0
  cv <- diag(1.0, nrow = nd, ncol = nd)
  chain <- array(dim = c(iterations+1,nd))
  chain[1,] <- startvalue
  nochange <- 0
  cv_old <- cv
  for (i in 1:iterations){
    #if(pracma::mod(i,100)==0){
    #  cat("iteration = ", i, "\n")
    #}
    if(i<iter_update_par){
      proposal <- proposalfunction(chain[i,], proposal.sigma)
    }
    else{
      mn <- chain[i,]
      if(mod(i,iter_update_par) == 0){
        begin_pos <- i-iter_update_par+1
        end_pos <- i
        current_chain <- chain[begin_pos:end_pos,]
        mn <- apply(current_chain,2,mean)
        cv <- cov(current_chain)
        #print(chain[begin_pos:end_pos,1])
        #print(c(i,100*(1-nochange/iter_update_par)))
        #print(cv)
        if(nochange != iter_update_par & sum(cv) > 0.0){
          if(!is.positive.definite(cv)){
            #print(c(as.character(i),'No PD MATRIX'))
            #print(cv)
            if(cov.corr){
              cv <- nearPD(cv, keepDiag = TRUE, posd.tol = 1e-5, maxit = 300)
              cv <- as.matrix(cv$mat)
            }
            else {
              cv <- cv_old
            }
          }
        } else {
          cv <- cv_old
        }
        nochange <- 0
        cv_old <- cv
      }
      proposal <- proposalfunction_cov(mn,cv)
    }
    #print(c('update',as.character(i)))
    for(ipar in 1:length(proposal)){
      if(prior.pdf[ipar] == 'unif'){
        if(proposal[ipar] < prior.parameters[ipar,1]){
          proposal[ipar] <- prior.parameters[ipar,1]
        }
        if(proposal[ipar] > prior.parameters[ipar,2]){
          proposal[ipar] <- prior.parameters[ipar,2]
        }
      }
      if(prior.pdf[ipar] == 'norm'){
        if(proposal[ipar] < prior.parameters[ipar,1]-3*prior.parameters[ipar,2]){
          proposal[ipar] <- prior.parameters[ipar,1]-3*prior.parameters[ipar,2]
        }
        if(proposal[ipar] > prior.parameters[ipar,1]+3*prior.parameters[ipar,2]){
          proposal[ipar] <- prior.parameters[ipar,1]+3*prior.parameters[ipar,2]
        }
      }
    }
    #proposal[pos_invalid] <- prior.parameters[pos_invalid,1]
    #print(proposal)
    #print(chain[i,])
    pnew <- posterior(proposal, ptest, model, prior.pdf, prior.parameters)
    pold <- posterior(chain[i,], ptest, model, prior.pdf, prior.parameters)
    #print(pnew)
    #print(pold)
    probab <- exp(pnew-pold)
    #print(probab)
    
    if(is.na(probab)){
      probab <- 1.0e-10
    }
    if(probab < 1.0e-15){
      probab <- 1e-3
    }
    
    if (runif(1) < probab){
      chain[i+1,] <- proposal
    } else {
      chain[i+1,] <- chain[i,]
      nochange <- nochange + 1
    }
  }
  return(chain)
}