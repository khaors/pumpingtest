#' @section Objective Function functions:
#' residual_sum_squares, loglikelihood, bc_transf, loglikelihood_bc, mean_absolute_deviation,
#' max_absolute_deviation, mean_square_error
#' @docType package
#' @name pumpingtest
NULL
#' @title
#' residual_sum_squares
#' @description
#' Function to calculate the residual sum of squares between the measured and simulated drawdown.
#' @param par A list with the parameters of the model.
#' @param ptest A pumping_test object.
#' @param model A character string with the name of the model
#' @return
#' The numeric value of the residual sum of squares
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @family objective_function functions
#' @export
residual_sum_squares <- function(par, ptest, model){
  if(class(ptest) != 'pumping_test'){
    stop('A pumping_test object is required as input')
  }
  smodel <- NULL
  solution_fn <- paste('smodel<-', model,'_solution(ptest,', sep = '')
  npar <- length(par)
  if(npar == 0){
    stop('Incorret number of parameters')
  }
  for(i in 1:npar){
    solution_fn <- paste(solution_fn, 'par[', as.character(i) , '],', sep = '')
  }
  solution_fn <- paste(solution_fn, 'ptest$t)', sep = '')
  eval(parse(text = solution_fn))
  res <- sum((smodel-ptest$s)^2)
  return(res)
}
#' @title
#' loglikelihood
#' @description
#' Function to calculate the loglikelihood function of the pumping test data assuming that the
#' residuals between the measured and simulated drawdown follow a normal distribution.
#' @param par A list with the parameters of the model.
#' @param ptest A pumping_test object.
#' @param model A character string with the name of the model.
#' @importFrom stats dnorm
#' @return
#' The numeric value of the loglikelihood function.
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @family objective_function functions
#' @export
loglikelihood <- function(par, ptest, model){
  if(class(ptest) != 'pumping_test'){
    stop('A pumping_test object is required as input')
  }
  smodel <- NULL
  solution_fn <- paste('smodel<-', model,'_solution(ptest,', sep = '')
  npar <- length(par)
  # Correct for negative parameters
  pos_invalid <- par < 0
  par[pos_invalid] <- 1e-3
  if(npar == 0){
    stop('Incorret number of parameters')
  }
  for(i in 1:(npar-1)){
    solution_fn <- paste(solution_fn, 'par[', as.character(i) , '],', sep = '')
  }
  solution_fn <- paste(solution_fn, 'ptest$t)', sep = '')
  #print(solution_fn)
  eval(parse(text = solution_fn))
  # Likelihood function for a normal distribution
  loglik <- dnorm(ptest$s, mean = smodel, sd = par[npar], log = TRUE)
  pos_invalid <- is.na(loglik)
  if(sum(pos_invalid)>=1){
    print(smodel)
    print(par)
    print(loglik)
    stop('NA')
  }
  res <- sum(loglik)
  return(res)
}
#' @title
#' box_cox_transform
#' @description
#' Function to apply the Box-Cox transform to a vector.
#' @param y A numeric vector with the values to be transformed.
#' @param lambda A numeric value with the lambda value that defines the transformation.
#' @return
#' A numeric vector with the transformed values
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @family objective_function functions
#' @export
box_cox_transf <- function(y, lambda){
  y1 <- (y^lambda-1)/lambda
  return(y1)
}
#' @title
#' loglikelihood_bc
#' @description
#' Function to calculate the loglikelihood function of the pumping test data using the Box-Cox
#' transform assuming that the residuals between the measured and simulated drawdown follow
#' a normal distribution.
#' @param par A list with the parameters of the model.
#' @param ptest A pumping_test object.
#' @param model A character string with the name of the model.
#' @return
#' The numeric value of the loglikelihood function.
#' @importFrom stats dnorm var
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @family objective_function functions
#' @export
loglikelihood_bc <- function(par, ptest, model){
  if(class(ptest) != 'pumping_test'){
    stop('A pumping_test object is required as input')
  }
  smodel <- NULL
  solution_fn <- paste('smodel<-', model,'_solution(ptest,', sep = '')
  npar <- length(par)
  if(npar == 0){
    stop('Incorret number of parameters')
  }
  for(i in 1:(npar-1)){
    solution_fn <- paste(solution_fn, 'par[', as.character(i) , '],', sep = '')
  }
  solution_fn <- paste(solution_fn, 'ptest$t)', sep = '')
  #print(solution_fn)
  eval(parse(text = solution_fn))
  # find lambda of box-cox transform
  ve <- vector('numeric', length = 30)
  dlambda <- 1.0/(length(ve)+1)
  for(i in 1:length(ve)){
    clambda <- dlambda*i
    obs <- box_cox_transf(ptest$s, clambda)
    sim <- box_cox_transf(smodel, clambda)
    e <- (obs-sim)/clambda
    ve[i] <- var(e)
  }
  #print(ve)
  p <- which.min(ve)
  #print(ve[p])
  clambda <- dlambda*p
  obs <- box_cox_transf(ptest$s, clambda)
  sim <- box_cox_transf(smodel, clambda)
  #loglik <- dnorm(ptest$s, mean = smodel, sd = par[npar], log = TRUE)
  loglik <- dnorm(obs, mean = sim, sd = par[npar], log = TRUE)
  res <- sum(loglik)+(clambda-1)*sum(log(ptest$s))
  return(res)
}
#' @title
#' mean_absolute_deviation
#' @description
#' Function to calculate the mean absolute deviation between the measured and simulated drawdown.
#' @param par A list with the parameters of the model.
#' @param ptest A pumping_test object.
#' @param model A character string with the name of the model.
#' @return
#' The numeric value of the mean absolute deviation.
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @family objective_function functions
#' @export
mean_absolute_deviation <- function(par, ptest, model){
  if(class(ptest) != 'pumping_test'){
    stop('A pumping_test object is required as input')
  }
  smodel <- NULL
  solution_fn <- paste('smodel<-', model,'_solution(ptest,', sep = '')
  npar <- length(par)
  if(npar == 0){
    stop('Incorret number of parameters')
  }
  for(i in 1:npar){
    solution_fn <- paste(solution_fn, 'par[', as.character(i) , '],', sep = '')
  }
  solution_fn <- paste(solution_fn, 'ptest$t)', sep = '')
  #print(solution_fn)
  eval(parse(text = solution_fn))
  res <- mean(abs(smodel-ptest$s))
  return(res)
}
#' @title
#' max_absolute_deviation
#' @description
#' Function to calculate the max absolute deviation between the measured and simulated drawdown.
#' @param par A list with the parameters of the model.
#' @param ptest A pumping_test object.
#' @param model A character string with the name of the model.
#' @return
#' The numeric value of the max absolute deviation.
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @family objective_function functions
#' @export
max_absolute_deviation <- function(par, ptest, model){
  if(class(ptest) != 'pumping_test'){
    stop('A pumping_test object is required as input')
  }
  smodel <- NULL
  solution_fn <- paste('smodel<-', model,'_solution(ptest,', sep = '')
  npar <- length(par)
  if(npar == 0){
    stop('Incorret number of parameters')
  }
  for(i in 1:npar){
    solution_fn <- paste(solution_fn, 'par[', as.character(i) , '],', sep = '')
  }
  solution_fn <- paste(solution_fn, 'ptest$t)', sep = '')
  #print(solution_fn)
  eval(parse(text = solution_fn))
  res <- max(abs(smodel-ptest$s))
  return(res)
}
#' @title
#' mean_square_error
#' @description
#' This function calculates the mean square error between the measured and calculated
#' drawdown
#' @param par A list with the parameters of the model.
#' @param ptest A pumping_test object.
#' @param model A character string with the name of the model.
#' @return
#' The numeric value of the mean squared error
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @family objective_function functions
#' @export
#' @examples
#' data(theis)
mean_square_error <- function(par, ptest, model) {
  if(class(ptest) != 'pumping_test'){
    stop('A pumping_test object is required as input')
  }
  smodel <- NULL
  solution_fn <- paste('smodel<-', model,'_solution(ptest,', sep = '')
  npar <- length(par)
  if(npar == 0){
    stop('Incorret number of parameters')
  }
  for(i in 1:npar){
    solution_fn <- paste(solution_fn, 'par[', as.character(i) , '],', sep = '')
  }
  solution_fn <- paste(solution_fn, 'ptest$t)', sep = '')
  #print(solution_fn)
  eval(parse(text = solution_fn))
  res <- sqrt(mean((smodel-ptest$s)^2))
  return(res)
}
