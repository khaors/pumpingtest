#' @section Equivalent_time functions:
#'
#' These functions are used in the estimation of the aquifer parameters for the recovery and variable rate tests.
#'
#' The functions included in this section are:
#' recovery_equivalent_time, variable_rate_equivalent_time
#'
#' @docType package
#' @name pumpingtest
NULL
#' @title
#' recovery_equivalent_time
#' @description
#' Function to calculate the equivalent time of the recovery tests
#' @param t_recovery Vector with the times since recovery started
#' @param s_recovery Vector with the drawdowns since recovery started
#' @param pump_history Vector with the pump history of the test
#' @return
#' A list with
#' \itemize{
#' \item ta equivalent time
#' \item sa drawdown corresponding to equivalent time
#' }
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @export
#' @family equivalent_time functions
#' @references
#' Agarwal, R. (1980). A New Method To Account For Producing Time Effects When Drawdown
#' Type Curves Are Used To Analyze Pressure Buildup And Other Test Data. SPE Paper 9289.
#' @examples
#' data(agarwal_recovery)
#' pump_history <- 14400
#' agarwal_corr <- recovery_equivalent_time(agarwal_recovery$t,
#'                                          agarwal_recovery$s,
#'                                          pump_history)
recovery_equivalent_time <- function(t_recovery, s_recovery, pump_history){
  if(class(t_recovery) != 'numeric'){
    stop('Numeric vector required for time')
  }
  if(class(s_recovery) != 'numeric'){
    stop('Numeric vector required for drawdown')
  }
  if(class(pump_history) != 'numeric'){
    stop('Numeric information required for pump history')
  }
  n <- length(pump_history);
  if( n == 1 ){
    ta <- pump_history*t_recovery/(pump_history+t_recovery)
    #print(ta)
  }
  else {
    t_pumping <- pump_history[,1]
    q_pumping <- pump_history[,2]
    ta <- t_pumping[(n-1)]/(t_recovery+t_pumping[(n-1)])^(q_pumping[1]/(q_pumping[n-1]-q_pumping[n]))
    for (j in 2:(n-1)){
      ta <- ta*((t_pumping[n-1]-t_pumping[j-1])/(t_recovery+t_pumping[n-1]-t_pumping[j-1]))^
        ((q_pumping[j]-q_pumping[j-1])/(q_pumping[n-1]-q_pumping[n]))
    }
    ta <- ta*t_recovery;
  }
  pos_valid <- !duplicated(ta);
  ta <- ta[pos_valid]
  sa <- s_recovery[pos_valid];
  results <- list(ta = ta, sa = sa)
  return(results)
}
#' @title
#' variable_rate_equivalent_time
#' @description
#' Function to calculate the equivalent time required in the variable-rate tests.
#' @param t Vector with the times since pumping started
#' @param s Vector with the drawdowns since pumping started
#' @param pump_history Vector with the pump history of the test
#' @return
#' A list with
#' \itemize{
#' \item ta equivalent time
#' \item sa drawdown corresponding to equivalent time
#' }
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @export
#' @family equivalent_time functions
#' @references
#' Birsoy, Y.K. and W.K. Summers, 1980. Determination of aquifer parameters from step tests
#' and intermittent pumping, Ground Water, vol. 18, no. 2, pp. 137-146.
variable_rate_equivalent_time <- function(t, s, pump_history){
  if(class(t) != 'numeric'){
    stop('Numeric vector required for time')
  }
  if(class(s) != 'numeric'){
    stop('Numeric vector required for drawdown')
  }
  if(class(pump_history) != 'matrix'){
    stop('Numeric information required for pump history')
  }
  Q <- pump_history
  dQ <- diff(pump_history[,2])
  dQ <- c(pump_history[1,2], dQ)
  nperiod <- nrow(pump_history)
  ntimes <- length(t)
  starting_times <- vector("numeric", length = (nperiod+1))
  starting_times[1] <- 0
  starting_times[2:(nperiod+1)] <- pump_history[,1]
  pump_history1 <- vector('numeric', length = length(t))
  equiv_time <- vector('numeric', length = length(t))
  equiv_time[1:length(t)] <- 1
  for(i in 1:nperiod){
    if(i == 1){
      current_t <- pump_history[i,1]
      current_pos <- t <= current_t
    }
    else{
      ct1 <- pump_history[(i-1),1]
      ct2 <- pump_history[i,1]
      current_pos <- t > ct1 & t <= ct2
    }
    pump_history1[current_pos] <- i
  }
  #
  for(j in 1:ntimes){
    n1<- pump_history1[j]
    for(i in 1:n1){
      equiv_time[j] <- equiv_time[j]*(t[j]-starting_times[i])^(dQ[i]/Q[n1,2])
    }
  }
  #
  pos_valid <- !duplicated(equiv_time)
  equiv_time <- equiv_time[pos_valid]
  s <- s[pos_valid]
  equiv_time.sort <- sort(equiv_time, index.return = TRUE)
  equiv_time <- equiv_time.sort$x
  lQ <- Q[pump_history1,2]
  sa <- s[equiv_time.sort$ix]/lQ[equiv_time.sort$ix]
  results <- list(ta = equiv_time, sa = sa)
  return(results)
}