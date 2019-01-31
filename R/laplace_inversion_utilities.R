#' @section Laplace Inversion functions:
#' stehfest_coeffs, stehfest_inversion
#' @docType package
#' @name pumpingtest
NULL
#' @title
#' stehfest_coeffs
#' @description
#' Function to calculate the coefficients used in the stehfest inversion. These coefficients
#' \eqn{c_{i}}are given by:
#' \deqn{c_{i}=(-1)^{i+n/2} \sum_{k=\left[ \frac{i+1}{2} \right]}^{\min{(i.n/2)}} \frac{k^{n/2}(2k)!}{\left(\frac{n}{2}-k\right)!k!(k-1)!(i-k)!(2k-i)!}}
#'@param n Number of coefficients
#'@return
#' Numeric vector with the coefficients
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @family laplace functions
#' @export
#' @examples
#' coeffs8 <- stehfest_coeffs(8)
#' coeffs20 <- stehfest_coeffs(20)
stehfest_coeffs <- function(n) {
  coeffs <- vector("numeric", length = n)
  nhlf <- n/2
  for (i in 1:n){
    coeffs[i] <- 0.0
    k1 <- floor((i+1)/2)
    k2 <- min(i,nhlf)
    for(k in k1:k2){
      num <- (k**nhlf)*factorial(2*k)
      den <- factorial(nhlf-k)*factorial(k)*factorial(k-1)*factorial(i-k)*factorial(2*k-i)
      coeffs[i] <- coeffs[i]+num/den
    }
    coeffs[i] <- ((-1)**(i+nhlf))*coeffs[i]
  }
  return(coeffs)
}
#' @title
#' stehfest_inversion
#' @description
#' Function to calculate the Laplace inverse transform using the
#' stehfest algorithm
#' @param t Numeric vector with the time
#' @param coeffs Numeric vector with the coefficients
#' @param fun function to be used
#' @param arg1,arg2,arg3 Optional arguments to be used by fun
#' @return
#' Numeric vector
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @family laplace functions
#' @export
#' @examples
#' test <- function(p,...){1/p}
#' coeffs8 <- stehfest_coeffs(8)
#' td <- logseq(-2,4,30)
#' stehfest_inversion(td, coeffs8,test)
stehfest_inversion <- function(t, coeffs, fun, arg1 = NULL, arg2 = NULL, arg3= NULL){
  nt <- length(t)
  ft <- vector("numeric", length = nt)
  for(it in 1:nt){
    ft[it] <- 0.0
    a <- log(2.0)/t[it]
    for(ic in 1:length(coeffs)){
      p <- ic*a
      if(missing(arg3)){
        par <- list("p" = p, "arg1" = arg1, "arg2" = arg2)
      }
      else {
        par <- list("p" = p, "arg1" = arg1, "arg2" = arg2, "arg3" = arg3)
      }
      
      ft[it] <- ft[it] + coeffs[ic]*do.call("fun", par)
    }
    ft[it] <- a*ft[it]
  }
  return(ft)
}