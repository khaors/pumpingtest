#include <Rcpp.h>
using namespace Rcpp;
//
// [[Rcpp::export]]
double factorial_cpp(const int& n){
  double fac=1.0;
  if(n == 1){
    return(n);
  }
  else{
    for(int j = 2;j <= n;j++){
      fac = fac * double(j);
    }
  }
  return(fac);
}

// [[Rcpp::export]]
NumericVector stehfest_coefficients_cpp(const int& n){
  NumericVector coeffs(n);
  int nhlf,k1,k2;
  double num, den;
  nhlf = n/2;
  for(int i = 1; i <= n;i++){
    coeffs[i-1] = 0.0;
    k1 = std::floor((i+1)/2);
    k2 = std::min(i, nhlf);
    //std::cout << i << " " << k1 << " " << k2 << std::endl;
    for(int k = k1; k <= k2; k++){
      num = pow(k,nhlf)*factorial_cpp(2*k);
      den = factorial_cpp(nhlf-k)*factorial_cpp(k)*factorial_cpp(k-1)*factorial_cpp(i-k)*factorial_cpp(2*k-i);
      coeffs[i-1] = coeffs[i-1]+num/den;
    }
    coeffs[i-1] = (pow(-1,(i+nhlf)))*coeffs[i-1];
  }
  return(coeffs);
}
// /*** R
// library(pumpingtest)
//   library(microbenchmark)
//   microbenchmark(
//     stehfest_coeffs(8),
//     stehfest_coefficients_cpp(8)
//   )
//   */
// [[Rcpp::export]]
NumericVector theis_WF_LT_cpp(NumericVector p, 
                              const double& par1=0.0, 
                              const double& par2=0.0, 
                              const double& par3=0.0){
  int n=p.size();
  NumericVector L(n);
  //
  for(int i = 0;i<n;i++){
    L[i] = (2.0/p[i])*R::bessel_k(2.0*sqrt(p[i]), 0.0,1);
  }
  return(L);
}
// [[Rcpp::export]]
NumericVector papadopulos_cooper_WF_LT_cpp(NumericVector p, 
                                           const double& alpha, 
                                           const double& rho, 
                                           const double& par=0.0){
  int n=p.size();
  NumericVector L(n);
  double num,den1,den2,currentp,var1,var2,den;
  //
  for(int i =0;i<n;i++){
    currentp = p[i];
    num = R::bessel_k(2.0*sqrt(currentp),0.0,1);
    var1 = (sqrt(currentp)/rho);
    den1 = R::bessel_k(var1,1.0,1);
    var2 = currentp/(alpha*rho*rho);
    den2 = R::bessel_k(var1,0.0,1);
    den = currentp*(var1*den1+var2*den2);
    L[i] = num/den;
  }
  return(L);
}
// [[Rcpp::export]]
NumericVector boulton_WF_LT_cpp(NumericVector p,
                                const double& par1, 
                                const double& par2, 
                                const double& par3=0.0){
  int n=p.size();
  NumericVector L(n);
  double sigma = fabs(par1);
  double phi = fabs(par2);
  NumericVector p1 = abs(p);
  double sp;
  for(int i = 0;i<n;i++){
    sp = sqrt( p1[i] + (phi*p1[i])/(sigma*(p1[i]+phi)) );
    L[i] = R::bessel_k(sp, 0.0, 1)/p1[i];
  }
  return(L);
}

// [[Rcpp::export]]
NumericVector hantush_jacob_WF_LT_cpp(NumericVector p,
                                      const double& par1, 
                                      const double& par2 = 0.0, 
                                      const double& par3 = 0.0){
  int n=p.size();
  NumericVector L(n);
  double pm, sp;
  double beta = par1;
  for(int i =0;i<n;i++){
    pm = fabs(4.0*p[i]+pow(beta,2));
    sp = sqrt(pm);
    L[i] = (2.0/p[i])*R::bessel_k(sp, 0.0, 1);
  }
  return(L);
}

// [[Rcpp::export]]
NumericVector agarwal_recovery_WF_LT_cpp(NumericVector p, 
                                         const double& par1 = 0.0,
                                         const double& par2 = 0.0, 
                                         const double& par3 = 0.0){
  NumericVector W= theis_WF_LT_cpp(p, par1, par2, par3);
  return(W);
}

// [[Rcpp::export]]
NumericVector agarwal_skin_WF_LT_cpp(NumericVector p, 
                                     const double& par1,
                                     const double& par2,
                                     const double& par3){
  double cd, rd, sigma;
  cd = par1;
  rd = par2;
  sigma = par3;
  int n = p.size();
  NumericVector W(n);
  double sp, k0, k1;
  for(int i = 0;i < n;i++){
    sp = sqrt(p[i]);
    k0 = R::bessel_k(sp, 0.0, 1);
    k1 = R::bessel_k(sp, 1.0, 1);
    //std::cout << sp << " " << k0 << " " << k1 << std::endl;
    W[i] = R::bessel_k(rd*sp, 0.0, 1)/(p[i]*(((1.0+p[i]*cd*sigma)*sp*k1)+(cd*p[i]*k0)));
  }
  return(W);
}
// [[Rcpp::export]]
NumericVector general_radial_flow_WF_LT_cpp(NumericVector p, 
                                            const double& par1 = 0.0, 
                                            const double& par2 = 0.0, 
                                            const double& par3 = 0.0){
  double n, rd, sp;
  n = par1;
  rd = par2;
  int np;
  np = p.size();
  NumericVector L(np);
  double term1, term2, cexp, currentp;
  cexp = (float(n)/2.0-1.0);
  //std::cout << "gamma= " << R::gammafn(n/2.0) << std::endl;
  for(int i =0; i < np; i++){
    currentp = p[i];
    if(currentp < 0.0)
      currentp = 1.0e-3;
    sp = sqrt(currentp);
    term1 = R::bessel_k(rd*sp, cexp, 1)/currentp/R::gammafn(n/2.0);
    term2 = std::pow(rd, 2.0-n)*std::pow((std::pow(rd, 2.0)*currentp/4), n/4.0-0.5);
    L[i] = term1*term2;
  }
  return(L);
}
// [[Rcpp::export]]
NumericVector cooper_WF_LT_cpp(NumericVector p, 
                               const double& par1 = 0.0, 
                               const double& par2 = 0.0, 
                               const double& par3 = 0.0){
  double cd;//,rd;
  cd = par1;
  //rd = par2;
  int n = p.size();
  NumericVector L(n);
  double sp, currentp, num, den;
  for(int i = 0; i < n; i++){
    currentp = p[i];
    if(currentp < 0.0)
      currentp = 1e-3;
    sp = sqrt(currentp);
    num = cd*R::bessel_k(sp, 0.0, 1);
    den = (cd*currentp*R::bessel_k(sp, 0.0, 1) + sp*R::bessel_k(sp, 0.0, 1));
    L[i] = num/den;
  }
  return(L);
}
// [[Rcpp::export]]
NumericVector neuzil_WF_LT_cpp(NumericVector p, 
                               const double& par1 = 0.0, 
                               const double& par2 = 0.0, 
                               const double& par3 = 0.0){
  double cd;
  cd = par1;
  int n = p.size();
  NumericVector L(n);
  double sp, currentp, k0, k1, num, den;
  for(int i = 0; i < n; i++){
    currentp = p[i];
    if(currentp < 0.0)
      currentp = 1e-3;
    sp = sqrt(currentp);
    k0 = R::bessel_k(sp, 0.0, 1);
    k1 = R::bessel_k(sp, 1.0, 1);
    num = cd*k0;
    den = (cd*currentp*k0 + sp*k1);
    L[i] = num/den;
  }
  return(L);
}

// [[Rxpp::export]]
NumericVector warren_root_WF_LT_cpp(NumericVector p, 
                                    const double& par1 = 0.0, 
                                    const double& par2 = 0.0, 
                                    const double& par3 = 0.0){
  double sigma, lambda;
  sigma = par1;
  lambda = par2;
  int n = p.size();
  NumericVector L(n);
  double currentp,sp;
  for(int i = 0; i < n; i++){
    currentp = p[i];
    if(currentp < 0.0)
      currentp = 1e-3;
    sp = sqrt(currentp + (lambda*sigma*currentp)/(sigma*currentp+lambda));
    L[i] = (1.0/currentp)*R::bessel_k(sp, 0.0, 1);
  }
  return(L);
  //s <- (1./p)*besselK(sqrt(p+(lambda*sigma*p)/(sigma*p+lambda)), 0)
}

// /*** R
// library(microbenchmark)
//   td <- pumpingtest::logseq(-2,3,100)
//   microbenchmark(
//     pumpingtest::theis_WF_LT(1/td),
//     theis_WF_LT_cpp(1/td)
//   )
//   */
// 
// 
// /*** R
// library(pumpingtest)
//   td <- logseq(-2,3,100)
//   microbenchmark(
//     papadopulos_cooper_F_WF_LT(1/td,0.5,0.5),
//     papadopulos_cooper_WF_LT_cpp(1/td,.5,.5)
//   )
//   */
// 
// /*** R
// td <- pumpingtest::logseq(-2,3,100)
//   microbenchmark(
//     pumpingtest::boulton_WF_LT(1/td,0.5,0.5),
//     boulton_WF_LT_cpp(1/td,.5,.5)
//   )
//   */
// 
// /*** R
// td <- pumpingtest::logseq(-2,3,100)
//   microbenchmark(
//     pumpingtest::hantush_jacob_WF_LT(1/td, 0.5, 0.0, 0.0),
//     hantush_jacob_WF_LT_cpp(1/td,.5)
//   )
//   */
// 
// /*** R
// par <- list(coeffs = stehfest_coeffs(8), rd = .1, cd = .1, sigma = .1)
//   microbenchmark(
//     agarwal_skin_WF_LT_cpp(1/td, .1, .1, .1),
//     agarwal_skin_WF_LT(1/td, .1, .1, .1)
//   )
//   LT <- agarwal_skin_WF_LT_cpp(1/td, .1, .1, .1)
//   LT1 <- agarwal_skin_WF_LT(1/td, .1, .1, .1)
//   print(mean(abs(LT-LT1)))
//   */
typedef NumericVector (*funcPtr)(NumericVector p, 
                       const double& par1,
                       const double& par2,
                       const double& par3);

XPtr<funcPtr> putFunPtrInXPtr(std::string fstr) {
  if (fstr == "theis")
    return(XPtr<funcPtr>(new funcPtr(&theis_WF_LT_cpp)));
  else if (fstr == "boulton")
    return(XPtr<funcPtr>(new funcPtr(&boulton_WF_LT_cpp)));
  else if (fstr == "hantush_jacob")
    return(XPtr<funcPtr>(new funcPtr(&hantush_jacob_WF_LT_cpp)));
  else if(fstr == "papadopulos_cooper")
    return(XPtr<funcPtr>(new funcPtr(&papadopulos_cooper_WF_LT_cpp)));
  else if(fstr == "general_radial_flow" || fstr == "grf")
    return(XPtr<funcPtr>(new funcPtr(&general_radial_flow_WF_LT_cpp)));
  else if(fstr == "agarwal_skin")
    return(XPtr<funcPtr>(new funcPtr(&agarwal_skin_WF_LT_cpp)));
  else if(fstr == "agarwal_recovery")
    return(XPtr<funcPtr>(new funcPtr(&agarwal_recovery_WF_LT_cpp)));
  else if(fstr == "cooper")
    return(XPtr<funcPtr>(new funcPtr(&cooper_WF_LT_cpp)));
  else if(fstr == "neuzil")
    return(XPtr<funcPtr>(new funcPtr(&neuzil_WF_LT_cpp)));
  else if(fstr == "warren_root")
    return(XPtr<funcPtr>(new funcPtr(&warren_root_WF_LT_cpp)));
  else
    return XPtr<funcPtr>(R_NilValue); // runtime error as NULL no XPtr
}
//
// [[Rcpp::export]]
NumericVector callViaString(NumericVector x, 
                            const double& arg1, 
                            const double& arg2,
                            const double& arg3,
                            std::string funname) {
  XPtr<funcPtr> xpfun = putFunPtrInXPtr(funname);
  funcPtr fun = *xpfun;
  NumericVector y = fun(x, arg1, arg2, arg3);
  return (y);
}

// /*** R
// callViaString(1/td, 0.0, 0.0, 0.0, "theis")
//   callViaString(1/td, 0.5, 0.5, 0.0, "boulton")
//   */


// [[Rcpp::export]]
NumericVector callViaXPtr(NumericVector x, 
                          const double& arg1, 
                          const double& arg2,
                          const double& arg3,
                          SEXP xpsexp) {
  XPtr<funcPtr> xpfun(xpsexp);
  funcPtr fun = *xpfun;
  NumericVector y = fun(x, arg1, arg2, arg3);
  return (y);
}

// /*** R
// fun1 <- putFunPtrInXPtr("theis")
//   r0 <- callViaXPtr(1/td, 0.0, 0.0, 0.0, fun1)
//   r1 <- theis_WF_LT(1/td)
//   print(sum(r0-r1))
//   fun2 <- putFunPtrInXPtr("boulton")
//   r0 <- callViaXPtr(1/td, 0.5, 0.5, 0.0, fun2)
//   r1 <- boulton_WF_LT(1/td, .5, .5)
//   print(sum(r0-r1))
//   */

// [[Rcpp::export]]
NumericVector stehfest_inversion_vector_cpp(NumericVector t, 
                                            NumericVector coeffs, 
                                            std::string funname, 
                                            NumericVector arg1, 
                                            NumericVector arg2, 
                                            NumericVector arg3){
  int n=t.length();
  int nc = coeffs.length();
  NumericVector ft(n);
  NumericVector y(1),p1(1);
  double a, arg1v,arg2v,arg3v;
  arg1v=arg1[0];
  arg2v=arg2[0];
  arg3v=arg3[0];
  for(int it = 1;it<=n;it++){
    ft[it-1] = 0.0;
    a = log(2.0)/t[it-1];
    if(arg1.length() > 1){
      arg1v=arg1[it];
    }
    if(arg2.length() > 1){
      arg2v=arg2[it];
    }
    if(arg3.length() > 1){
      arg3v=arg3[it];
    }
    for(int ic = 1; ic <= nc; ic++){
      p1(0) = ic*a;
      y = callViaString(p1, arg1v, arg2v, arg3v, funname); 
      ft[it-1] = ft[it-1] + coeffs[ic-1]*y[0]; 
    }
    ft[it-1] = a*ft[it-1];
  }
  return(ft);
}

// [[Rcpp::export]]
NumericVector stehfest_inversion_cpp(NumericVector t, 
                                     NumericVector coeffs, 
                                     std::string funname,
                                     const double& arg1, 
                                     const double& arg2, 
                                     const double& arg3){
  int n=t.length();
  int nc = coeffs.length();
  NumericVector ft(n);
  NumericVector y(1),p1(1);
  double a;//,p;
  for(int it = 1;it<=n;it++){
    ft[it-1] = 0.0;
    a = log(2.0)/t[it-1];
    for(int ic = 1; ic <= nc; ic++){
      p1(0) = ic*a;
      y = callViaString(p1, arg1, arg2, arg3, funname); 
      ft[it-1] = ft[it-1] + coeffs[ic-1]*y[0]; 
    }
    ft[it-1] = a*ft[it-1];
  }
  return(ft);
}
// /*** R
// td <- pumpingtest::logseq(-2, 4, 50)
//   coeffs <- stehfest_coeffs(8)
//   coeffs1 <- stehfest_coefficients_cpp(8)
//   res.cpp <- stehfest_inversion_cpp(td, coeffs1, "theis", 0.0, 0.0, 0.0)
//   res <- stehfest_inversion(td, coeffs, theis_WF_LT, 0.0, 0.0)
//   
//   microbenchmark(
//     stehfest_inversion(td, coeffs, theis_WF_LT, 0.0, 0.0),
//     stehfest_inversion_cpp(td, coeffs1, "theis", 0.0, 0.0, 0.0)
//   )
//   */

// [[Rcpp::export]]
NumericVector theis_well_function_cpp(NumericVector td, 
                                      const double& par1 = 0.0,
                                      const double& par2 = 0.0,
                                      const double& par3 = 0.0){
  int n=td.size();
  NumericVector W(n);
  NumericVector coeffs = stehfest_coefficients_cpp(8);
  W = stehfest_inversion_cpp(td, coeffs, "theis", par1, par2, par3);
  return(W);
}
// [[Rcpp::export]]
double erf_cpp(double x)
{
  // constants
  double a1 =  0.254829592;
  double a2 = -0.284496736;
  double a3 =  1.421413741;
  double a4 = -1.453152027;
  double a5 =  1.061405429;
  double p  =  0.3275911;
  
  // Save the sign of x
  int sign = 1;
  if (x < 0)
    sign = -1;
  x = fabs(x);
  
  // A&S formula 7.1.26
  double t = 1.0/(1.0 + p*x);
  double y = 1.0 - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t*exp(-x*x);
  
  return sign*y;
}

// /*** R
// W <- theis_well_function_cpp(td, 0.0, 0.0, 0.0)
//   W1 <- theis_well_function(td)
//   print(mean(abs(W-W1)))
//   */

// [[Rcpp::export]]
NumericVector boulton_well_function_cpp(NumericVector td, 
                                        const double& par1, 
                                        const double& par2, 
                                        const double& par3 = 0.0){
  int n=td.size();
  NumericVector W(n);
  NumericVector coeffs = stehfest_coefficients_cpp(8);
  W = stehfest_inversion_cpp(td, coeffs, "boulton", par1, par2, par3);
  return(W);
}
//
// [[Rcpp::export]]
NumericVector boulton_well_function_vector_cpp(NumericVector td, 
                                               NumericVector par1, 
                                               NumericVector par2, 
                                               NumericVector par3){
  int n=td.size();
  NumericVector W(n);
  NumericVector coeffs = stehfest_coefficients_cpp(8);
  W = stehfest_inversion_vector_cpp(td, coeffs, "boulton", par1, par2, par3);
  return(W);
}
// [[Rcpp::export]]
NumericVector hantush_jacob_well_function_cpp(NumericVector td, 
                                              const double& par1, 
                                              const double& par2, 
                                              const double& par3){
  int n=td.size();
  NumericVector W(n);
  NumericVector coeffs = stehfest_coefficients_cpp(8);
  W = stehfest_inversion_cpp(td, coeffs, "hantush_jacob", par1, par2, par3);
  return(W);
}
//
// [[Rcpp::export]]
NumericVector hantush_jacob_well_function_vector_cpp(NumericVector td, 
                                                     NumericVector par1, 
                                                     NumericVector par2, 
                                                     NumericVector par3){
  int n=td.size();
  NumericVector W(n);
  NumericVector coeffs = stehfest_coefficients_cpp(8);
  W = stehfest_inversion_vector_cpp(td, coeffs, "hantush_jacob", par1, 
                                    par2, par3);
  return(W);
}
//[[Rcpp::export]]
NumericVector general_radial_flow_well_function_cpp(NumericVector td,
                                                    const double& par1, 
                                                    const double& par2, 
                                                    const double& par3){
  int n = td.size();
  NumericVector W(n);
  NumericVector coeffs = stehfest_coefficients_cpp(8);
  W = stehfest_inversion_cpp(td, coeffs, "general_radial_flow", par1, par2, par3);
  return(W);
}
// [[Rcpp::export]]
NumericVector papadopulos_cooper_well_function_cpp(NumericVector td, 
                                                   const double& par1, 
                                                   const double& par2, 
                                                   const double& par3){
  int n = td.size();
  NumericVector W(n);
  NumericVector coeffs = stehfest_coefficients_cpp(8);
  W = stehfest_inversion_cpp(td, coeffs, "papadopulos_cooper", par1, par2, 
                             par3);
  return(W);
}
//
// [[Rcpp::export]]
NumericVector papadopulos_cooper_well_function_vector_cpp(NumericVector td, 
                                                   NumericVector par1, 
                                                   NumericVector par2, 
                                                   NumericVector par3){
  int n = td.size();
  NumericVector W(n);
  NumericVector coeffs = stehfest_coefficients_cpp(8);
  W = stehfest_inversion_vector_cpp(td, coeffs, "papadopulos_cooper", par1, 
                                    par2, par3);
  return(W);
}
//
// [[Rcpp::export]]
NumericVector cooper_well_function_cpp(NumericVector td, 
                                       const double& par1, 
                                       const double& par2, 
                                       const double& par3){
  int n = td.size();
  NumericVector W(n);
  NumericVector coeffs = stehfest_coefficients_cpp(8);
  W = stehfest_inversion_cpp(td, coeffs, "cooper", par1, par2, par3);
  return(W);
}
// [[Rcpp::export]]
NumericVector neuzil_well_function_cpp(NumericVector td, 
                                       const double& par1, 
                                       const double& par2, 
                                       const double& par3){
  int n = td.size();
  NumericVector W(n);
  NumericVector coeffs = stehfest_coefficients_cpp(8);
  W = stehfest_inversion_cpp(td, coeffs, "neuzil", par1, par2, par3);
  return(W);
}
// [[Rcpp::export]]
NumericVector agarwal_skin_well_function_cpp(NumericVector td,
                                             const double& par1, 
                                             const double& par2, 
                                             const double& par3){
  int n = td.size();
  NumericVector W(n);
  NumericVector coeffs = stehfest_coefficients_cpp(8);
  W = stehfest_inversion_cpp(td, coeffs, "agarwal_skin", par1, par2, par3);
  return(W);
}

// [[Rcpp::export]]
NumericVector warren_root_well_function_cpp(NumericVector td, 
                                            const double& par1, 
                                            const double& par2, 
                                            const double& par3){
  int n = td.size();
  NumericVector W(n);
  NumericVector coeffs = stehfest_coefficients_cpp(8);
  W = stehfest_inversion_cpp(td, coeffs, "warren_root", par1, par2, par3);
  return(W);
}                    
// /*** R
// par <- list(coeffs = stehfest_coeffs(8), sigma =0.1, phi = 0.1)
//   microbenchmark(
//     boulton_well_function_cpp(td, .1, .1, 0.0),
//     boulton_well_function(td, par)
//   )
//   W <- boulton_well_function_cpp(td, .1, .1, 0.0)
//   W1 <-boulton_well_function(td, par)
//   print(mean(abs(W-W1)))
//   */

// /*** R
// par <- list(coeffs = stehfest_coeffs(8), beta = 0.5)
//   microbenchmark(
//     hantush_jacob_well_function_cpp(td, 0.5, 0.0, 0.0),
//     hantush_jacob_well_function(td, par)
//   )
// */
// 
// 
// /*** R
// W <- hantush_jacob_well_function_cpp(td, .5,0.0,0.0)
// W1 <-hantush_jacob_well_function(td, par)
//   print(mean(abs(W-W1)))
// */
// 
// /*** R
// par <- list(coeffs = stehfest_coeffs(8), n = 1.5, rd = 10)
// microbenchmark(
//   general_radial_flow_well_function_cpp(td, 1.5, 10.0, 0.0),
//   general_radial_flow_well_function(td, par)
// )
// */
// 
// /*** R
// par <- list(coeffs = stehfest_coeffs(8), cd = 0.1)
// microbenchmark(
//   cooper_well_function_cpp(td, 0.1, 0.0, 0.0),
//   cooper_well_function(td, par)
// )
// */
// 
// /*** R
// par <- list(coeffs = stehfest_coeffs(8), cd = 0.1)
// microbenchmark(
//   neuzil_well_function_cpp(td, 0.1, 0.0, 0.0),
//   neuzil_well_function(td, par)
// )
// */
// 
// //[[Rcpp::export]]
// double test_bessel(const double& p){
//   double k0 = R::bessel_k(p,0.0,1.0);
//   double k1 = R::bessel_k(p,1.0,1.0);
//   std::cout << k0 << " " << k1 << std::endl; 
//   return(k0);
// }
// 
// /*** R
// res <- test_bessel(0.057893)
//   print(res)
//   */

//
// [[Rcpp::export]]
NumericMatrix calculate_aquifer_coordinates(const int& nx, 
                                            const int& ny, 
                                            const double& xmin, 
                                            const double& xmax, 
                                            const double& ymin, 
                                            const double& ymax){
  NumericMatrix coords(nx*ny, 2);
  double dx,dy,xcoord,ycoord;
  dx = (xmax-xmin)/(nx);
  dy = (ymax-ymin)/(ny);
  int pos = 0;
  for(int i=0;i<nx;i++){
    for(int j=0;j<ny;j++){
      xcoord = xmin + float(i)*dx;
      ycoord = ymin + float(j)*dy;
      coords(pos,0) = xcoord;
      coords(pos,1) = ycoord;
      pos++;
    }
  }
  return(coords);
}
// [[Rcpp::export]]
NumericVector calculate_distance_well(const double& x0, const double& y0, 
                                      const int& nx, const int& ny, 
                                      const double& xmin, 
                                      const double& xmax, 
                                      const double& ymin, 
                                      const double& ymax){
  NumericVector dist(nx*ny);
  int pos = 0;
  double dx,dy;
  dx = (xmax-xmin)/(nx);
  dy = (ymax-ymin)/(ny);
  double xcoord,ycoord,dinx,diny;
  for(int i=0;i<nx;i++){
    for(int j=0;j<ny;j++){
      xcoord = xmin + float(i)*dx;
      ycoord = ymin + float(j)*dy;
      //Rcout << xcoord << ',' << ycoord << std::endl;
      dinx = std::pow((xcoord-x0), 2.0);
      diny = std::pow((ycoord-y0), 2.0);
      dist[pos]= std::sqrt(dinx + diny);
      pos++;
    }
  }
  return(dist);
}
// [[Rcpp::export]]
NumericVector theis_solution_space(const double& Q, const double& x0, 
                                   const double& y0, const double&t,
                                   NumericVector hydrpar,
                                   const int& nx, const int& ny, 
                                   const double& xmin, 
                                   const double& xmax, 
                                   const double& ymin, 
                                   const double& ymax){
  NumericVector dist(nx*ny), u(nx*ny), td(nx*ny), drawdown(nx*ny);
  dist = calculate_distance_well(x0,y0,nx,ny,xmin,xmax,ymin,ymax);
  double Tr, Ss;
  Tr = hydrpar[0];
  Ss = hydrpar[1];
  u = (Ss*dist*dist)/(4.0*Tr*t);
  td = 1.0/u;
  drawdown=(Q/(4.0*M_PI*Tr))*theis_well_function_cpp(td, 0.0, 0.0, 0.0);
  return(drawdown);
}
//
// [[Rcpp::export]]
NumericVector boulton_solution_space(const double& Q, const double& x0, 
                                     const double& y0, const double&t,
                                     NumericVector hydrpar,
                                     const int& nx, const int& ny, 
                                     const double& xmin, 
                                     const double& xmax, 
                                     const double& ymin, 
                                     const double& ymax){
  NumericVector dist(nx*ny), u(nx*ny), td(nx*ny), drawdown(nx*ny);
  dist = calculate_distance_well(x0,y0,nx,ny,xmin,xmax,ymin,ymax);
  double Tr, Ss, Sy, alpha1;
  NumericVector W(nx*ny), phi(nx*ny), sigma(1), par3(1);
  Tr = hydrpar[0];
  Ss = hydrpar[1];
  Sy = hydrpar[2];
  alpha1 = hydrpar[3];
  sigma(1) = Ss/Sy;
  phi = (alpha1*dist*dist*Ss)/(Tr);
  u = (Ss*dist)/(4.0*Tr*t);
  td = 1.0/u;
  par3(1)=0.0;
  W = boulton_well_function_vector_cpp(td, phi, sigma, par3);
  drawdown = (Q/(4.0*M_PI*Tr))*W;
  return(drawdown);
}
//
// [[Rcpp::export]]
NumericVector papadopulos_solution_space(const double& Q, const double& x0, 
                                         const double& y0, const double&t,
                                         NumericVector hydrpar,
                                         const int& nx, const int& ny, 
                                         const double& xmin, 
                                         const double& xmax, 
                                         const double& ymin, 
                                         const double& ymax){
  NumericVector dist(nx*ny), u(nx*ny), td(nx*ny), drawdown(nx*ny);
  dist = calculate_distance_well(x0,y0,nx,ny,xmin,xmax,ymin,ymax);
  double Tr, Ss, rw, rc;
  NumericVector W(nx*ny),rho(nx*ny),cd(1);
  Tr = hydrpar[0];
  Ss = hydrpar[1];
  rw = hydrpar[2];
  rc = hydrpar[3];
  u = (Ss*dist)/(4.0*Tr*t);
  cd(1) = std::pow(rw, 2)*Ss/std::pow(rc, 2);
  rho = dist/rw;
  td = 1.0/u;
  W = papadopulos_cooper_well_function_vector_cpp(td, cd, rho, 0.0);
  drawdown=(Q/(4.0*M_PI*Tr))*W;
  return(drawdown);
}
//
// [[Rcpp::export]]
NumericVector hantush_jacob_solution_space(const double& Q, const double& x0, 
                                         const double& y0, const double&t,
                                         NumericVector hydrpar,
                                         const int& nx, const int& ny, 
                                         const double& xmin, 
                                         const double& xmax, 
                                         const double& ymin, 
                                         const double& ymax){
  NumericVector dist(nx*ny), u(nx*ny), td(nx*ny), drawdown(nx*ny);
  dist = calculate_distance_well(x0,y0,nx,ny,xmin,xmax,ymin,ymax);
  double Tr, Ss, rw, rc;
  NumericVector W(nx*ny),rho(nx*ny),cd(1);
  Tr = hydrpar[0];
  Ss = hydrpar[1];
  rw = hydrpar[2];
  rc = hydrpar[3];
  u = (Ss*dist)/(4.0*Tr*t);
  cd(1) = std::pow(rw, 2)*Ss/std::pow(rc, 2);
  rho = dist/rw;
  td = 1.0/u;
  W=hantush_jacob_well_function_vector_cpp(td, cd, rho, 0.0);
  drawdown=(Q/(4.0*M_PI*Tr))*W;
  return(drawdown);
}
//
//
typedef NumericVector (*spacefuncPtr)(const double& Q,
                       const double& x0, 
                       const double& y0,
                       const double& t,
                       NumericVector hydrpar, 
                       const int &nx, 
                       const int& ny, 
                       const double& xmin, 
                       const double& xmax, 
                       const double& ymin,
                       const double& ymax);
//
XPtr<spacefuncPtr> putSpaceFunPtrInXPtr(std::string fstr) {
  if (fstr == "theis")
    return(XPtr<spacefuncPtr>(new spacefuncPtr(&theis_solution_space)));
  if (fstr == "boulton")
    return(XPtr<spacefuncPtr>(new spacefuncPtr(&boulton_solution_space)));
  else
    return XPtr<spacefuncPtr>(R_NilValue); // runtime error as NULL no XPtr
}
//
// [[Rcpp::export]]
NumericVector space_calculation_via_string(std::string model, 
                                           const double& Q,
                                           const double& x0, 
                                           const double& y0,
                                           const double& t,
                                           NumericVector hydrpar, 
                                           const int &nx, 
                                           const int& ny, 
                                           const double& xmin, 
                                           const double& xmax, 
                                           const double& ymin,
                                           const double& ymax){
  XPtr<spacefuncPtr> xpfun = putSpaceFunPtrInXPtr(model);
  spacefuncPtr fun = *xpfun;
  NumericVector y = fun(Q, x0, y0, t, hydrpar, nx, ny, xmin, xmax, ymin, ymax);
  return (y);
}
//
// [[Rcpp::export]]
NumericMatrix infinite_aquifer_calculate_drawdown_cpp(std::string model, 
                                                      NumericVector Q, 
                                                      NumericVector x0, 
                                                      NumericVector y0, 
                                                      NumericVector t,
                                                      NumericVector hydrpar,
                                                      const int& nx, 
                                                      const int& ny, 
                                                      const double& xmin, 
                                                      const double& xmax, 
                                                      const double& ymin, 
                                                      const double& ymax){
  
  int ntimes = t.length();
  int nwells = Q.length();
  NumericMatrix drawdown(nx*ny,ntimes);
  double current_x0,current_y0,currentQ,current_time;
  //
  NumericVector current_drawdown(nx*ny),current_drawdown1(nx*ny);
  for(int itime = 0; itime < ntimes; itime++){
    current_time = t[itime];
    current_drawdown1.fill(0.0);
    for(int iwell = 0;iwell<nwells;iwell++){
      currentQ = Q[iwell];
      current_x0 = x0[iwell];
      current_y0 = y0[iwell];
      current_drawdown = space_calculation_via_string(model,currentQ, 
                                                      current_x0, 
                                                      current_y0,
                                                      current_time,
                                                      hydrpar, 
                                                      nx, 
                                                      ny, 
                                                      xmin, 
                                                      xmax, 
                                                      ymin,
                                                      ymax);
      current_drawdown1=current_drawdown1+current_drawdown;
    }
    drawdown(_,itime)=current_drawdown1;
  }
  return(drawdown);
}
// 
// [[Rcpp::export]]
NumericVector constant_head_theis_well_function_cpp(const double& Q, 
                                                    const double& r1, 
                                                    const double& r2,
                                                    NumericVector t,
                                                    NumericVector hydrpar){
  int ntimes = t.length();
  NumericVector u1(ntimes),u2(ntimes), W1(ntimes), W2(ntimes);
  NumericVector t1(ntimes), t2(ntimes), drawdown(ntimes);
  double Tr = hydrpar[0];
  double Ss = hydrpar[1];
  u1 = (r1*r1*Ss)/(4.0*Tr*t);
  u2 = (r2*r2*Ss)/(4.0*Tr*t);
  t1 = 1.0/u1;
  t2 = 1.0/u2;
  W1 = theis_well_function_cpp(t1, 0.0, 0.0, 0.0);
  W2 = theis_well_function_cpp(t2, 0.0, 0.0, 0.0);
  drawdown = (Q/(4.0*M_PI*Tr))*(W1-W2);
  return(drawdown);
}
//
// [[Rcpp::export]]
NumericVector no_flow_theis_well_function_cpp(const double& Q, 
                                              const double& r1, 
                                              const double& r2,
                                              NumericVector t,
                                              NumericVector hydrpar){
  int ntimes = t.length();
  NumericVector u1(ntimes),u2(ntimes), W1(ntimes), W2(ntimes);
  NumericVector t1(ntimes), t2(ntimes), drawdown(ntimes);
  double Tr = hydrpar[0];
  double Ss = hydrpar[1];
  u1 = (r1*r1*Ss)/(4.0*Tr*t);
  u2 = (r2*r2*Ss)/(4.0*Tr*t);
  t1 = 1.0/u1;
  t2 = 1.0/u2;
  W1 = theis_well_function_cpp(t1, 0.0, 0.0, 0.0);
  W2 = theis_well_function_cpp(t2, 0.0, 0.0, 0.0);
  drawdown = (Q/(4.0*M_PI*Tr))*(W1+W2);
  return(drawdown);
}
//
// [[Rcpp::export]]
NumericVector no_flow_boulton_well_function_cpp(const double& Q, 
                                                const double& r1, 
                                                const double& r2,
                                                NumericVector t,
                                                NumericVector hydrpar){
  int ntimes = t.length();
  NumericVector u1(ntimes),u2(ntimes), W1(ntimes), W2(ntimes);
  NumericVector t1(ntimes), t2(ntimes), drawdown(ntimes);
  double Tr = hydrpar[0];
  double Ss = hydrpar[1];
  u1 = (r1*r1*Ss)/(4.0*Tr*t);
  u2 = (r2*r2*Ss)/(4.0*Tr*t);
  t1 = 1.0/u1;
  t2 = 1.0/u2;
  W1 = boulton_well_function_cpp(t1, 0.0, 0.0, 0.0);
  W2 = boulton_well_function_cpp(t2, 0.0, 0.0, 0.0);
  drawdown = (Q/(4.0*M_PI*Tr))*(W1-W2);
  return(drawdown);
}
//

