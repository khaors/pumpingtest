#' pumpingtest: A package for the analysis and evaluation of aquifer tests
#'
#' The pumpingtest package provides functions to analyze and evaluate aquifer tests.
#'
#' The package includes the solution to the most common
#' pumping tests applied in aquifers such as:
#' \itemize{
#' \item Theis and Cooper-Jacob solutions (confined aquifer)
#' \item Hantush-Jacob solution (leaky aquifer)
#' \item Jacob-Lohman solution (artesian confined aquifer)
#' \item Boulton solution (phreatic aquifer)
#' \item Cooper solution (slug tests)
#' \item Agarwal solution (recovery tests)
#' \item Agarwal solution (tests with skin effects)
#' \item General Radial Flow (fractured aquifer)
#' \item Neuzil solution (pulse tests)
#' \item Papadopoulous-Cooper solution (large diameter wells)
#' \item Warren and Root solution (dual porosity aquifer)
#' \item Gringarten solution (single fracture)
#' \item Hvorslev solution (slug tests)
#' \item Bower-Rice solution (slug tests)
#' }
#'
#' The estimation of the hydraulic parameter is achieved by different optimization algorithms
#' including
#' \itemize{
#' \item Nonlinear least squares with the Levenberg-Marquart algorithm (package minpack.lm)
#' \item L-BFGS-B (package optim)
#' \item Simulated Annealing (package GenSA)
#' \item Genetic Algorithms(package GA)
#' \item Particle Swarm Optimization (package PSO)
#' \item Differential Evolution (package DEoptim)
#' \item Estimation of Distribution Algorithm using Copulas (package copulaedas)
#' }
#' with four possible objective functions:
#' \itemize{
#' \item Sum of Squares Residuals
#' \item Mean Absolute Deviation
#' \item Maximum Absolute Deviation
#' \item Maximum Likelihood (under the assumption that the residuals follow a normal distribution).
#' }
#'
#' @section Base functions:
#'  The base function includes the contructor of the S3 class pumping_test and associated functions
#'  to display summaries, print on the console, create diagnostic and estimation plots with
#'  the information from a pumping test, estimate the hydraulic parameters, and predict.
#'
#'  The functions in this section are:
#'
#'  pumping_test, summary, print, plot, fit, evaluate, simulate, confint
#'
#' @docType package
#' @name pumpingtest
NULL
#' @title
#' pumping_test
#' @description
#' Function to create a pumping_test object
#' @param id  a character string defining the Pumping test ID
#' @param Q   rate measured at pumping well (m3/s)
#' @param r  distance to observation well (m)
#' @param t  Numeric vector with the times at which measurements were taken (s)
#' @param s  Numeric vector with the measured drawdown (m)
#' @return A pumping_test object
#' @usage pumping_test(id = character(0), Q = 0.0, r = 0.0, t = NULL, s = NULL)
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @export
#' @family base functions
#' @examples
#' # Define a pumping test
#' tt <- logseq(0, 5, 30)
#' Tr <- 1e-3
#' Ss <- 0.0001
#' r <- 1e3
#' Q <- 1e-3
#' s <- Q/(4*pi*Tr)*theis_well_function(tt)
#' ptest <- pumping_test("Test", Q = Q, r = r, t = tt, s = s)
#' ptest
#' plot(ptest)
pumping_test <- function(id = character(0), Q = 0.0, r = 0.0, t = NULL, s = NULL){
  coeffs <- stehfest_coefficients_cpp(8)
  nt <- length(t)
  ns <- length(s)
  nq <- length(Q)
  if(nt != ns & nt!=nq){
    stop('The t and s/q vectors are not of the same size')
  }
  if(nt == 1){
    stop('t is a vector of a single element')
  }
  if(ns == 1){
    warning('WARNING: s is a vector of a single element. Constant.drawdown test assumed.')
  }
  #
  pos_valid <- t > 1e-12
  if(sum(pos_valid) != nt){
    warning('WARNING: some t values are less than 1e-12. Removing times and drawdown.')
    t <- t[pos_valid]
    if(ns > 1){
      s <- s[pos_valid]
    }
    if(nq > 1){
      Q <- Q[pos_valid]
    }
  }
  # Define aquifer test type
  current.type <- "constant.rate"
  if(class(Q) == "data.frame"){
    current.type <- "variable.rate"
  }
  if(class(Q) != "data.frame" & abs(Q) < 1e-12){
    current.type <- "slug"
  }
  if(class(Q) == "numeric" & length(s) == 0){
    current.type <- "constant.drawdown"
  }
  #
  test <- list(id = id, Q = Q, r = r, t = t, s = s,
               additional_parameters = NULL,
               hydraulic_parameters = NULL,
               parameters = NULL,
               coeffs = coeffs,
               model = character(0),
               estimated = FALSE,
               hydraulic_parameters_names = NULL,
               test.type = current.type)
  class(test) <- "pumping_test"
  invisible(test)
}
#' @title
#' hydraulic.parameters<-
#' @description
#' Function to assign the hydraulic parameters to a pumping_test object
#' @param x A pumping_test object
#' @param value A list or matrix with the hydraulic parameters.
#' @return
#' The pumping_test object with the hydraulic parameters
#' @family base functions
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @export
#' @examples
#' data(theis)
#' ptest.theis <- pumping_test('Well1', Q = 1.38e-2, r = 250, t = theis$t, s = theis$s)
#' ptest.theis.fit <- fit(ptest.theis, 'theis')
#' hydraulic.parameters(ptest.theis) <- ptest.theis.fit$hydraulic_parameters
`hydraulic.parameters<-` <- function(x, value) {
  x$hydraulic_parameters <- value
  return(x)
}
#' @title
#' additional.parameters<-
#' @description
#' Function to assign the additional parameters (well radius, aquifer thickness, ) to a
#' pumping_test object
#' @param x A pumping_test object
#' @param value A list with the additional parameters
#' @return
#' The pumping_test object with the additional parameters
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @family base functions
#' @export
#' @examples
#' data(theis)
`additional.parameters<-` <- function(x, value) {
  x$additional_parameters <- value
  return(x)
}
#' @title
#' fit.parameters<-
#' @description
#' Function to assign the fit parameters (a, t0 and others depending on the model) to a
#' pumping_test object
#' @param x A pumping_test object
#' @param value A list with the parameters obtained by the fit procedure
#' @return
#' The pumping_test object with the fit parameters
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @family base functions
#' @export
#' @examples
#' data(theis)
#' ptest.theis <- pumping_test('Well1', Q = 1.38e-2, r = 250, t = theis$t, s = theis$s)
#' ptest.theis.fit <- fit(ptest.theis, 'theis')
#' hydraulic.parameters(ptest.theis) <- ptest.theis.fit$hydraulic_parameters
#' fit.parameters(ptest.theis) <- ptest.theis.fit$parameters
`fit.parameters<-` <- function(x, value) {
  x$parameters <- value
  return(x)
}
#' @title
#' model<-
#' @description
#' Function to assign the model for estimation to a  pumping_test object
#' @param x A pumping_test object
#' @param value Character sting wit the model used in the intepretation
#' @return
#' The pumping_test object with the model assigned
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @family base functions
#' @export
#' @examples
#' data(theis)
#' ptest.theis <- pumping_test('Well1', Q = 1.38e-2, r = 250, t = theis$t, s = theis$s)
#' ptest.theis.fit <- fit(ptest.theis, 'theis')
#' hydraulic.parameters(ptest.theis) <- ptest.theis.fit$hydraulic_parameters
#' model(ptest.theis) <- 'theis'
`model<-` <- function(x, value) {
  x$model <- value
  return(x)
}
#' @title
#' estimated<-
#' @description
#' Function to define an intepreted pumping_test object
#' @param x A pumping_test object
#' @param value Logical value
#' @return
#' The pumping_test object with estimated variable set to TRUE
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @family base functions
#' @export
#' @examples
#' data(theis)
#' ptest.theis <- pumping_test('Well1', Q = 1.38e-2, r = 250, t = theis$t, s = theis$s)
#' ptest.theis.fit <- fit(ptest.theis, 'theis')
#' hydraulic.parameters(ptest.theis) <- ptest.theis.fit$hydraulic_parameters
#' model(ptest.theis) <- 'theis'
#' estimated(ptest.theis) <- TRUE
`estimated<-` <- function(x, value) {
  x$estimated <- value
  return(x)
}
#' @title
#' hydraulic.parameter.names <-
#' @description
#' Function to assign name of the hydraulic parameters to a  pumping_test object. This is
#' useful in plotting results
#' @param x A pumping_test object
#' @param value A list with the names of the hydraulic parameters
#' @return
#' The pumping_test object with the names of the hydraulic parameters
#' @family base functions
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @export
`hydraulic.parameter.names<-` <- function(x, value) {
  x$hydraulic_parameters_names <- value
  return(x)
}
#' @title
#' summary.pumping_test
#' @description
#' Function to display a short summary of the drawdown data
#' @param object A pumping_test object
#' @param ... additional parameters for the data.frame summary function
#' @return
#' This function displays on the console a summary of the drawdown data
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @family base functions
#' @import stringi
#' @export
#' @examples
#' # Define a pumping test
#' tt <- logseq(0, 5, 30)
#' Tr <- 1e-3
#' Ss <- 0.0001
#' r <- 1e3
#' Q <- 1e-3
#' s <- Q/(4*pi*Tr)*theis_well_function(tt)
#' ptest <- pumping_test("Test", Q = Q, r = r, t = tt, s = s)
#' summary(ptest)
summary.pumping_test <- function(object, ...){
  ptest <- object
  cat('Pumping Test: '%s+%ptest$id%s+%'\n')
  cat('Test Type: '%s+%ptest$test.type%s+%'\n')
  pdata <- as.data.frame(cbind(ptest$t, ptest$s))
  names(pdata) <- c("time","drawdown")
  summary(pdata, ...)
}
#' @title
#' print
#' @description
#' Function to print on the screen a pumping_test object
#' @param x  A pumping_test object
#' @param ... Additional parameters to the data.frame print function
#' @return
#' This function prints the information from the pumping test on the screen
#' @usage \\method{print}{pumping_test}(x, ...)
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @family base functions
#' @export
#' @examples
#' # Define a pumping test
#' tt <- logseq(0, 5, 30)
#' Tr <- 1e-3
#' Ss <- 0.0001
#' r <- 1e3
#' Q <- 1e-3
#' s <- Q/(4*pi*Tr)*theis_well_function(tt)
#' ptest <- pumping_test("Test", Q = Q, r = r, t = tt, s = s)
#' print(ptest)
print.pumping_test <- function(x, ...){
  ptest <- x
  cat('Pumping Test: '%s+%ptest$id%s+%'\n')
  cat('Test Type: '%s+%ptest$test.type%s+%'\n')
  if(ptest$test.type != "slug"){
    if(ptest$test.type == "constant.rate"){
      cat('Q = '%s+%ptest$Q%s+%'\n')
    }
    else if(ptest$test.type == "variable.rate"){
      cat("Pumping history\n")
      print(ptest$Q)
    }
  }
  cat('r= '%s+%ptest$r%s+%'\n')
  pdata <- as.data.frame(cbind(ptest$t, ptest$s))
  names(pdata) <- c("time","drawdown")
  print(pdata, ...)
}
#' @title
#' plot.pumping_test
#' @description
#' Function to plot the pumping test data. This function can create two different
#' types of plots: diagnostic and estimation. The diagnostic plot includes the
#' drawdown vs time plot and the derivative of drawdown with respect to the log of
#' time. This derivative can help in the identification of the flow regime that
#' occurred when the data was acquired.
#' @param x  A pumping_test object
#' @param type  Type of plot. Current options include
#' \itemize{
#' \item diagnostic
#' \item estimation
#' \item model.diagnostic
#' \item uncertainty
#' \item mcmc.trace
#' \item mcmc.run_mean
#' \item mcmc.compare
#' \item mcmc.autocorr
#' \item sample.influence
#' }
#' @param d  Derivative parameter. If method is bourdet then d is a parameter to specify
#' the number of lags in the derivative. If method is spline then d is the number of points
#' used to calculate the derivative.
#' @param dmethod  Method to calculate the derivative (central, horner, bourdet, spline)
#' @param scale  Option to define a loglog or semilog diagnostic plot
#' @param y.intersp  Numeric value to define the interspacing between lines in the legend
#' @param slug  Logical flag to indicate a slug test
#' @param legend Logical flag to indicate if legend is included (only for estimation plot)
#' @param results Logical flag to indicate if the estimation results are going to be included in the estimation plot
#' @param cex character expansion factor relative to current par("cex"). This is a parameter of the plot functions.
#' @param ...  Additional parameters for the plot, points and lines functions.
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @importFrom ggplot2 ggplot geom_point geom_line labs xlab ylab ggtitle aes theme_bw
#' @importFrom ggplot2 scale_x_log10 scale_y_log10 scale_shape_manual scale_color_manual
#' @importFrom dplyr full_join
#' @importFrom stringr str_to_upper str_pad str_length str_c
#' @family base functions
#' @export
#' @examples
#' # Define a pumping test
#' data(theis)
#' ptest <- pumping_test("Test", Q = 1.388e-2, r = 250, t = theis$t, s = theis$s)
#' # Diagnostic plot using default parameters
#' plot(ptest)
#' # Diagnostic plot with Horner derivative
#' plot(ptest, dmethod = 'horner')
#' # Diagnostic plot with Bourdet derivative d = 3
#' plot(ptest, dmethod = 'bourdet', d = 3)
#' # Diagnostic plot with Spline derivative
#' plot(ptest, dmethod = 'spline', d = 20)
#' # Diagnostic plot with semilog scale
#' plot(ptest, scale = 'slog')
#' #estimation Plot
#' ptest.fit <- fit(ptest, "theis")
#' hydraulic.parameters(ptest) <- ptest.fit$hydraulic_parameters
#' fit.parameters(ptest) <- ptest.fit$parameters
#' model(ptest) <- "theis"
#' estimated(ptest) <- TRUE
#' plot(ptest, type = 'estimation', dmethod = "spline", d = 30, results = FALSE)
#' # Model Diagnostic plot
#' plot(ptest, type = 'model.diagnostic')
#' # Uncertainty plot (bootstrap)
#' ptest.confint <- confint(ptest, level = c(0.025, 0.975), method = 'bootstrap', n = 30, neval = 100)
#' hydraulic.parameters(ptest) <- ptest.confint$hydraulic.parameters
#' hydraulic.parameter.names(ptest) <- ptest.confint$hydraulic.parameters.names
#' plot(ptest, type = 'uncertainty')
plot.pumping_test <- function(x, type = c('diagnostic','estimation',
                                          'model.diagnostic', 'uncert',
                                          'mcmc.trace', 'mcmc.run_mean',
                                          'mcmc.compare', 'mcmc.autocorr',
                                          'sample.influence'),
                              dmethod = 'central', d = 2, scale = 'loglog',
                              y.intersp = 0.5, slug = FALSE, legend = TRUE,
                              results = FALSE, cex = 1, ...){
  if(class(x) != "pumping_test"){
    stop('ERROR: a pumpingtest object is required as input')
  }
  ptest <- x
  size <- cex*2.5
  variable <- NULL
  shape <- NULL
  colors <- NULL
  title <- NULL
  # create dataframe with all information
  t <- ptest$t
  s <- ptest$s
  ndat <- length(t)
  variable1 <- vector('character', length = ndat)
  variable1[1:ndat] <- 'Drawdown'
  df.drawdown <- data.frame(t= t, s = s, variable = variable1)
  # create dataframe with derivative
  if(slug){
    ds <- log_derivative(t, -s, d = d, method = dmethod)
  }
  else {
    ds <- log_derivative(t, s, d = d, method = dmethod)
  }
  dt <- ds$x
  ds <- ds$y
  nder <- length(dt)
  variable2 <- vector('character', length = nder)
  variable2[1:nder] <- str_c('Derivative(',dmethod,')')
  df.derivative <- data.frame(t = dt, s = ds, variable = variable2)
  # Join data.frames
  suppressWarnings(
    ptest.def <- full_join(df.drawdown, df.derivative)
  )
  #
  #print(ptest.def)
  # Define title
  type <- type[1]
  if(type == 'diagnostic'){
    if(is.null(title)){
      ctitle <- str_c("Diagnostic Plot: ", ptest$id)
    }
    else {
      ctitle <- title
    }
  }
  else if(type == 'estimation'){
    if(is.null(title)){
      ctitle <- str_c("Estimation Plot: ", ptest$id)
    }
    else {
      ctitle <- title
    }
  }
  #
  if(is.null(colors)){
    colors <- c('#00BFC4', '#F8766D')
  }
  #
  if(is.null(shape)){
    shape <- c(17, 16)
  }
  #
  # Create diagnostic plot
  p1 <- NULL
  if(type == 'diagnostic'){
    p1 <- ggplot(data = ptest.def) + geom_point(aes(x = t, y = s, group = variable,
                                                    color = variable,
                                                    shape = variable), size = size) +
      scale_x_log10()+
      labs(title = ctitle) +
      xlab("Time(s)") +
      ylab("Drawdown(m)") +
      scale_shape_manual(values = rev(shape)) +
      scale_color_manual(values = rev(colors))
    if(scale == 'loglog'){
      p1 <- p1 + scale_y_log10()
    }
    p1 <- p1 + theme_bw()
  }
  # Create Estimation plot
  if(type == 'estimation'){
    # Calculate drawdown from hydraulic parameters
    if(ptest$estimated){
      #tpred <- logseq(from = 0.99*log10(min(ptest$t)), to = 1.01*log10(max(ptest$t)), 100)
      ptest.pred <- evaluate(ptest, ptest$model, FALSE, n = 200)
      ptest.pred.df <- data.frame(t = ptest.pred$t[1:199],
                                  s = ptest.pred$s[1:199],
                                  variable = "estimation")
      dptest.pred.df <- data.frame(t = ptest.pred$t[1:199],
                                   s = ptest.pred$dsdlnt[1:199])
    }
    else {
      stop('ERROR: A estimated hydraulic parameters are required for estimation plot')
    }
    #
    p1 <- ggplot(data = ptest.def) + geom_point(aes(x = t, y = s, group = variable,
                                                    color = variable,
                                                    shape = variable), size = size) +
      scale_x_log10() +
      labs(title = ctitle) +
      xlab("Time(s)") +
      ylab("Drawdown(m)") +
      geom_line(aes(x = t, y = s), ptest.pred.df, color = colors[1]) +
      geom_line(aes(x = t, y = s), dptest.pred.df, color = colors[2]) +
      scale_shape_manual(values = rev(shape)) +
      scale_color_manual(values = rev(colors)) +
      theme_bw()
    if(scale == 'loglog'){
      p1 <- p1 + scale_y_log10()
    }
    p1 <- p1 + theme_bw()
  }#Estimation
  else if(type == 'model.diagnostic'){
    p1 <- plot_model_diagnostic(ptest, cex = cex, ...)
  }
  else if(type == 'uncertainty'){
    p1 <- plot_uncert(ptest)
  }
  #
  return(p1)
}
#' @title
#' plot_model_diagnostic
#' @description
#' Function to plot the residuals of an estimated pumping test
#' @param ptest A pumping_test object. It must be estimated.
#' @param ... Additional parameters to the plot function used in scatter.smooth
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @family base functions
#' @importFrom stats  sd
#' @importFrom ggplot2 ggplot geom_point coord_equal ggtitle geom_smooth geom_qq
#' @importFrom gridExtra arrangeGrob grid.arrange
#' @export
#' @examples
#' data(theis)
#' ptest <- pumping_test("Test", Q = 1.388e-2, r = 250, t = theis$t, s = theis$s)
#' ptest.fit <- fit(ptest, "theis")
#' hydraulic.parameters(ptest) <- ptest.fit$hydraulic_parameters
#' fit.parameters(ptest) <- ptest.fit$parameters
#' model(ptest) <- "theis"
#' estimated(ptest) <- TRUE
#' plot(ptest, type = 'estimation', dmethod = "spline", d = 30)
#' plot_model_diagnostic(ptest)
plot_model_diagnostic <- function(ptest, ...){
  if(class(ptest) != 'pumping_test'){
    stop('A pumping_test object is required as input')
  }
  if(!ptest$estimated){
    stop('An estimated pumping_test object is required as input')
  }
  res <- evaluate(ptest, ptest$model)
  s.pred <- res$sd
  s.residuals <- s.pred - ptest$s
  s.std.residuals <- s.residuals/sd(s.residuals)
  s.std.residuals1 <- abs(s.residuals) #sqrt(abs(s.std.residuals))
  # Define global variables
  measured <- NULL
  calculated <- NULL
  residuals <- NULL
  abs.residuals <- NULL
  # Plot measured rho vs calculated rho
  df1 <- data.frame(measured = ptest$s, calculated = s.pred,
                    residuals = s.residuals,
                    abs.residuals = s.std.residuals1)
  p1 <- ggplot() + geom_point(aes(x = measured, y = calculated), data = df1,
                              color = "red") +
    coord_equal() +
    ggtitle("a) Measured vs Calculated") +
    theme_bw()
  #
  p2 <- ggplot(data = df1, aes(x = calculated, y = residuals)) +
    geom_point(color = "red") +
    geom_smooth() +
    ggtitle("b) Residuals") +
    theme_bw()
  #
  p3 <- ggplot(data = df1, aes(x = calculated, y = abs.residuals)) +
    geom_point(color = "red") +
    geom_smooth() +
    ggtitle("c) Absolute Residuals") +
    theme_bw()
  #
  p4 <- ggplot(data = df1, aes(sample = residuals)) + geom_qq(color = "red") +
    ggtitle("d) QQ plot") +
    theme_bw()
  #
  #ptot <- arrangeGrob(p1, p2, p3, p4, ncol = 2)
  ptot <- grid.arrange(p1, p2, p3, p4, ncol = 2)
  return(ptot)
}
#' @title
#' plot_sample_influence
#' @description
#' Function to create a plot with an influence measure for each sample of the pumping test.
#' This plot is helpful in the identification of the samples that are influential in the
#' estimation of the hydraulic parameters.
#' @param res A list with the results of the jackniffe CI estimation.
#' @param ... Additional parameters to the plot function
#' @return
#' A plot with the influence measure of each sample
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @family base functions
#' @importFrom grDevices colors
#' @importFrom graphics plot abline lines points legend title text par mtext
#' @export
plot_sample_influence <- function(res, ...){
  if(res$method != 'jackniffe'){
    stop('ERROR the results of confint_jackniffe are required to create this plot.')
  }
  par.label.names <- list(Tr = 'Transmissivity(m2/s)',
                          Ss = 'Storage Coefficient',
                          radius_influence = 'Radius Influence(m)',
                          r = 'Distance to well(m)',
                          Q = 'Discharge(m3/s)',
                          omegad = 'Drainage Porosity',
                          cd = 'Wellbore storage',
                          rho = 'Dimensionless radius',
                          Ka = 'K aquitard (m/s)',
                          B = 'Aquitard Thickness(m)',
                          rw = 'Well radius(m)',
                          rc = 'Casing radius(m)',
                          rd = 'Dimensionless radius',
                          Sxf2 = 'Sxf2',
                          n = 'n (Flow Dimension)')
  label <- c('a', 'b', 'c', 'd', 'e', 'f')
  dfbeta <- res$dfbeta
  npar <- ncol(dfbeta)
  ndat <- nrow(dfbeta)
  nrowp <- ceiling(npar/2)
  hydr.par.names <- colnames(dfbeta)
  #
  par0 <- par(no.readonly = TRUE)
  par(mfrow=c(nrowp, 2))
  pos <- seq(1, nrow(dfbeta), by = 1)
  for(ipar in 1:npar){
    mx <- max(max(dfbeta[,ipar]), 2.5)
    plot(dfbeta[,ipar], type = "p", xlab = "Sample", ylab = "DFBETA",
         main = paste0(label[ipar],') ', par.label.names[[hydr.par.names[ipar]]]),
         ylim = c(0, mx), ...)
    lines(c(1,nrow(dfbeta)),c(2.,2.), lty = 3, col = "red", ...)
    pos_influential <- dfbeta[,ipar] > 2
    if(sum(pos_influential) > 0){
      points(pos[pos_influential], dfbeta[pos_influential,ipar], col = "red",
             pch = 7, ...)
    }
    for(i in 1:ndat){
      lines(c(i,i),c(0,dfbeta[i,ipar]), col = "black")
    }
  }
  par(par0)
}
#' @title
#' plot_uncert
#' @description
#' Function to plot the distributions of the confidence intervals estimated using bootstrapp
#' from a model fitted using nonlinear regression.
#' @param ptest A pumping_test object.
#' @param cex Character expansion
#' @param ... Additional parameters for the plot function
#' @importFrom GGally ggpairs
#' @importFrom ggplot2 theme_bw
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @family base functions
#' @export
plot_uncert <- function(ptest, cex =1, ...){
  if(class(ptest) != 'pumping_test'){
    stop('A pumping_test object is required as input')
  }
  #
  if(!ptest$estimated){
    stop('An estimated pumping_test object is required')
  }
  #
  if(class(ptest$hydraulic_parameters) != 'matrix'){
    stop('The values of the hydraulic parameters must be stored in a matrix')
  }
  #
  if(nrow(ptest$hydraulic_parameters) < 30){
    stop('There are not enough values of the hydraulic parameters')
  }
  #
  if(is.null(ptest$hydraulic_parameters_names)){
    stop('ERROR the names of the hydraulic parameters have not been assigned')
  }
  #
  par.label.names <- list(Tr = 'Transmissivity(m2/s)',
                          Ss = 'Storage Coefficient',
                          radius_influence = 'Radius Influence(m)',
                          r = 'Distance to well(m)',
                          Q = 'Discharge(m3/s)',
                          omegad = 'Drainage Porosity',
                          cd = 'Wellbore storage',
                          rho = 'Dimensionless radius',
                          Ka = 'K aquitard (m/s)',
                          B = 'Aquitard Thickness(m)',
                          rw = 'Well radius(m)',
                          rc = 'Casing radius(m)',
                          rd = 'Dimensionless radius',
                          Sxf2 = 'Sxf2',
                          n = 'n (Flow Dimension)')
  #
  hydr.parameters.names <- ptest$hydraulic_parameters_names
  hydr.parameters.names1 <- vector('character', length(hydr.parameters.names))
  for(ipar in 1:length(hydr.parameters.names)){
    current_par <- hydr.parameters.names[[ipar]]
    hydr.parameters.names1[ipar] <- par.label.names[[current_par]]
  }
  #
  hydraulic_parameters.df <- as.data.frame(ptest$hydraulic_parameters)
  names(hydraulic_parameters.df) <- hydr.parameters.names1
  p1 <- ggpairs(hydraulic_parameters.df,
                columnLabels = hydr.parameters.names1) +
    theme_bw()
  return(p1)
}
#' @title
#' fit
#' @description
#' Generic function to estimate the aquifer parameters from a pumping test. This function uses
#' nonlinear least squares to estimate these parameters.
#' @param ptest  A pumping_test object
#' @param model  Character string specifying the model used in the estimation
#' @param control.par A list with parameters of the parameter estimation using nonlinear
#' regression
#' @param trace A logical flag indicating if the results of the nonlinear regression are
#' printed on the screen
#' @return
#' A list with the following entries:
#' \itemize{
#' \item hydraulic_parameters: hydraulic parameters of the model (includes transmissivity, storage
#' coefficient and radius of influence)
#' \item paraemters: fitted parameters (includes a and t0)
#' \item resfit: List with the results of the nonlinear regression
#' \item value: Value of the residual sum of squares
#' }
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @import minpack.lm
#' @importFrom  stats as.formula
#' @family base functions
#' @export
#' @examples
#' #Fit test from confined aquifer
#' data(theis)
#' ptest <- pumping_test('Well1', Q = 1.388e-2, r = 250, t = theis$t, s = theis$s)
#' res_th <- fit(ptest, 'theis')
#' print(res_th)
#' #Fit test from confined aquifer
#' res_cj <- fit(ptest, 'cooper_jacob')
fit <- function(ptest, model, control.par, trace = F){
  drawdown_curve <- as.data.frame(cbind(ptest$t,ptest$s))
  names(drawdown_curve) <- c('t','s')
  #print(names(drawdown_curve))
  initial_solution_fn <- model%s+%'_solution_initial'
  initial_solution <- do.call(initial_solution_fn,list(ptest))
  parameters <- names(initial_solution)
  #npar <- length(parameters)
  #lower <- rep(1e-3, npar)
  current_formula <- 's ~ '%s+%model%s+%'_solution(ptest,'
  for(ipar in parameters){
    current_formula <- current_formula%s+%ipar%s+%','
  }
  current_formula <- current_formula%s+%'t'%s+%')'
  current_formula <- as.formula(current_formula)
  #print(current_formula)
  Q <- ptest$Q
  r <- ptest$r
  t <- ptest$t
  nprint <- 10
  if(missing(control.par)){
    if(!trace){
      nprint <- 0
    }
    control.par <- nls.lm.control(ftol = 1e-10,
                                  ptol = 1e-10,
                                  maxiter = 200,
                                  nprint = nprint,
                                  maxfev = 1000)
  }
  resfit <- nlsLM(current_formula, data = drawdown_curve,
                  start = initial_solution, trace = trace,
                  control =  control.par) #, lower = lower, upper = rep(Inf,npar)
  pars <- resfit$m$getAllPars()
  pars <- as.list(pars)
  fn <- model%s+%'_calculate_parameters'
  args <- list('ptest' = ptest, 'par' = pars)
  res <- do.call(fn, args = args)
  res1 <- list(hydraulic_parameters = res, parameters = pars, resfit = resfit,
               value = resfit$m$deviance())
  return(res1)
}
#' @title
#' fit.optimization
#' @description
#' Function to estimate the aquifer parameters from a pumping test using several optimization
#' functions.
#' @param ptest A pumping_test object.
#' @param model A character string specifying the model used in the parameter estimation.
#' @param obj.fn A character string specifying the objective function used in the parameter estimation.
#' Currently the following objective functions are included:
#' \itemize{
#' \item 'rss': Residual sum of squares
#' \item 'mnad': Mean absolute deviation
#' \item 'mxad': Maximum absolute deviation
#' \item 'loglik': Loglikelihood function
#' }
#' @param opt.method A character string specifying the optimization method used in the parameter estimation.
#' Currently the following methologies are included:
#' \itemize{
#' \item 'nls': Nonlinear regression
#' \item 'sa': Simulated Annealing (GenSA package)
#' \item 'ga': Genetic Algorithms (GA package)
#' \item 'l-bfgs-b': using optim function (stats package)
#' \item 'pso': Particle Swarm Optimization (pso package)
#' \item 'copulaedas': Estimation of Distribution Algorithms Based on Copulas (copulaedas package)
#' \item 'de': Differential Evolution (DEoptim package)
#' }
#' @param lower A numeric vector with the lower values of the search region
#' @param upper A numeric vector with the upper values of the search region
#' @param control.par A list with the parameters of the optimization method
#' @param seed A random seed
#' @return
#' A list with the following entries:
#' \itemize{
#' \item hydraulic_parameters: hydraulic parameters of the model (includes transmissivity, storage
#' coefficient and radius of influence, or others)
#' \item parameters: fitted parameters (including a and t0 and other depending on the model)
#' \item resfit: The list or object returned by the optimization driver of each method.
#' \item value: The value of the objective function reached at the end of the optimization run.
#' }
#' @importFrom stats optim
#' @importFrom GenSA GenSA
#' @importFrom GA ga
#' @importFrom pso psoptim
#' @importFrom copulaedas edaRun edaTerminateMaxGen edaReportSimple VEDA
#' @importFrom DEoptim DEoptim DEoptim.control
#' @importFrom methods setMethod
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @family base functions
#' @export
#' @examples
#' \dontrun{
#' # Define pumping_test object
#' data("boulton")
#' ptest.boulton <- pumping_test("Well1", Q = 0.03, r = 20,
#'                               t = boulton$t, s = boulton$s)
#' # Parameter estimation using L-BFGS-B
#' ptest.boulton.bfgs.rss <- fit.optimization(ptest.boulton,
#' "boulton", obj.fn = "rss", opt.method = "l-bfgs-b",
#' seed = 54321)
#' # Parameter estimation using Simulated Annealing
#' ptest.boulton.sa.rss <- fit.optimization(ptest.boulton,
#' "boulton", obj.fn = "rss", opt.method = "sa", seed = 54321)
#' # Parameter estimation using Genetic Algorithms
#' ptest.boulton.ga.rss <- fit.optimization(ptest.boulton,
#' "boulton", obj.fn = "rss", opt.method = "ga", seed = 54321)
#' # Parameter estimation using Differential Evolution
#' ptest.boulton.de.rss <- fit.optimization(ptest.boulton,
#'                         "boulton", obj.fn = "rss", opt.method = "de", seed = 54321)
#' # Parameter estimation using Particle Swarm Optimization
#' ptest.boulton.pso.rss <- fit.optimization(ptest.boulton,
#'                          "boulton", obj.fn = "rss", opt.method = "pso", seed = 54321)
#' }
fit.optimization <- function(ptest, model, obj.fn = 'rss', opt.method = 'nls',
                             lower = 1e-9, upper = Inf, control.par, seed = 12345){
  if(class(ptest) != 'pumping_test'){
    stop('A pumping_test object is required as input')
  }
  drawdown_curve <- as.data.frame(cbind(ptest$t,ptest$s))
  names(drawdown_curve) <- c('t','s')
  # Determine initial solution using the corresponding function for each model
  initial_solution_fn <- paste(model, '_solution_initial', sep = '')
  initial_solution <- do.call(initial_solution_fn,list(ptest))
  parameters <- names(initial_solution)
  par <- as.numeric(initial_solution)
  # Initialize results list
  res1 <- list()
  # Set Random seed
  set.seed(seed)
  # Add the sigma parameter for loglikelihood estimation
  if(obj.fn == 'loglik' | obj.fn == 'loglik-bc'){
    par[(length(par)+1)] <- 0.1
  }
  # Define lower and upper limits of the search region
  if(missing(lower)){
    lower <- 1e-3*par
  }
  if(missing(upper)){
    upper <- 20.0*par
  }
  
  #print(c(obj.fn, opt.method))
  #res1 <- list()
  # Check optimization methods for loglikelihood estimation
  if(opt.method == 'nls' & obj.fn == 'loglik' |
     opt.method == 'nls' & obj.fn == 'loglik-bc'){
    stop('LogLikelihood function cannot be minimized with nls')
  }
  # Define the corresponding objective functions
  par.obj.fn <- NULL
  if(obj.fn == 'rss'){
    if(opt.method == 'ga'){
      par.obj.fn <- function(par, ptest, model){-1*residual_sum_squares(par, ptest, model) }
    }
    else if(opt.method == 'copulaedas'){
      par.obj.fn1 <- function(par){
        residual_sum_squares(par, ptest, model)
      }
    }
    else{
      par.obj.fn = residual_sum_squares
    }
  }
  else if(obj.fn == 'loglik'){
    if(opt.method == 'ga' ){
      par.obj.fn <- loglikelihood
    }
    else if(opt.method == 'copulaedas'){
      par.obj.fn1 <- function(par){
        -1*loglikelihood(par, ptest, model)
      }
    }
    else{
      par.obj.fn <- function(par, ptest, model){ -1* loglikelihood(par, ptest, model)}
    }
  }
  else if(obj.fn == 'loglik-bc'){
    if(opt.method == 'ga'){
      par.obj.fn <- loglikelihood_bc
    }
    else if(opt.method == 'copulaedas'){
      par.obj.fn1 <- function(par){
        -1*loglikelihood_bc(par, ptest, model)
      }
    }
    else{
      par.obj.fn <- function(par, ptest, model){ -1* loglikelihood_bc(par, ptest, model)}
    }
  }
  else if(obj.fn == 'mnad'){
    if(opt.method == 'ga'){
      par.obj.fn <- function(par, ptest, model){-1*mean_absolute_deviation(par, ptest, model)}
    }
    else if(opt.method == 'copulaedas'){
      par.obj.fn1 <- function(par){
        mean_absolute_deviation(par, ptest, model)
      }
    }
    else{
      par.obj.fn <- mean_absolute_deviation
    }
  }
  else if(obj.fn == 'mxad'){
    if(opt.method == 'ga'){
      par.obj.fn <- function(par, ptest, model){-1*max_absolute_deviation(par, ptest, model)}
    }
    else if(opt.method == 'copulaedas'){
      par.obj.fn1 <- function(par){
        max_absolute_deviation(par, ptest, model)
      }
    }
    else{
      par.obj.fn <- max_absolute_deviation
    }
  }
  #
  if(opt.method == 'nls' & obj.fn == 'rss'){
    current_formula <- paste('s ~ ', model, '_solution(ptest,', sep = '')
    for(ipar in parameters){
      current_formula <- paste(current_formula, ipar, ',', sep = '')
    }
    current_formula <- paste(current_formula, 't', ')', sep = '')
    current_formula <- as.formula(current_formula)
    #print(current_formula)
    Q <- ptest$Q
    r <- ptest$r
    t <- ptest$t
    if(missing(control.par)){
      control.par <- nls.lm.control(ftol = 1e-10,
                                    ptol = 1e-10,
                                    maxiter = 200,
                                    nprint = 0,
                                    maxfev = 1000)
    }
    resfit <- nlsLM(current_formula, data = drawdown_curve,
                    start = initial_solution, trace = F,
                    control = control.par )
    pars <- resfit$m$getAllPars()
    pars <- as.list(pars)
    fn <- paste(model, '_calculate_parameters', sep = '')
    args <- list('ptest' = ptest, 'par' = pars)
    res <- do.call(fn, args = args)
    current_par <- pars
    res1$resfit <- resfit
    res1$value <- resfit$m$deviance()
  }
  else if(opt.method == 'l-bfgs-b'){
    if(missing(control.par)){
      control.par <- list(trace = -1, maxit = 500)
    }
    res.opt <- optim(par = par, fn = par.obj.fn, ptest = ptest, model = model,
                     method = 'L-BFGS-B', lower = lower, upper = upper, control = control.par)
    current_par <- res.opt$par
    res1$resfit <- res.opt
    res1$value <- res.opt$value
  }
  else if(opt.method == 'sa'){
    if(missing(control.par)){
      simple.function.flag <- TRUE
      if(obj.fn == "mxad"){
        simple.function.flag <- FALSE
      }
      control.par <- list( verbose = FALSE, maxit = 100,
                           simple.function = simple.function.flag,
                           max.time = 120)
    }
    res.sa <- GenSA(lower = lower, upper = upper, fn = par.obj.fn, control = control.par,
                    ptest = ptest, model = model)
    current_par <- res.sa$par
    res1$resfit <- res.sa
    res1$value <- res.sa$value
  }
  else if(opt.method == 'ga'){
    set.seed(seed)
    optimArgs.local <- list(method = "L-BFGS-B", potim = 0.4, pressel = 0.6)
    res.ga <- ga(type = "real-valued", fitness = par.obj.fn, ptest = ptest,  model = model,
                 min = lower, max = upper, popSize = 100, maxiter = 100, run = 50,
                 pcrossover = 0.8, pmutation = 0.1, parallel = FALSE,  optim = TRUE,
                 optimArgs = optimArgs.local)
    current_par <- res.ga@solution
    res1$resfit <- res.ga
    res1$value <- -res.ga@fitnessValue
  }
  else if(opt.method == 'pso'){
    if(missing(control.par)){
      if(obj.fn != 'loglik'){
        control.par <- list(abstol=1e-8, trace = 1, maxit.stagnate = 30)
      }
      else {
        control.par <- list(trace = 0, maxit.stagnate = 30, maxit = 300)
      }
    }
    #
    res.pso<- psoptim(par = par, fn = par.obj.fn, ptest = ptest, model = model, lower = lower,
                      upper= upper, control = control.par)
    current_par <- res.pso$par
    res1$resfit <- res.pso
    res1$value <- res.pso$value
  }
  else if(opt.method == 'copulaedas'){
    setMethod("edaTerminate", "EDA", edaTerminateMaxGen)
    setMethod("edaReport", "EDA", edaReportSimple)
    UMDA <- VEDA(vine="CVine",indepTestSigLevel=0.01,
                 copulas = c("normal"),margin = "norm")
    UMDA@name=paste('pumpingtest_', ptest$id, model, sep = '')
    res.eda<- edaRun(UMDA, f = par.obj.fn1, lower = lower, upper = upper)
    current_par <- res.eda@bestSol
    res1$resfit <- res.eda
    res1$value <- res.eda$value
  }
  else if(opt.method == 'de'){
    if(missing(control.par)){
      control.par <- DEoptim.control(strategy=1, NP=100, itermax=300, CR=0.9,
                                     F=0.8, trace=0, storepopfrom=1, storepopfreq=1)
    }
    res.de <- suppressWarnings(DEoptim(par.obj.fn, lower = lower,
                                       upper = upper, control=control.par, ptest, model))
    current_par <- res.de$optim$bestmem
    res1$resfit <- res.de
    res1$value <- res.de$optim$bestval
  }
  pars <- NULL
  for(ipar in 1:length(parameters)){
    par_fn <- paste('pars$', parameters[ipar], '=', as.character(current_par[ipar]), sep='')
    eval(parse(text = par_fn))
  }
  fn <- paste(model, '_calculate_parameters', sep = '')
  args <- list('ptest' = ptest, 'par' = pars)
  hydr.par <- do.call(fn, args = args)
  res1$parameters <- pars
  res1$hydraulic_parameters <- hydr.par
  if(obj.fn == "loglik" | obj.fn == "loglik-bc"){
    res1$parameters1 <- current_par
  }
  return(res1)
}
#' @title
#' fit.sampling
#' @description
#' Function to estimate the aquifer parameters from a pumping test using Markov Chain Monte Carlo
#' simulation.
#' @param ptest A pumping_test object.
#' @param model A character string specifying the model used in the parameter estimation.
#' @param method A character string specifying the sampling method used in the parameter
#' estimation. Currently the following methods are supported:
#' \itemize{
#' \item amcmc: Adaptative Markov Chain Monte Carlo proposed by Brooks and Rosenthal(2009).
#' This method requires the specification of the variances of the proposal distributions (
#' see proposal.sigma below).
#' \item amcmc1: Adaptative Markov Chain Monte Carlo proposed by Bai(2009). This method does
#' not require the specification of any additional parameter.
#' \item twalk: transverse-walk method proposed by Christen and Fox (2010). This methods does
#' not require any additional parameter.
#' }
#' @param prior.pdf A character vector with the names of the probability density functions used
#' for each hydraulic parameter. Currently only 'unif' and 'norm' are implemented.
#' @param prior.parameters A matrix with the parameters of the probability density functions of the
#' hydraulic parameters used in the drawdown calculations. For the uniform distribution, this parameters
#' are the minimun and maximum, whereas for the normal distribution the mean and standard deviation
#' must be specified.
#' @param iterations An integer specifying the number of iterations required for sampling.
#' @param burnIn An integer specifying the length of the burn-in period of the Markov Chain.
#' @param seed A random seed.
#' @param proposal.sigma A numeric vector with the standard deviation of the normal distributions used
#' as proposal distribution. These standard deviations are only used during the first stage of the
#' sampling where the proposals come from a set of independent normal distributions.
#' @param iter_update_par An integer specifying the number of iterations to update the covariance matrix
#' of the proposal distribution.
#' @param cov.corr A logical flag indicating if the covariance matrix needs to be corrected for
#' positive definiteness.
#' @return
#' A list with the hydraulic parameters of the model (includes transmissivity, storage
#' coefficient and radius of influence) and the fitted parameters (includes a and t0)
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @family base functions
#' @export
#' @references
#'  Bai. An Adaptive Directional Metropolis-within-Gibbs algorithm. Technical Report,
#'  Department of Statistics, University of Toronto, 2009. \url{http://www.utstat.toronto.edu/wordpress/WSFiles/technicalreports/0903.pdf}
#'
#'  Roberts, G. O. & Rosenthal, J. S. Examples of adaptive MCMC Journal of Computational
#'  and Graphical Statistics, 2009, 18, 349-367.
#'
#'  Christen, J. A. & Fox, C. A general purpose sampling algorithm for continuous distributions (the t-walk).
#'  Bayesian Analysis, 2010, 5, 2, 263-281.
#' @examples
#' \dontrun{
#' # Bayesian estimation of hydraulic parameters for a confined aquifer
#' data("theis")
#' # Create a pumping_test object
#' ptest.theis <- pumping_test("Well1", Q = 0.01388, r = 250,
#'                             t = theis$t, s = theis$s)
#' # Define prior distributions of hydraulic parameters
#' prior.pdf <- rep("unif", 3)
#' prior.parameters <- matrix(0, 3, 2)
#' # # Transmissivity PDF par
#' prior.parameters[1, ] <- c(1e-5, 1e-2)
#' # Storage coefficient PDF par
#' prior.parameters[2, ] <- c(1e-7, 1e-3)
#' # Std deviation of residuals PDF par
#' prior.parameters[3, ] <- c(1e-05, 0.1)
#' proposal.sigma <- c(0.1, 5, 0.001)
#' # Parameter estimation using AMCMC
#' ptest.theis.mcmc <- fit.sampling(ptest = ptest.theis,
#'                     model = "theis", method = "amcmc",
#'                     prior.pdf = prior.pdf,
#'                     prior.parameters = prior.parameters,
#'                     iterations = 30000, burnIn = 25000, seed = 12355,
#'                     proposal.sigma = proposal.sigma,
#'                     iter_update_par = 300, cov.corr = TRUE)
#' # Assign parameters
#' hydraulic.parameters(ptest.theis) <-
#' ptest.theis.mcmc$hydraulic.parameters[25000:30000,]
#' hydraulic.parameter.names(ptest.theis) <-
#' ptest.theis.mcmc$hydraulic.parameters.names
#' estimated(ptest.theis) <- TRUE
#' # Plot Results
#' cex <- 1.0
#' plot(ptest.theis, type = "uncertainty", cex = cex, cex.axis = cex,
#'      cex.lab = cex, cex.main = cex)
#' # Bayesian estimation using twalk
#' ptest.theis1 <- ptest.theis
#' ptest.theis.mcmc1 <- fit.sampling(ptest = ptest.theis1,
#'                                  model = "theis", method = "twalk",
#'                                  prior.pdf = prior.pdf,
#'                                  prior.parameters = prior.parameters,
#'                                  iterations = 30000)
#' # Assign parameters
#' hydraulic.parameters(ptest.theis1) <-
#' ptest.theis.mcmc1$hydraulic.parameters[25000:30000,]
#' hydraulic.parameter.names(ptest.theis1) <-
#' ptest.theis.mcmc1$hydraulic.parameters.names
#' estimated(ptest.theis1) <- TRUE
#' # Plot Results
#' cex <- 1.0
#' plot(ptest.theis1, type = "uncertainty", cex = cex, cex.axis = cex,
#'      cex.lab = cex, cex.main = cex)
#'
#' # Bayesian estimation of hydraulic parameters for Phreatic aquifer
#' data("boulton")
#' ptest.boulton <- pumping_test("Well1", Q = 0.03, r = 20, t = boulton$t, s = boulton$s)
#' ptest.boulton.fit <- fit(ptest.boulton, 'boulton')
#' prior.pdf <- rep("unif", 5)
#' prior.parameters <- matrix(0, 5, 2)
#' # Transmissivity
#' prior.parameters[1, ] <- c(1e-5, 1e-1)
#' # Storage coefficient
#' prior.parameters[2, ] <- c(1e-6, 1e-1)
#' # Omegad (drainage porosity)
#' prior.parameters[3,] <- c(1e-4, 1e-1)
#' # phi
#' prior.parameters[4, ] <- c(1e-05, 0.1)
#' # Sigma
#' prior.parameters[5, ] <- c(1e-05, 0.1)
#' proposal.sigma <- c(0.1, 5, 5, 0.001, 0.001)
#' ptest.boulton.mcmc <- fit.sampling(ptest = ptest.boulton, model = 'boulton', method = 'amcmc',
#'                                   prior.pdf = prior.pdf, prior.parameters = prior.parameters,
#'                                   iterations = 10000, burnIn = 9000, seed = 12345,
#'                                   proposal.sigma = proposal.sigma,
#'                                   iter_update_par = 200, cov.corr = TRUE)
#' # Assign results
#' hydraulic.parameters(ptest.boulton) <-
#' ptest.boulton.mcmc$hydraulic.parameters[9000:10000,]
#' hydraulic.parameter.names(ptest.boulton) <-
#' ptest.boulton.mcmc$hydraulic.parameters.names
#' estimated(ptest.boulton) <- TRUE
#' # Plot Results
#' cex <- 1.0
#' plot(ptest.boulton, type = "uncertainty", cex = cex, cex.axis = cex,
#'      cex.lab = cex, cex.main = cex)
#'}
fit.sampling <- function(ptest, model, method = 'amcmc', prior.pdf, prior.parameters,
                         iterations = 10000,  burnIn = 0.8*iterations,
                         seed = 12345, proposal.sigma=NULL, iter_update_par = 100,
                         cov.corr = FALSE){
  if(class(ptest) != 'pumping_test'){
    stop('A pumping_test object is required as input')
  }
  set.seed(seed)
  ptest.fit <- fit(ptest, model)
  fn_sol0 <- paste(model, '_solution_initial', sep = '')
  args <- list(ptest = ptest)
  s0 <- do.call(fn_sol0, args)
  startvalue <- as.numeric(s0)
  nfitpar <- length(startvalue)
  fit.parameter.names <- names(s0)
  hydr.parameter.names <- model.parameters(model)
  nhydrpar <- length(hydr.parameter.names)
  #print(fit.parameter.names)
  #print(hydr.parameter.names)
  #print(nhydrpar)
  # Define initial value for variance
  startvalue[(nfitpar+1)] <- 0.1
  # Define prior distributions of self-starter function parameters
  prior.parameters1 <- matrix(0.0, (nfitpar+1), 2)
  fn_par0 <- paste0(model, '_calculate_parameters')
  #print(class(prior.parameters))
  #print(prior.parameters[1:nhydrpar,1])
  par1 <- as.list(prior.parameters[1:nhydrpar,1])
  names(par1) <- hydr.parameter.names
  par2 <- as.list(prior.parameters[1:nhydrpar,2])
  names(par2) <- hydr.parameter.names
  args1 <- list(ptest = ptest, par = par1, hydraulic = F)
  args2 <- list(ptest = ptest, par = par2, hydraulic = F)
  p1 <- do.call(fn_par0, args = args1)
  p2 <- do.call(fn_par0, args = args2)
  p1 <- unlist(p1)
  p2 <- unlist(p2)
  p3 <- cbind(p1, p2)
  #print(p3)
  mn <- unname(apply(p3, 1, min))
  mx <- unname(apply(p3, 1, max))
  prior.parameters1[1:nhydrpar,1] <- mn
  prior.parameters1[1:nhydrpar,2] <- mx
  prior.parameters1[(nhydrpar+1),] <- prior.parameters[(nfitpar+1),]
  # Check if initial value belongs to domain prior PDFS
  inside <- vector("logical", length = (nhydrpar-1))
  for(ipar in 1:nhydrpar){
    if(startvalue[ipar] > prior.parameters1[ipar,1] &
       startvalue[ipar] <= prior.parameters1[ipar,2]){
      inside[ipar] <- TRUE
    }
  }
  if(sum(inside) != nhydrpar){
    #print(prior.parameters1)
    startvalue <- apply(prior.parameters1, 1, mean)
    #print(startvalue)
  }
  #
  chain <- NULL
  acceptance <- 0.0
  if(method == 'amcmc'){
    # run the adaptative MCMC
    chain <- run_adap_metropolis_MCMC(startvalue, iterations, iter_update_par, ptest,
                                      model, prior.pdf, prior.parameters1,
                                      proposal.sigma,
                                      cov.corr)
    # Calculate acceptance rate
    acceptance <-  1.0-mean(duplicated(chain[-(1:burnIn),]))
  }
  else if(method == 'twalk'){
    # run the traverse-walk
    startvalue <- apply(prior.parameters1, 1, mean)
    startvalue2 <- 1.01*startvalue
    #print(startvalue)
    #print(startvalue2)
    restwalk <- run_twalk_MCMC(startvalue1 = startvalue, startvalue2 = startvalue2,
                               ptest = ptest, model = model, prior.pdf = prior.pdf,
                               prior.parameters = prior.parameters1,
                               iterations = iterations)
    acceptance <- restwalk$acceptance
    chain <- restwalk$xxp[,1:(nfitpar+1)]
  }
  else{
    stop('ERROR: The requested method is not implemented')
  }
  # Calculate hydraulic parameters
  fit.par <- list()
  pos <- 0
  for(ipar in fit.parameter.names){
    #print(ipar)
    pos <- pos + 1
    current_cmd <- paste('fit.par$', ipar, '<-chain[,',as.character(pos),']', sep = '')
    #print(current_cmd)
    eval(parse(text = current_cmd))
  }
  fn <- paste(model, '_calculate_parameters', sep = '')
  #print(fn)
  args <- list(ptest = ptest, par = fit.par)
  hydr.par <- do.call(fn, args)
  #print(names(hydr.par))
  #print(hydr.par$Tr)
  npar1 <- length(names(hydr.par))
  hydr.par.values <- matrix(0.0, nrow = (iterations+1), ncol = npar1)
  for(ipar in 1:npar1){
    #print(length(hydr.par[[ipar]]))
    hydr.par.values[,ipar] <- hydr.par[[ipar]]
  }
  hydr.par.names <- names(hydr.par)
  #print(hydr.par.names)
  # Create list with results
  res <- list(chain = chain, acceptance = acceptance, initial = startvalue,
              final = chain[iterations,], prior.pdf = prior.pdf,
              prior.parameters = prior.parameters, burnIn = burnIn,
              hydraulic.parameters = hydr.par.values,
              hydraulic.parameters.names = hydr.par.names)
  return(res)
}
#' @title
#' model.parameters
#' @description
#' This function returns a list with the hydraulic or fit parameters of a specific analytical
#' solution.
#' @param model A character string specifying the analytical model.
#' @param hydraulic A logical flag indicating the type of parameters to return (TRUE for
#' hydraulic, FALSE for fit parameters)
#' @return
#' A list with the corresponding hydraulic or fit parameters.
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @family base functions
#' @export
#' @examples
#' model <- 'theis'
#' # Hydraulic parameters of the Theis model
#' print(model.parameters(model))
#' # Fit parameters of the Theis model
#' print(model.parameters(model, FALSE))
model.parameters <- function(model, hydraulic = TRUE){
  model.par <- list("theis" = c("a", "t0"),
                    "cooper_jacob" = c("a", "t0"),
                    "hantush_jacob" = c("a", "t0", "lambda"),
                    "boulton" = c("a", "t0", "t1", "phi"),
                    "agarwal_recovery" = c("a", "t0"),
                    "general_radial_flow" = c("a", "t0", "n"),
                    "papadopulos_cooper" = c("a", "t0"),
                    "warren_root" = c("a", "t0", "t1", "tm"),
                    "gringarten" = c("a", "t0"),
                    "cooper" = c("cd", "t0")
  )
  #
  base.par <- c("Tr", "Ss")
  hydr.par <- list("theis" = base.par,
                   "cooper_jacob" = c(base.par),
                   "hantush_jacob" = c(base.par, "Ka"),
                   "boulton" = c(base.par, "omegad", "phi"),
                   "general_radial_flow" = c(base.par, "n"),
                   "papadopulos_cooper" = c(base.par, "cd"),
                   "agarwal_skin" = c(base.par, "cd"),
                   "agarwal_recovery" = c(base.par),
                   "warren_root" = c(base.par),
                   "gringarten" = c(base.par),
                   "cooper" = c(base.par))
  if(hydraulic){
    return(hydr.par[[model]])
  }
  else{
    return(model.par[[model]])
  }
}
#' @title
#' evaluate
#' @description
#' Function to calculate the drawdown in time or space from a set of parameters determined using
#' the \code{fit} function.
#' @param ptest A pumping_test object
#' @param model Character string specifying the model used to calculate the drawdown
#' @param test.time A logical flag indicating if the evaluation of the pumping test data needs to
#' be calculated at the drawdown times
#' @param n Number of points used to calculate the drawdown if test.time = FALSE
#' @return
#' A list with the entries
#' \itemize{
#' \item t: Vector with the time
#' \item s: Vector with the calculated drawdown
#' \item dsdlnt: Vector with the derivative of drawdown with respect to the log of time
#' }
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @family base functions
#' @export
#' @examples
#' data(theis)
#' ptest <- pumping_test("Well1", Q = 1.3e-3, r = 200, t = theis$t, s = theis$s)
#' ptest.fit <- fit(ptest, "theis")
#' hydraulic.parameters(ptest) <- ptest.fit$hydraulic_parameters
#' fit.parameters(ptest) <- ptest.fit$parameters
#' model(ptest) <- "theis"
#' estimated(ptest) <- TRUE
#' res <- evaluate(ptest, "theis")
evaluate <- function(ptest, model, test.time = TRUE, n = 50){
  if(class(ptest) != 'pumping_test'){
    stop('A pumping_test object is expected as input')
  }
  if(!ptest$estimated){
    stop('An estimated pumping_test is required')
  }
  #
  hydraulic_parameters <- ptest$hydraulic_parameters
  parameters <- ptest$parameters
  parameter_names <- names(parameters)
  #
  t <- ptest$t
  if(!test.time){
    trange <- range(t)
    t <- logseq(log10(trange[1]), log10(trange[2]), n)
  }
  #
  sd1 <- dsdlnt <- NULL
  #
  current_fn <- 'sd1<-'%s+%model%s+%'_solution(ptest,'
  current_dfn <- 'dsdlnt<-'%s+%model%s+%'_solution_dlogt(ptest,'
  #fn_par <- list(ptest = ptest)
  for(ipar in parameter_names){
    current_fn <- current_fn%s+%ipar%s+%','
    current_dfn <- current_dfn%s+%ipar%s+%','
    current_par <- ipar%s+%'<-parameters$'%s+%ipar
    eval(parse(text = current_par))
    #print(current_par)
  }
  current_fn <- current_fn%s+%'t)'
  current_dfn <- current_dfn%s+%'t)'
  #print(current_dfn)
  eval(parse(text = current_fn))
  eval(parse(text = current_dfn))
  #print(dsdlnt)
  res <- list(t = t, sd = sd1, dsdlnt = dsdlnt)
  return(res)
}
#' @title
#' simulate
#' @description
#' Function to generate a drawdown curve using a given model with given hydraulic parameters
#' and noise
#' @param model Model used to generate the drawdown curve
#' @param Q Pumping rate in m3/s
#' @param r Distance from pumping well to observation well in m
#' @param t Numeric vector with the values of time to simulate the drawdown
#' @param par List with parameters of the model (Transmissivity, Storage Coefficient, etc)
#' @param add.par List with the additional parameter required to define the pumping_test object
#' @param sigma Standard deviation of the noise added to the simulated drawdown.
#' @param seed Seed for the random number generator
#' @return A pumping_test object with the simulated drawdown
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @family base functions
#' @export
#' @examples
#'  par <- list(Tr = 1.4e-3, Ss = 2.1e-5)
#'  t <- logseq(2, 5, 50)
#'  ptest.sim <- simulate('theis', Q = 1.388e-2, r =250, t = t, par)
simulate <- function(model, Q, r, t, par, add.par = NULL, sigma = 0.01, seed = 12345){
  sd1 <- vector('numeric', length = length(t))
  ptest <- pumping_test('W1', Q = Q, r = r, t = t, s = sd1)
  if(!is.null(add.par)){
    ptest$additional_parameters <- add.par
  }
  fit.par <- NULL
  fn <- paste('fit.par<-', model, '_calculate_parameters(ptest, par, FALSE)', sep = '')
  eval(parse(text = fn))
  par.names <- names(fit.par)
  fn1 <- paste('sd1 <- ', model, '_solution(ptest,', sep='')
  for(ipar in par.names){
    current_par <- paste('fit.par$', ipar, ',', sep='')
    fn1 <- paste(fn1, current_par, sep = '')
  }
  fn1 <- paste(fn1,'ptest$t)', sep = '')
  eval(parse(text = fn1))
  set.seed(seed)
  sd1 <- sd1*(1 + rnorm(length(t), mean = 0, sd = sigma))
  pos <- sd1 > 1e-2
  ptest <- pumping_test('SimulatedCurve', Q = Q, r = r, t = t, s = sd1)
  return(ptest)
}
#' @title
#' confint.pumping_test
#' @description
#' Calculates the confindence intervals of the hydraulic parameters, drawdown and
#' its derivative by bootstrapping the residuals
#' @param object An estimated pumping_test object
#' @param parm Name of parameters to be analyzed
#' @param level A numeric vector with the probabilities used to calculate the confidence
#' intervals.
#' @param ... additional argument(s) for methods.
#' @param slug Logical flag to indicate a slug test.
#' @param method A character string with the method to estimate the confidence intervals.
#' The methods currently supported are:
#' \itemize{
#' \item wald
#' \item bootstrap
#' \item jackniffe
#' }.
#' @param d Number of points used in the calculation of the derivative.
#' @param neval Number of bootstrap realizations used in the estimation of confidence
#' intervals.
#' @param seed Random seed.
#' @return
#' In the case of method=wald the function returns a list with the following entries:
#' \itemize{
#' \item hydraulic.parameters.names: A character vector with the names of the hydraulic parameters
#' \item hydraulic.parameters.ci: Confidence intervals of the hydraulic parameters
#' }
#' In the case of method=bootstrap the function returns a list with the following entries:
#' \itemize{
#' \item fit.parameters: A matrix with the values of the parameters of the nonlinear
#' regression model
#' \item total.drawdown: A matrix with the values of the drawdown
#' \item total_dsdlogt: A matrix with the values of the derivative of drawdown
#' \item sci: A matrix with the confidence intervals of drawdown
#' \item dsci: A matrix with the confidence intervals of the derivative of drawdown
#' \item hydraulic.parameters: A matrix with the hydraulic parameter estimated from the boostrapp
#' realizations
#' \item parameter.names: A vector with the name of the fitting parameters
#' \item hydraulic.parameters.names: A character vector with the names of the hydraulic parameters
#' \item hydraulic.parameters.ci: Confidence intervals of the hydraulic parameters
#' }
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @family base functions
#' @importFrom stats qt
#' @importFrom stats vcov
#' @importFrom stats coef
#' @export
#' @examples
#' data(theis)
#' ptest.theis <- pumping_test("Well1", Q = 1.388e-2, r = 250, t= theis$t, s = theis$s)
#' ptest.theis.fit <- fit(ptest.theis, "theis")
#' hydraulic.parameters(ptest.theis) <- ptest.theis.fit$hydraulic_parameters
#' fit.parameters(ptest.theis) <- ptest.theis.fit$parameters
#' model(ptest.theis) <- "theis"
#' estimated(ptest.theis) <- TRUE
#' ptest.theis.ci1 <- confint(ptest.theis, level = 0.975)
#' ptest.theis.ci2 <- confint(ptest.theis, level = c(0.025, 0.975), method = 'bootstrap',
#'  d = 20, neval = 100)
confint.pumping_test <- function(object, parm, level = 0.95, ..., slug = FALSE,
                                 method = 'wald', d = 20, neval = 100,
                                 seed = 12345){
  ptest <- object
  if(class(ptest) != "pumping_test"){
    stop('A pumping_test object is required as input')
  }
  if(!ptest$estimated){
    stop('An estimated pumping_test is required as input')
  }
  res <- NULL
  if(method == 'wald'){
    res <- confint_wald(ptest, level)
    return(res)
  }
  if (method == 'bootstrap'){
    res <- confint_bootstrap(ptest, level = level, slug = slug, d = d, neval = neval,
                             seed = seed)
  }
  if(method == 'jackniffe'){
    res <- confint_jackniffe(ptest, level, slug = slug)
  }
  return(res)
}
#' @title
#' confint_wald
#' @description
#' Function to calculate the Wald's confidence intervals of the hydraulic parameters
#' @param ptest A pumping_test object
#' @param level A numeric vector with the confidence levels required
#' @return
#' A list with the following entries:
#' \itemize{
#' \item method: A character string specifying the method used in the CI estimation (wald)
#' \item hydraulic.parameters.names: A character vector with the names of the hydraulic
#' parameters.
#' \item hydraulic.parameters.ci: A matrix with the requested confidence of the hydraulic
#' parameters.
#' }
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @family base functions
#' @importFrom stats vcov qt
#' @export
confint_wald <- function(ptest, level = 0.95){
  if(class(ptest) != 'pumping_test'){
    stop('ERROR a pumping_test object is required as input')
  }
  if(length(level) != 1){
    stop('ERROR too many confidence levels requested fro the analysis')
  }
  model <- ptest$model
  res.fit <- fit(ptest, model)
  fn <- paste0(model, '_calculate_parameters')
  parameters <- unname(ptest$parameters)
  npars <- length(parameters)
  nhpars <- length(ptest$hydraulic_parameters)
  fit.pars <- ptest$parameters
  par_names <- names(fit.pars)
  # var-covar matrix of the parameters
  cov.pars <- vcov(res.fit$resfit)
  # std error of the parameters
  se.pars <- sqrt(diag(cov.pars))
  # Define degrees of freedom
  df <- length(ptest$t)-npars
  # Define significance level
  p <- max(level)
  # Calculate t.value
  t.val <- qt(p, df)
  # Initialize the CI matrices
  pars.ci <- matrix(0.0, nrow = npars, ncol = 2)
  hydr.par.ci <- matrix(0.0, nrow = (nhpars), ncol = 2)
  # Actual calculation of CI
  par.low <- ptest$parameters
  par.upp <- ptest$parameters
  for(i in 1:npars){
    current.ci <- c(parameters[[i]]-t.val*se.pars[i],
                    parameters[[i]]+t.val*se.pars[i])
    pars.ci[i,] <- current.ci
    par.low[[i]] <- pars.ci[i,1]
    par.upp[[i]] <- pars.ci[i,2]
  }
  hydrpar.low <- do.call(fn, args = list(ptest = ptest, par = par.low,
                                         hydraulic = TRUE))
  hydrpar.upp <- do.call(fn, args = list(ptest = ptest, par = par.upp,
                                         hydraulic = TRUE))
  low <- matrix(unlist(hydrpar.low), nrow = (nhpars), ncol = 1)
  upp <- matrix(unlist(hydrpar.upp), nrow = (nhpars), ncol = 1)
  hydr.par.ci[,1] <- low
  hydr.par.ci[,2] <- upp
  low <- apply(hydr.par.ci, 1, min)
  upp <- apply(hydr.par.ci, 1, max)
  hydr.par.ci[,1] <- as.numeric(format(low, digits = 5, scientific = TRUE))
  hydr.par.ci[,2] <- as.numeric(format(upp, digits = 5, scientific = TRUE))
  level1 <- as.character(c(1-level, level)*100)
  for(i in 1:2){
    level1[i] <- paste0(level1[i],'%')
  }
  row.names(hydr.par.ci) <- names(ptest$hydraulic_parameters)
  colnames(hydr.par.ci) <- level1
  #
  res <- list(method = 'wald', hydraulic.parameters.names = names(ptest$hydraulic_parameters),
              hydraulic.parameters.ci = hydr.par.ci)
  return(res)
}
#' @title
#' confint_bootstrap
#' @description
#' Function to calculate the confidence intervals of the hydraulic parameters using
#' bootstrapping.
#' @param ptest A pumping_test object
#' @param level A numeric vector with the significance levels of the confidence intervals.
#' @param slug A logical flag indicating if a slug test is used in the analysis.
#' @param d An integer value specifying the number of points used in the calculation of the
#' derivative.
#' @param neval An integer value specifying the number of realizations of the resampling procedure.
#' @param seed random seed
#' @return
#' This function returns a list with the following entries:
#' \itemize{
#' \item fit.parameters: A matrix with the values of the parameters of the nonlinear
#' regression model
#' \item method: A character string specifying the method used in the CI estimation (bootstrap)
#' \item total.drawdown: A matrix with the values of the drawdown
#' \item total_dsdlogt: A matrix with the values of the derivative of drawdown
#' \item sci: A matrix with the confidence intervals of drawdown
#' \item dsci: A matrix with the confidence intervals of the derivative of drawdown
#' \item hydraulic.parameters: A matrix with the hydraulic parameter estimated from the boostrapp
#' realizations
#' \item parameter.names: A vector with the name of the fitting parameters
#' \item hydraulic.parameters.names: A character vector with the names of the hydraulic parameters
#' \item hydraulic.parameters.ci: Confidence intervals of the hydraulic parameters
#' }
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @family base functions
#' @importFrom stats quantile
#' @importFrom minpack.lm nlsLM
#' @export
confint_bootstrap <- function(ptest, level = c(0.025, 0.975), slug = F, d = 20,
                              neval = 100, seed = 12345){
  if(class(ptest) != 'pumping_test'){
    stop('ERROR a pumping_test object is required as input')
  }
  model <- ptest$model
  # Call function to generate values of parameters: a and t0
  # depending on model used in estimation
  shat <- NULL
  initial_sol <- NULL
  fn <- paste('shat <- ',model, '_solution(ptest,', sep='')
  fn_parameters <- paste(model,'_calculate_parameters',sep='')
  model_formula <- paste('s~', model, '_solution(ptest,', sep='')
  fit.pars <- ptest$parameters
  par_names <- names(fit.pars)
  for(ipar in par_names){
    fn <- paste(fn, ipar, ',')
    model_formula <- paste(model_formula, ipar, ',', sep = '')
    current_par <- paste(ipar, '<-ptest$parameters$', ipar, sep ='')
    eval(parse(text = current_par))
  }
  model_formula <- paste(model_formula, 't)', sep='')
  model_formula <- as.formula(model_formula)
  #print(model_formula)
  fn <- paste(fn,' ptest$t)', sep = '')
  eval(parse(text = fn))
  # Calculate initial solution
  initial_fn <- paste('initial_sol<-', model, '_solution_initial(ptest)', sep = '')
  eval(parse(text = initial_fn))
  # calculate residuals from original fit
  residuals <- shat - ptest$s
  #print(residuals)
  nres <- length(residuals)
  nt <- length(ptest$t)
  set.seed(seed)
  if(slug){
    s_dlogt <- log_derivative(ptest$t, -ptest$s, d = d, method = 'spline')
  } else {
    s_dlogt <- log_derivative(ptest$t, ptest$s, d = d, method = 'spline')
  }
  n_dlogt <- length(s_dlogt$y)
  npar <- length(names(initial_sol))
  fit.parameters <- matrix(0.0, nrow = neval, ncol = npar)
  total_drawdown <- matrix(0.0, nrow = neval, ncol = nres)
  total_dsdlogt <- matrix(0.0, nrow = neval, ncol = n_dlogt)
  hydr.parameters <- NULL
  nhydrpar <- -1
  lower <- vector('numeric', length = npar)
  lower[1:npar] <- 1.0e-15
  for(ieval in 1:neval){
    #print(ieval)
    # resample residuals
    current_residuals <- sample(residuals, size = nres, replace = TRUE)
    # define current drawdown
    current_drawdown <- shat + current_residuals
    #print(c(min(current_drawdown),max(current_drawdown)))
    #print(c('s',length(current_drawdown)))
    total_drawdown[ieval,] <- current_drawdown
    current_drawdown_dlogt <- log_derivative(ptest$t, current_drawdown,
                                             method = 'spline', d = d)
    current_ndlogt <- length(current_drawdown_dlogt$y)
    if(current_ndlogt == n_dlogt){
      total_dsdlogt[ieval,] <- current_drawdown_dlogt$y
    } else {
      total_dsdlogt[ieval,] <- NA
    }
    #print(length(current_drawdown_dlogt$y))
    
    current_df <- data.frame(t = ptest$t, s = current_drawdown)
    model.nls1 <- nlsLM(model_formula, data = current_df,
                        start = initial_sol, trace = F,
                        control = nls.lm.control(ftol = 1e-10,
                                                 ptol = 1e-10,
                                                 maxiter = 200,
                                                 nprint = -1,
                                                 maxfev = 1000),
                        lower = lower)
    #
    fit.parameters[ieval,] <- coef(model.nls1)
    current_par <- as.list(coef(model.nls1))
    current_hydr.parameters <- do.call(fn_parameters, list(ptest = ptest, par = current_par))
    if(ieval == 1){
      nhydrpar <- length(names(current_hydr.parameters))
      hydr.parameters <- matrix(0,0, nrow = neval, ncol = nhydrpar)
    }
    for(ipar in 1:nhydrpar){
      hydr.parameters[ieval,ipar] <- current_hydr.parameters[[ipar]]
    }
  }
  #
  nlevel <- length(level)
  ssvalues_ci <- matrix(0.0, nrow = nt, ncol = nlevel)
  dsvalues_ci <- matrix(0.0, nrow = n_dlogt, ncol = nlevel)
  for(it in 1:nt){
    current_drawdown_ci <- quantile(total_drawdown[,it],
                                    probs = level,
                                    na.rm = TRUE)
    ssvalues_ci[it,1:nlevel] <- unname(current_drawdown_ci)
  }
  for(it in 1:n_dlogt){
    current_dsdlogt <- total_dsdlogt[,it]
    pos_valid <- !is.na(current_dsdlogt)
    if(sum(pos_valid) > 0){
      current_dsdlogt <- current_dsdlogt[pos_valid]
    }
    current_ddrawdown_ci <- quantile(current_dsdlogt,
                                     probs = level,
                                     na.rm = TRUE)
    dsvalues_ci[it,] <- unname(current_ddrawdown_ci)
  }
  hydr.par.ci <- matrix(0.0, nrow = nhydrpar, ncol = length(level))
  for(ipar in 1:nhydrpar){
    hydr.par.ci[ipar,] <- quantile(hydr.parameters[,ipar], prob = level,
                                   na.rm = TRUE)
  }
  row.names(hydr.par.ci) <- names(current_hydr.parameters)
  level1 <- as.character(level*100)
  for(ipar in 1:length(level)){
    level1[ipar] <- paste(level1[ipar],'%',sep = '')
  }
  colnames(hydr.par.ci) <- level1
  #
  res <- list(method = 'bootstrap', fit.parameters = fit.parameters,
              total_drawdown = total_drawdown,
              total_dsdlogt = total_dsdlogt, sci = ssvalues_ci, dsci = dsvalues_ci,
              t_dsdlogt = s_dlogt$x, hydraulic.parameters = hydr.parameters,
              parameter.names= par_names,
              hydraulic.parameters.names = names(current_hydr.parameters),
              hydraulic.parameters.ci = hydr.par.ci)
  return(res)
}
#' @title
#' confint_jackniffe
#' @description
#' Function to calculate the
#' @param ptest A pumping_test object
#' @param level A numeric vector with the significance levels of the confidence intervals.
#' @param slug A logical flag indicating if a slug test is used in the analysis.
#' @return
#' A list with the following entries:
#' \itemize{
#' \item method: A character string specifying the method used in the CI estimation (jackniffe)
#' \item hydraulic.parameters.names: A character vector with the names of the hydraulic
#' parameters.
#' \item hydraulic.parameters.ci: A matrix with the requested confidence of the hydraulic
#' parameters.
#' }
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @family base functions
#' @importFrom stats vcov quantile
#' @export
confint_jackniffe <- function(ptest, level = c(0.025, 0.975), slug = F){
  if(class(ptest) != 'pumping_test'){
    stop('ERROR a pumping_test object is required')
  }
  #
  if(!ptest$estimated){
    stop('ERROR the pumping test must be estimated')
  }
  #
  model <- ptest$model
  ptest.fit <- fit(ptest, model)
  cov.fitpar <- vcov(ptest.fit$resfit)
  se.fitpar <- diag(vcov(ptest.fit$resfit))
  rfitpar <- rmvnorm(1000, mean = unlist(ptest.fit$parameters), sigma = cov.fitpar)
  fn <- paste0(model,'_calculate_parameters')
  par <- list()
  for(i in 1:ncol(rfitpar)){
    par[[i]] <- rfitpar[,i]
  }
  names(par) <- names(ptest.fit$parameters)
  args <- list(ptest = ptest, par = par, hydraulic = TRUE)
  hydr.par <- do.call(fn, args = args)
  hydr.par.se <- vector('numeric', length = length(hydr.par))
  for(i in 1:length(hydr.par)){
    hydr.par.se[i] <- sqrt(var(hydr.par[[i]]))
  }
  #
  current.hydrpar <- ptest$hydraulic_parameters
  current.names <- names(current.hydrpar)
  current.hydrpar <- unlist(current.hydrpar)
  #
  current.fitpar <- unlist(ptest.fit$parameters)
  #
  t <- ptest$t
  s <- ptest$s
  ndat <- length(s)
  pos <- seq(1, ndat, by = 1)
  ptest1 <- ptest
  dfbeta <- matrix(0.0, nrow = ndat, ncol = length(current.hydrpar))
  hydr.par.res <- matrix(0.0, nrow = ndat, ncol = length(current.hydrpar))
  for(i in 1:ndat){
    pos_valid <- pos != i
    current_t <- t[pos_valid]
    current_s <- s[pos_valid]
    ptest1 <- pumping_test("Well1", Q = ptest$Q, r = ptest$r,
                           t = current_t, s = current_s)
    ptest1.fit <- fit(ptest1, model)
    current.hydrpar1 <- unlist(ptest1.fit$hydraulic_parameters)
    current.fitpar1 <- unlist(ptest1.fit$parameters)
    #print(current.fitpar1)
    hydr.par.res[i,] <- current.hydrpar1
    #rel.err <- 100*abs(current.hydrpar1-current.hydrpar)/current.hydrpar
    dfbeta[i,] <- abs(current.hydrpar1-current.hydrpar)/hydr.par.se
    #print(dfbeta)
  }
  #
  hydr.par.ci <- matrix(0.0, nrow = ncol(hydr.par.res), ncol = length(level))
  for(ipar in 1:ncol(hydr.par.res)){
    hydr.par.ci[ipar,] <- quantile(hydr.par.res[,ipar], probs = level)
  }
  rownames(hydr.par.ci) <- current.names
  level.names <- vector('character', length = length(level))
  for(ilev in 1:length(level.names)){
    level.names[ilev] <- paste0(as.character(100*level[ilev]),"%")
  }
  colnames(hydr.par.ci) <- level.names
  #
  colnames(dfbeta) <- current.names
  res <- list(method = 'jackniffe', hydraulic.parameters.ci = hydr.par.ci,
              hydraulic.parameters.names = current.names,
              dfbeta = dfbeta)
  return(res)
}
