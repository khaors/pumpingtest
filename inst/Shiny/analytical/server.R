#
# This is the server logic of a Shiny web application. You can run the
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(pumpingtest)
library(rhandsontable)
library(ggplot2)
###############################################################################
##                      Global variables
###############################################################################
# Max file input size :
options(shiny.maxRequestSize=30*1024^2)
#
model.types <- list( "Pumping" =
                       c("None", "theis", "cooper_jacob", "hantush_jacob",
                         "boulton", "general_radial_flow",
                         "papadopulos_cooper", "agarwal_skin",
                         "no_flow_boundary", "constant_head_boundary"), #, "agarwal_recovery","warren_root",
                         #"gringarten"),
                     "Slug" = c("None", "cooper", "neuzil"),
                     "Bailer" = c("None", "mcdonald"))
#
colors <- c('#00BFC4', '#F8766D')
# Define server logic required to draw a histogram
shinyServer(function(input, output) {
  output$uptc.logo <- renderImage(list(src = "uptc_jpg.jpg"),
                                  deleteFile = FALSE)
  #########################################################################################
  #                         Curve Type Explorer Tab
  #########################################################################################
  output$curve1 <- renderUI({
    selectInput(inputId = "curve_model", label = "Model", choices = model.types[[input$type]])
  })
  #
  output$curve15 <- renderUI({
    tmp <- NULL
    #if(is.null(input$curve_model)) return(NULL)
    #if(input$curve_model == "None") return(NULL)
    tmp <- selectInput(inputId = "curve_scale", label = "Scale", 
                       choices = c("LogLog", 'SemiLog'), selected = "LogLog")
  })
  #
  output$curve2 <- renderText({
    tmp <- NULL
    if(is.null(input$curve_model)) return(NULL)
    if(input$curve_model == "None") return(NULL)
    if(input$curve_model == "boulton"){
      if(is.null(input$curve_sigma) | is.null(input$curve_phi)) return(NULL)
    }
    if(input$curve_model == "hantush_jacob"){
      if(is.null(input$curve_beta)) return(NULL)
    }
    if(input$curve_model == "papadopulos_cooper"){
      if(is.null(input$curve_cd) | is.null(input$curve_rho)) return(NULL)
    }
    if(input$curve_model == "general_radial_flow"){
      if(is.null(input$curve_n) | is.null(input$curve_rd)) return(NULL)
    }
    if(input$curve_model == "agarwal_skin"){
      if(is.null(input$curve_cd) | is.null(input$curve_rd) | is.null(input$curve_sigma)) 
        return(NULL)
    }
    #
    # if(input$curve_model == "theis") {
    #   tmp <- paste0(p(HTML("<h3>Theis Model</h3>")),
    #                 p(HTML("<h5>Theis (1935) presented an exact analytical solution
    #                        for the transient drawdown in an  infinite uniform confined
    #                        aquifer.  </h5>")))
    # }
    # #if(input$curve_model == "cooper_jacob")
    # if(input$curve_model == "boulton") {
    #   paste0(p(HTML("<h3>Boulton Model</h3><br>")),
    #          p(HTML("<h5>This model is chevere</h5>")))
    # }
    # if(input$curve_model == "hantush_jacob") {
    #   paste0(p(HTML("<h3>Hantush-Jacob Model</h3><br>")),
    #          p(HTML("<h5>This model is chevere</h5>")))
    # }
    # if(input$curve_model == "papadopulos_cooper") {
    #   paste0(p(HTML("<h3>Papadopulos-Cooper Model</h3><br>")),
    #          p(HTML("<h5>This model is chevere</h5>")))
    # }
    return(tmp)
    })
  #
  output$curve3 <- renderUI({
    tmp <- NULL
    if(is.null(input$curve_model)) return(NULL)
    if(input$curve_model == "boulton"){
      tmp <- textInput(inputId = "curve_phi", label = "Phi= ", value = "1.0")
    }
    if(input$curve_model == "hantush_jacob"){
      tmp <- textInput(inputId = "curve_beta", label = "Beta= ", value = "1.0")
    }
    if(input$curve_model == "papadopulos_cooper"){
      tmp <- textInput(inputId = "curve_cd", label = "Cd (Well Storage)= ", value = "1.0")
    }
    if(input$curve_model == "general_radial_flow"){
      tmp <- textInput(inputId = "curve_n", label = "Flow Dimension= ", value = "1.0")
    }
    if(input$curve_model == 'agarwal_skin'){
      tmp <- textInput(inputId = "curve_cd", label = "Cd (Well Storage)", value = "1.0")
    }
    if(input$curve_model == "warren_root"){
      tmp <- textInput(inputId = "curve_sigma", label = "Sigma= ", value = "1.0")
    }
    return(tmp)
  })
  #
  output$curve4 <- renderUI({
    tmp <- NULL
    if(is.null(input$curve_model)) return(NULL)
    if(input$curve_model == "boulton"){
      tmp <- textInput(inputId = "curve_sigma", label = "Sigma= ", value = "1.0e-3")
    }
    if(input$curve_model == "papadopulos_cooper"){
      tmp <- textInput(inputId = "curve_rho", label = "rho= ", value = "1.0")
    }
    if(input$curve_model == "general_radial_flow"){
      tmp <- textInput(inputId = "curve_rd", label = "rd (Dimensionless radius)= ",
                       value = "1.0")
    }
    if(input$curve_model == "agarwal_skin"){
      tmp <- textInput(inputId = "curve_rd", label = "rd (Dimensionless radius)= ",
                       value = "1.0")
    }
    if(input$curve_model == "warren_root"){
      tmp <- textInput(inputId = "curve_lambda", label = "lambda= ", value = "1.0")
    }
    #
    return(tmp)
  })
  #
  output$curve4.5 <- renderUI({
    tmp <- NULL
    if(is.null(input$curve_model)) return(NULL)
    if(input$curve_model == 'agarwal_skin'){
      tmp <- textInput(inputId = "curve_sigma", label = "sigma= ", value = "1.0")
    }
    return(tmp)
  })
  #
  output$curve5 <- renderUI({
    if(is.null(input$curve_model)) return(NULL)
    tmp <- textInput(inputId = "t_range", label = "t (min,max)", value = "0.5,1e3")
  })
  #
  output$plot_curve <- renderPlot({
    coeffs <- stehfest_coeffs(8)
    if(is.null(input$curve_model)) return(NULL)
    if(input$curve_model == "None") return(NULL)
    if(input$curve_model == "boulton"){
      if(is.null(input$curve_sigma) | is.null(input$curve_phi)) return(NULL)
    }
    if(input$curve_model == "hantush_jacob"){
      if(is.null(input$curve_beta)) return(NULL)
    }
    if(input$curve_model == "papadopulos_cooper"){
      if(is.null(input$curve_cd) | is.null(input$curve_rho)) return(NULL)
    }
    if(input$curve_model == "general_radial_flow"){
      if(is.null(input$curve_n) | is.null(input$curve_rd)) return(NULL)
    }
    t_range <-  as.numeric(unlist(strsplit(input$t_range,",")))
    t_vector <- logseq(from = log10(t_range[1]), to = log10(t_range[2]), n = 100)
    current.model <- input$curve_model
    current.fn <- paste0(current.model,"_well_function")
    current.par <- NULL
    if(current.model == "theis"){
      current.par <- NULL
    }
    if(current.model == "cooper_jacob"){
      current.par <- NULL
    }
    if(current.model == "boulton"){
      current.par <- list(coeffs = coeffs,
                          sigma = as.numeric(input$curve_sigma),
                          phi = as.numeric(input$curve_phi))
    }
    if(input$curve_model == "hantush_jacob"){
      current.par <- list(coeffs = coeffs,
                          beta = as.numeric(input$curve_beta))
    }
    if(input$curve_model == "papadopulos_cooper"){
      current.par <- list(coeffs = coeffs,
                          cd = as.numeric(input$curve_cd),
                          rho = as.numeric(input$curve_rho))
    }
    if(input$curve_model == "general_radial_flow"){
      current.par <- list(coeffs = coeffs,
                          n = as.numeric(input$curve_n),
                          rd = as.numeric(input$curve_rd))
    }
    if(input$curve_model == "agarwal_skin"){
      #print(coeffs)
      if(is.null(input$curve_cd)) return(NULL)
      #print(c(as.numeric(input$curve_cd),as.numeric(input$curve_rd),as.numeric(input$curve_sigma)))
      current.par <- list(coeffs = coeffs, 
                          cd = as.numeric(input$curve_cd), 
                          rd = as.numeric(input$curve_rd),
                          sigma = as.numeric(input$curve_sigma))
      #print(current.par)
    }
    #print(current.fn)
    W <- NULL
    if(current.model != "cooper_jacob" & current.model != "no_flow_boundary" & 
       current.model != "constant_head_boundary"){
      W <- do.call(current.fn, args = list(td = t_vector, par = current.par))
    }
    else if(current.model == "cooper_jacob"){
      W <- do.call(current.fn, args = list(u = 1/t_vector))
    }
    else if(current.model == "no_flow_boundary"){
      rd <- 1.0
      Rd <- 10.0
      td <- t_vector
      term1 <- pracma::expint_E1(rd^2/(4*td))
      term2 <- pracma::expint_E1((rd^2)*(Rd^2)/(4*td))
      W <- 0.5*(term1+term2)#theis_well_function(td1) + theis_well_function(td2)
    }
    else if(current.model == "constant_head_boundary"){
      rd <- 1.0
      Rd <- 20.0
      td <- t_vector
      term1 <- pracma::expint_E1(rd^2/(4*td))
      term2 <- pracma::expint_E1((rd^2)*(Rd^2)/(4*td))
      W <- 0.5*(term1-term2)
    }
    #print(W)
    dW <- log_derivative_central(t_vector, W)
    Wr <- range(W)
    dWr <- range(dW)
    mn <- abs(min(Wr[1], dWr[1]))
    mx <- Wr[2]
    #print(mn)
    #print(mx)
    ctitle <- paste0("Well Function - ", current.model)
    W.df <- data.frame(t = t_vector, W = W)
    dW.df <- data.frame(t = dW$x, dW = dW$y)
    plt <- ggplot() + geom_line(aes(x = t, y = W, colour = colors[1]), data = W.df) + 
      geom_line(aes(x = t, y = dW, colour = colors[2]), data = dW.df) +
      ggtitle(paste0("Diagnostic Plot: ",ctitle)) + 
      scale_color_manual(labels = c("Well Function", "Derivative"), values = colors) +
      #scale_color_manual(values = rev(colors)) +
      theme_bw() + 
      theme(axis.text.x = element_text(size = 14), 
            axis.title.x = element_text(size = 16),
            axis.text.y = element_text(size = 14), 
            axis.title.y = element_text(size = 16),
            legend.title = element_text(size = 14),
            legend.text = element_text(size = 14), 
            plot.title = element_text(size = 18, face = "bold")) +
      scale_x_log10() +
      guides(color=guide_legend("Variables"))
    if(input$curve_scale == "LogLog")
      plt <- plt + scale_y_log10()
    #plot(t_vector, W, type = "l", log = "xy", ylim = c(mn, mx), main = ctitle,
    #     xlab = "1/u")
    #lines(dW$x, dW$y, col = "red")
    #legend("topleft", legend = c("Well function", "Derivative"), col = c("black", "red"),
    #       lty = c(1,1))
    return(plt)
  })
  #
  output$curve6 <- renderUI({
    tmp <- NULL
    if(is.null(input$curve_model)) return(NULL)
    if(input$curve_model == "None") return(NULL)
    if(input$curve_model == "boulton"){
      if(is.null(input$curve_sigma) | is.null(input$curve_phi)) return(NULL)
    }
    if(input$curve_model == "hantush_jacob"){
      if(is.null(input$curve_beta)) return(NULL)
    }
    if(input$curve_model == "papadopulos_cooper"){
      if(is.null(input$curve_cd) | is.null(input$curve_rho)) return(NULL)
    }
    if(input$curve_model == "general_radial_flow"){
      if(is.null(input$curve_n) | is.null(input$curve_rd)) return(NULL)
    }
    # if(input$curve_model == "theis"){
    #   tmp <- withMathJax(
    #     helpText('The equation that describes the variation in the drawdown in a confined 
    #              aquifer is given by: \n $$s(r,t)=\\frac{Q}{4\\pi Tr}W(u)$$ \n where W(u) is 
    #              called the well function and it is defined as:\n 
    #              $$W(u) = \\int_{u}^{+\\infty} \\frac{\\exp{(-u)}}{u}du$$
    #              and u is a dimensionless variable given by:\n
    #              $$u=\\frac{rS}{4\\;Tr\\; t}$$ where
    #              - r: distance between the pumping and 
    #              observation well\n - S: Storage coefficient(dimensionless)<br>
    #              - Tr: Transmisivity (m^2/s)\n - t: Time (in seconds)'))
    #}
    return(tmp)
  })
  #######################################################################################
  #                           Compare Solutions TabSet
  #######################################################################################
  output$compare_time <- renderUI({
    textInput(inputId = "compare_time1", label = "Time(s)", value = "0.5,1e3")
  })
  #
  output$compare_scale <- renderUI({
    tmp <- NULL
    tmp <- selectInput(inputId = "compare_scale1", label = "Scale", 
                       choices = c("LogLog", 'SemiLog'), selected = "LogLog")
  })
  #
  output$compare1a <- renderUI({
    selectInput(inputId = "compare_model1", label = "Model", 
                choices = model.types[[input$type]])
  })
  output$compare1b <- renderUI({
    tmp <- NULL
    if(is.null(input$compare_model1)) return(NULL)
    if(input$compare_model1 == "boulton"){
      tmp <- textInput(inputId = "curve_phi1", label = "Phi= ", value = "1.0")
    }
    if(input$compare_model1 == "hantush_jacob"){
      tmp <- textInput(inputId = "curve_beta1", label = "Beta= ", value = "1.0")
    }
    if(input$compare_model1 == "papadopulos_cooper"){
      tmp <- textInput(inputId = "curve_cd1", label = "Cd (Well Storage)= ", value = "1.0")
    }
    if(input$compare_model1 == "general_radial_flow"){
      tmp <- textInput(inputId = "curve_n1", label = "Flow Dimension= ", value = "1.0")
    }
    if(input$compare_model1 == "warren_root"){
      tmp <- textInput(inputId = "curve_sigma1", label = "Sigma= ", value = "1.0")
    }
    return(tmp)
  })
  #
  output$compare1c <- renderUI({
    tmp <- NULL
    if(is.null(input$compare_model1)) return(NULL)
    if(input$compare_model1 == "boulton"){
      tmp <- textInput(inputId = "curve_sigma1", label = "Sigma= ", value = "1.0e-3")
    }
    if(input$compare_model1 == "papadopulos_cooper"){
      tmp <- textInput(inputId = "curve_rho1", label = "rho= ", value = "1.0")
    }
    if(input$compare_model1 == "general_radial_flow"){
      tmp <- textInput(inputId = "curve_rd1", label = "rd (Dimensionless radius)= ",
                       value = "1.0")
    }
    if(input$compare_model1 == "warren_root"){
      tmp <- textInput(inputId = "curve_lambda1", label = "lambda= ", value = "1.0")
    }
    #
    return(tmp)
    
  })
  #
  output$compare2a <- renderUI({
    selectInput(inputId = "compare_model2", label = "Model", choices = model.types[[input$type]])
  })
  output$compare2b <- renderUI({
    tmp <- NULL
    if(is.null(input$compare_model2)) return(NULL)
    if(input$compare_model2 == "boulton"){
      tmp <- textInput(inputId = "curve_phi2", label = "Phi= ", value = "1.0")
    }
    if(input$compare_model2 == "hantush_jacob"){
      tmp <- textInput(inputId = "curve_beta2", label = "Beta= ", value = "1.0")
    }
    if(input$compare_model2 == "papadopulos_cooper"){
      tmp <- textInput(inputId = "curve_cd2", label = "Cd (Well Storage)= ", value = "1.0")
    }
    if(input$compare_model2 == "general_radial_flow"){
      tmp <- textInput(inputId = "curve_n2", label = "Flow Dimension= ", value = "1.0")
    }
    if(input$compare_model2 == "warren_root"){
      tmp <- textInput(inputId = "curve_sigma2", label = "Sigma= ", value = "1.0")
    }
    return(tmp)
  })
  #
  output$compare2c <- renderUI({
    tmp <- NULL
    if(is.null(input$compare_model2)) return(NULL)
    if(input$compare_model2 == "boulton"){
      tmp <- textInput(inputId = "curve_sigma2", label = "Sigma= ", value = "1.0e-3")
    }
    if(input$compare_model2 == "papadopulos_cooper"){
      tmp <- textInput(inputId = "curve_rho2", label = "rho= ", value = "1.0")
    }
    if(input$compare_model2 == "general_radial_flow"){
      tmp <- textInput(inputId = "curve_rd2", label = "rd (Dimensionless radius)= ",
                       value = "1.0")
    }
    if(input$compare_model2 == "warren_root"){
      tmp <- textInput(inputId = "curve_lambda2", label = "lambda= ", value = "1.0")
    }
    #
    return(tmp)
  })
  #
  output$plot_comparison <- renderPlot({
    coeffs <- stehfest_coeffs(8)
    #print(coeffs)
    #print(input$compare_model1)
    #
    # Model 1
    #
    if(is.null(input$compare_model1)) return(NULL)
    if(input$compare_model1 == "None") return(NULL)
    if(input$compare_model1 == "boulton"){
      if(is.null(input$curve_sigma1) | is.null(input$curve_phi1)) return(NULL)
    }
    if(input$compare_model1 == "hantush_jacob"){
      if(is.null(input$curve_beta1)) return(NULL)
    }
    if(input$compare_model1 == "papadopulos_cooper"){
      if(is.null(input$curve_cd1) | is.null(input$curve_rho1)) return(NULL)
    }
    if(input$compare_model1 == "general_radial_flow"){
      if(is.null(input$curve_n1) | is.null(input$curve_rd1)) return(NULL)
    }
    #
    # Model 2
    #
    if(is.null(input$compare_model2)) return(NULL)
    if(input$compare_model2 == "None") return(NULL)
    
    if(input$compare_model2 == "boulton"){
      if(is.null(input$curve_sigma2) | is.null(input$curve_phi2)) return(NULL)
    }
    if(input$compare_model2 == "hantush_jacob"){
      if(is.null(input$curve_beta2)) return(NULL)
    }
    if(input$compare_model2 == "papadopulos_cooper"){
      if(is.null(input$curve_cd2) | is.null(input$curve_rho2)) return(NULL)
    }
    if(input$compare_model2 == "general_radial_flow"){
      if(is.null(input$curve_n2) | is.null(input$curve_rd2)) return(NULL)
    }
    if(input$compare_model2 == "agarwal_skin"){
      if(is.null(input$curve_cd2) | is.null(input$curve_rd2)) return(NULL)
    }
    #
    #
    #
    t_range <-  as.numeric(unlist(strsplit(input$compare_time1,",")))
    t_vector <- logseq(from = log10(t_range[1]), to = log10(t_range[2]), n = 100)
    #
    current.par1 <- NULL
    current.par2 <- NULL
    current.model1 <- input$compare_model1
    current.model2 <- input$compare_model2
    print(current.model1)
    current.fn1 <- paste0(current.model1,"_well_function")
    current.fn2 <- paste0(current.model2,"_well_function")
    #
    # Define parameters for model 1
    #
    if(current.model1 == "theis" ){
      current.par1 <- NULL
    }
    if(current.model1 == "cooper_jacob"){
      current.par1 <- NULL
    }
    if(current.model1 == "boulton"){
      current.par1 <- list(coeffs = coeffs,
                          sigma = as.numeric(input$curve_sigma1),
                          phi = as.numeric(input$curve_phi1))
    }
    if(current.model1 == "hantush_jacob"){
      current.par1 <- list(coeffs = coeffs,
                          beta = as.numeric(input$curve_beta1))
    }
    if(current.model1 == "papadopulos_cooper"){
      current.par1 <- list(coeffs = coeffs,
                          cd = as.numeric(input$curve_cd1),
                          rho = as.numeric(input$curve_rho1))
    }
    if(current.model1 == "general_radial_flow"){
      current.par1 <- list(coeffs = coeffs,
                          n = as.numeric(input$curve_n1),
                          rd = as.numeric(input$curve_rd1))
    }
    if(current.model1== "agarwal_skin"){
      #print(coeffs)
      if(is.null(input$curve_cd)) return(NULL)
      current.par1 <- list(coeffs = coeffs, 
                          cd = as.numeric(input$curve_cd1), 
                          rd = as.numeric(input$curve_rd1),
                          sigma = as.numeric(input$curve_sigma1))
    }
    #
    # Define parameters for model 2
    #
    if(current.model2 == "theis" ){
      current.par2 <- NULL
    }
    if(current.model2 == "cooper_jacob"){
      current.par2 <- NULL
    }
    if(current.model2 == "boulton"){
      current.par2 <- list(coeffs = coeffs,
                           sigma = as.numeric(input$curve_sigma2),
                           phi = as.numeric(input$curve_phi2))
    }
    if(current.model2 == "hantush_jacob"){
      current.par2 <- list(coeffs = coeffs,
                           beta = as.numeric(input$curve_beta2))
    }
    if(current.model2 == "papadopulos_cooper"){
      current.par2 <- list(coeffs = coeffs,
                           cd = as.numeric(input$curve_cd2),
                           rho = as.numeric(input$curve_rho2))
    }
    if(current.model2 == "general_radial_flow"){
      current.par2 <- list(coeffs = coeffs,
                          n = as.numeric(input$curve_n2),
                          rd = as.numeric(input$curve_rd2))
    }
    if(current.model2 == "agarwal_skin"){
      #print(coeffs)
      if(is.null(input$curve_cd)) return(NULL)
      current.par2 <- list(coeffs = coeffs, 
                           cd = as.numeric(input$curve_cd2), 
                           rd = as.numeric(input$curve_rd2),
                           sigma = as.numeric(input$curve_sigma2))
    }
    #
    # Calculate well functions
    #
    W1 <- NULL
    W2 <- NULL
    if(current.model1 != "cooper_jacob" & current.model1 != "no_flow_boundary" & 
       current.model1 != "constant_head_boundary"){
      W1 <- do.call(current.fn1, args = list(td = t_vector, par = current.par1))
    }
    else if(current.model1 == "cooper_jacob"){
      W1 <- do.call(current.fn1, args = list(u = 1/t_vector))
    }
    else if(current.model1 == "no_flow_boundary"){
      rd <- 1.0
      Rd <- 10.0
      td <- t_vector
      term1 <- pracma::expint_E1(rd^2/(4*td))
      term2 <- pracma::expint_E1((rd^2)*(Rd^2)/(4*td))
      W1 <- 0.5*(term1+term2)#theis_well_function(td1) + theis_well_function(td2)
    }
    else if(current.model1 == "constant_head_boundary"){
      rd <- 1.0
      Rd <- 20.0
      td <- t_vector
      term1 <- pracma::expint_E1(rd^2/(4*td))
      term2 <- pracma::expint_E1((rd^2)*(Rd^2)/(4*td))
      W1 <- 0.5*(term1-term2)
    }
    #
    if(current.model2 != "cooper_jacob" & current.model2 != "no_flow_boundary" & 
       current.model2 != "constant_head_boundary"){
      W2 <- do.call(current.fn2, args = list(td = t_vector, par = current.par2))
    }
    else if(current.model2 == "cooper_jacob"){
      W2 <- do.call(current.fn2, args = list(u = 1/t_vector))
    }
    else if(current.model2 == "no_flow_boundary"){
      rd <- 1.0
      Rd <- 10.0
      td <- t_vector
      term1 <- pracma::expint_E1(rd^2/(4*td))
      term2 <- pracma::expint_E1((rd^2)*(Rd^2)/(4*td))
      W2 <- 0.5*(term1+term2)#theis_well_function(td1) + theis_well_function(td2)
    }
    else if(current.model2 == "constant_head_boundary"){
      rd <- 1.0
      Rd <- 20.0
      td <- t_vector
      term1 <- pracma::expint_E1(rd^2/(4*td))
      term2 <- pracma::expint_E1((rd^2)*(Rd^2)/(4*td))
      W2 <- 0.5*(term1-term2)
    }
    #
    ctitle <- paste0("Well Function - ", current.model1,' ',current.model2)
    W.df <- data.frame(t = t_vector, W1 = W1, W2 = W2)
    #print(W.df)
    p.comp <- ggplot() + geom_line(aes(x = t, y = W1, color = colors[1]), data = W.df) + 
      geom_line(aes(x = t, y = W2, color = colors[2]), data = W.df) + 
      ggtitle(paste0("Comparison: ",ctitle)) + 
      scale_color_manual(labels = c(current.model1, current.model2), values = colors) +
      theme_bw() + 
      theme(axis.text.x = element_text(size = 14), 
            axis.title.x = element_text(size = 16),
            axis.text.y = element_text(size = 14), 
            axis.title.y = element_text(size = 16),
            legend.title = element_text(size = 14),
            legend.text = element_text(size = 14), 
            plot.title = element_text(size = 18, face = "bold")) +
      scale_x_log10() +
      guides(color=guide_legend("Variables"))
    if(input$compare_scale1 == "LogLog")
      p.comp <- p.comp + scale_y_log10()
    return(p.comp)
  })
  #
})
