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
library(DT)
library(rhandsontable)
library(shinyalert)
########################################################################################
##                           Global variables
########################################################################################
# Max file input size :
options(shiny.maxRequestSize=30*1024^2)
#
model.types <- list( "Pumping" =
                       c("None", "theis", "cooper_jacob", "hantush_jacob",
                         "boulton", "general_radial_flow",
                         "papadopulos_cooper", "agarwal_skin",
                         "agarwal_recovery", "warren_root",
                         "gringarten"),
                     "Slug" = c("None", "cooper", "neuzil"),
                     "Bailer" = c("None", "mcdonald"))
#
smooth.types <- c("None", "spline", "kersmooth", "lokern", "lpridge")
#
derivative.types <-  c("None", "central", "horner", "bourdet", "spline", "spane",
                       "smoothspline", "kernelreg", "locpol", "lokern",
                       "lpridge", "wavelet")
#
model.par <- list("theis" = c("a", "t0"),
                  "cooper_jacob" = c("a", "t0"),
                  "hantush_jacob" = c("a", "t0", "lambda"),
                  "boulton" = c("a", "t0", "t1", "phi"),
                  "agarwal_recovery" = c("a", "t0"),
                  "general_radial_flow" = c("a", "t0", "n"),
                  "papadopulos_cooper" = c("a", "t0"),
                  "warren_root" = c("a", "t0", "t1", "tm"),
                  "gringarten" = c("a", "t0")
                  )
#
base.par <- c("Tr", "Ss")
hydr.par <- list("theis" = base.par,
                 "cooper_jacob" = c(base.par),
                 "hantush_jacob" = c(base.par, "Ka"),
                 "boulton" = c(base.par, "omegad"),
                 "general_radial_flow" = c(base.par),
                 "papadopulos_cooper" = c(base.par, "cd"),
                 "agarwal_skin" = c(base.par, "cd"),
                 "agarwal_recovery" = c(base.par),
                 "warren_root" = c(base.par),
                 "gringarten" = c(base.par))
#
hydr.par1 <- list("theis" = c(base.par,"Ri"),
                  "cooper_jacob" = c(base.par, "Ri"),
                  "hantush_jacob" = c(base.par, "Ka"),
                  "boulton" = c(base.par, "omegad"),
                  "general_radial_flow" = c(base.par),
                  "papadopulos_cooper" = c(base.par, "cd"),
                  "agarwal_skin" = c(base.par, "cd"),
                  "agarwal_recovery" = c(base.par, "Ri"),
                  "warren_root" = c(base.par),
                  "gringarten" = c(base.par))
#
hydr.par.mcmc <- list("theis" = c(base.par, "sigma"),
                 "cooper_jacob" = c(base.par, "sigma"),
                 "hantush_jacob" = c(base.par, "Ka", "sigma"),
                 "boulton" = c(base.par, "omegad", "sigma"),
                 "general_radial_flow" = c(base.par, "sigma"),
                 "papadopulos_cooper" = c(base.par, "cd", "sigma"),
                 "agarwal_skin" = c(base.par, "cd", "sigma"),
                 "agarwal_recovery" = c(base.par, "sigma"),
                 "warren_root" = c(base.par, "sigma"),
                 "gringarten" = c(base.par, "sigma"))
#
obj.fn <- list("None",
               "Residual Sum of Squares",
               "Mean Absolute Deviation",
               "Max. Absolute Deviation",
               "LogLikelihood")
#
obj.fn1 <- list("None" = NULL,
                "Residual Sum of Squares" = "residual_sum_squares",
                "Mean Absolute Deviation" = "mean_absolute_deviation",
                "Max. Absolute Deviation" = "max_absolute_deviation",
                "LogLikelihood" = "loglikelihood")
#
opt.methods <- list("None" = NULL,
                   "nls" = "nls", "l-bfgs-b" = "l-bfgs-b", "SimulatedAnnealing" = "sa",
                   "GeneticAlgorithms" = "ga", "ParticleSwarmOptimization" = "pso",
                   "DifferentialEvolution" = "de", "CopulasEdas" = "copulaedas")
#
pdf.models <- list("None" = NULL, "Uniform" = "unif", "Normal" = "norm")
# Define server
shinyServer(function(input, output, session) {
  ## Server variables
  server.env <- environment() # used to allocate in functions
  current.table <- NULL
  current.ptest <- NULL # this variable will contain the current pumping test
  initial.solution.flag <- FALSE
  initial.solution <- NULL
  first <- TRUE
  current.hydr.par <- NULL
  cint_jack <- NULL
  confint_wald <- NULL
  ## Panel 'About' (left hand side)
  # Output the uptc logo :
  output$uptc.logo <- renderImage(list(src = "uptc_jpg.jpg"),
                                      deleteFile = FALSE)
  ## Panel 'Import data'
  dInput <- reactive({
    in.file <- input$file1
    #
    validate(
      need(input$file1, 'Check if file is loaded')
    )
    #
    if (is.null(in.file))
      return(NULL)
    #
    fname <- strsplit(input$file1$name, "\\.")
    updateTextInput(session, inputId = "testname", value = fname[[1]][1])
    #
    the.sep <- switch(input$sep, "Comma"=",", "Semicolon"=";", "Tab"="\t",
                      "Space"="")
    the.quote <- switch(input$quote, "None"="","Double Quote"='"',
                        "Single Quote"="'")
    the.dec <- switch(input$dec, "Period"=".", "Comma"=",")
    if (input$rownames) {
      the.table <- read.table(in.file$datapath, header=input$header,
                              sep=the.sep, quote=the.quote, row.names=1,
                              dec=the.dec)
    } else {
      the.table <- read.table(in.file$datapath, header=input$header,
                              sep=the.sep, quote=the.quote, dec=the.dec)
    }
    if(!first){
      server.env$current.ptest <- NULL
      server.env$initial.solution.flag <- FALSE
      # Reset some values in the test info tab
      updateTextInput(session, inputId = "testname", value = "Test0")
      updateNumericInput(session, inputId = "pumprate", value = "0.0")
      updateNumericInput(session, inputId = "distance", value = "0.0")
      updateCheckboxInput(session, inputId = "selectB", value = FALSE)
      updateCheckboxInput(session, inputId = "selectrw", value = FALSE)
      updateCheckboxInput(session, inputId = "selectrc", value = FALSE)
      updateCheckboxInput(session, inputId = "selectCeff", value = FALSE)
      updateCheckboxInput(session, inputId = "selectV", value = FALSE)
      updateTextInput(session, inputId = "filter.min", value = min(the.table$t))
      updateTextInput(session, inputId = "filter.max", value = max(the.table$t))
      output$plot.simple <- renderPlot({
        t <- the.table$t
        s <- the.table$s
        plot(t, s, plot ="p", log = "x", main = "Current Pumping Test")
      })

      # Reset values in the diagnostic plot tab
      updateSelectInput(session, inputId = "derivative", selected = "None")
      # Reset values in the parameter estimation tab
      updateSelectInput(session, inputId = "parest_method", selected = "None")
      updateSelectInput(session, inputId = "parest_objfn1", selected = "None")
      updateSelectInput(session, inputId = "parest_model1", selected = "None")
      updateSelectInput(session, inputId = "parest_derivative1", selected = "None")
      # Reset model diagnostic tab
      updateSelectInput(session, inputId = "diagnostic.type", selected = "None")
      server.env$cint_jack <- NULL
      # Reset values in the UQ tab
      updateSelectInput(session, input = "uncertainty.method", selected = "None")
      updateTextInput(session, input = "uncert_nreal1", value= "100")
      # Reset values in the Drawdown Calculator tab
      updateTextInput(session, inputId = "calc_pumprate", value = "0.0")
      updateTextInput(session, inputId = "calc_distance", value = "0.0")
      updateTextInput(session, inputId = "calc_time", value = "0.0")
      updateTextInput(session, "calc_transmissivity", value = "1e-16")
      updateTextInput(session, "calc_storage", value = "1.0e-16")
      updateSelectInput(session, "calc_model", selected = "None")
    }
    if(first)
      first <- FALSE
    server.env$first <- first
    # return the table
    server.env$current.table <- the.table
    the.table
  })
  # data preview table
  output$view <- renderTable({
    d.input <- dInput()
    if (is.null(d.input))
      return(NULL)
    if (ncol(d.input)>input$ncol.preview)
      d.input <- d.input[,1:input$ncol.preview]
    head(d.input, n=input$nrow.preview)
  })
  #
  output$summary <- renderPrint({
    d.input <- dInput()
    if (is.null(d.input))
      return(NULL)
   the.table <- server.env$current.table
   if(is.null(the.table))
     return(NULL)
    summary(the.table)
  })
  # update prate
  output$prate <- renderUI({
    res <- NULL
    if(input$type == "Pumping"){
      res <- textInput(inputId = "pumprate", label = "Pumping rate (m3/s)",
                value = "0.0")
    }
    if(input$type == "VariableRate"){
      res <- HTML("<br><h4>Variable Rate Data</h4><br>")
    }
    return(res)
  })
  #
  output$vrate <- renderRHandsontable({
    res <- NULL
    if(input$type == "VariableRate"){
      DF <- data.frame(t = seq(60.0, 300.0, by = 60.0), Q = rep(0.0, 5))
      res <- rhandsontable(DF, useTypes = TRUE, stretchH = "all",
                           digits = 2,
                           colHeaders = c("Time[s]", "Discharge[m3/s]"),
                           width = 250)
    }
    return(res)
  })

  # update dist
  output$dist <- renderUI({
    res <- NULL
    if(input$type == "Pumping"){
      res <- textInput(inputId = "distance", label = "Distance(m)", value = "0.0")
    }
    else if(input$type == "VariableRate"){
      res <- textInput(inputId = "distance", label = "Distance(m)", value = "0.0")
    }
    return(res)
  })
  # update addpar1
  output$addpar1 <- renderUI({
    checkboxInput(inputId = "selectB", label = "Aquitard Thickness",
                           value = FALSE)
  })
  #
  output$setpar1 <- renderUI({
    tmp <- NULL
    if(is.null(input$selectB)) return(NULL)
    if(input$selectB){
      tmp <- textInput(inputId = "setB", "Aquitard Thickness(m)= ", value = "0.0")
    }
    return(tmp)
  })
  # update addpar2
  output$addpar2 <- renderUI({
    checkboxInput(inputId = "selectrw", label = "Radius Well", value = FALSE)
  })
  #
  output$setpar2 <- renderUI({
    tmp <- NULL
    if(is.null(input$selectrw)) return(NULL)
    if(input$selectrw){
      tmp <- textInput(inputId = "setrw", label = "Radius Well(m)= ", value = "0.0")
    }
    return(tmp)
  })
  # update addpar3
  output$addpar3 <- renderUI({
    checkboxInput(inputId = "selectrc", label = "Radius Casing", value = FALSE)
  })
  #
  output$setpar3 <- renderUI({
    tmp <- NULL
    if(is.null(input$selectrc)) return(NULL)
    if(input$selectrc){
      tmp <- textInput(inputId = "setrc", label = "Radius Casing(m)= ", value = "0.0")
    }
    return(tmp)
  })
  # update addpar4
  output$addpar4 <- renderUI({
    checkboxInput(inputId = "selectCeff", label = "Effective Compressibility", value = FALSE)
  })
  #
  output$setpar4 <- renderUI({
    tmp <- NULL
    if(is.null(input$selectCeff)) return(NULL)
    if(input$selectCeff){
      tmp <- textInput(inputId = "setCeff", label = "Effective Compressibility(Pa^-1)= ",
                       value = "0.0")
    }
    return(tmp)
  })
  # update addpar5
  output$addpar5 <- renderUI({
    checkboxInput(inputId = "selectV", label = "Volume Test Section", value = FALSE)
  })
  #
  output$setpar5 <- renderUI({
    tmp <- NULL
    if(is.null(input$selectV)) return(NULL)
    if(input$selectV){
      tmp <- textInput(inputId = "setV", label = "Volume Test Section(m3)= ", value = "0.0")
    }
    return(tmp)
  })
  #
  output$filter.par1 <- renderUI({
    the.table <- server.env$current.table
    if(is.null(the.table)){
      return(NULL)
    }
    tmp <- textInput(inputId = "filter.min", label = "Minimum", value = min(the.table$t))
    return(tmp)
  })
  #
  output$filter.par2 <- renderUI({
    the.table <- server.env$current.table
    if(is.null(the.table)){
      return(NULL)
    }
    tmp <- textInput(inputId = "filter.max", label = "Maximum", value = max(the.table$t))
    return(tmp)
  })

  # Create pumpingtest object
  create_ptest <- function(){
    the.table <- server.env$current.table
    input$defineptest
    if(input$type == "Pumping"){
      server.env$current.ptest <- isolate(
        pumping_test(id = input$testname,
                     Q = as.numeric(input$pumprate),
                     r = as.numeric(input$distance),
                     t = as.numeric(the.table$t),
                     s = as.numeric(the.table$s))
      )
      #
      add.par <- list()
      if(input$selectB){
        add.par$B <- as.numeric(input$setB)
      }
      #
      if(input$selectrw){
        add.par$rw <- as.numeric(input$setrw)
      }
      #
      if(input$selectrc){
        add.par$rc <- as.numeric(input$setrc)
      }
      #
      if(input$selectCeff){
        add.par$Ceff <- as.numeric(input$setCeff)
      }
      #
      if(input$selectV){
        add.par$Vs <- as.numeric(input$setV)
      }
      server.env$current.ptest$additional_parameters <- add.par
    }
    else if(input$type == "Slug"){
      server.env$current.ptest <- isolate(
        pumping_test(id = input$testname, Q = 0.0, r = 0.0,
                     t = as.numeric(the.table$t),
                     s = as.numeric(the.table$s))
      )
    }
    print("pumpingtest defined")
    server.env$current.ptest
  }
  #
  output$plot.simple <- renderPlot({
    the.table <- server.env$current.table
    if(is.null(the.table)) return(NULL)
    t <- the.table$t
    s <- the.table$s
    plot(t, s, type = "p", log = "x", main = "Current Pumping Test")
  })
  #
  observeEvent(input$defineptest, {
    create_ptest()
    shinyalert(title = "Pumpingtest object Defined!!!", type = "success")
  })
  #
  # Smoothing tab
  #
  output$smoothmethod0 <- renderUI({

  })

  #
  # Diagnostic Plot tab
  #
  output$diagnostic.table <- renderImage({
    tmp <- NULL
    current.ptest <- server.env$current.ptest
    if(is.null(current.ptest)) return(list(src=""))
    #if(is.null(input$scale)) return(NULL)
    if(input$type == "Pumping"){
      if(input$scale == "loglog"){
        tmp <- list(src = "diagnostic_plots_log.png", width = 600,
                    height = 450)
      }
      else if(input$scale == "slog"){
        tmp <- list(src = "diagnostic_plots_semilog.png", width = 600,
                    height = 450)
      }
    }
    else if(input$type == "Slug"){
      tmp <- list(src = "", width = 600,
                  height = 450)
    }
    return(tmp)
  },
    deleteFile = FALSE)
  #
  output$derivative0 <- renderUI({
    selectInput(inputId = "derivative", label = "Derivative Type",
                choices = switch(input$type,
                                 "Slug" = derivative.types,
                                 "Pumping" = derivative.types),
                selected = switch(input$type, "Slug"="None",
                                  "Pumping"="None"))
  })
  #
  output$diagnostic.plot <- renderPlot({
    current.table <- server.env$current.table
    current.ptest <- server.env$current.ptest
    #validate(
    #  need(current.ptest & current.table, "The Pumping Test has not been defined.")
    #)
    if(is.null(current.table))
      return(NULL)
    if(is.null(current.ptest))
      return(NULL)
    if(is.null(input$derivative) ||  input$derivative == "None")
      return(NULL)
    #
    tmp.plt <- NULL
    #
    if(input$type == "Pumping"){
      dd <- 2
      if(input$derivative == "bourdet"){
        dd <- as.numeric(input$points)
      }
      if(input$derivative == "spline"){
        dd <- as.numeric(input$sppoints)
      }
      if(input$derivative == "spane"){
        dd <- as.numeric(input$sppoints1)
      }
      if(input$derivative == "kernelreg"){
        dd <- as.numeric(input$krbw)
      }
      tmp.plt <- plot(current.ptest, type = "diagnostic",
                      dmethod = input$derivative, d = dd,
                      scale = input$scale,
                      slug = FALSE,
                      results = FALSE, y.intersp = 1)
    }
    else if(input$type == "Slug"){
      dd <- 2
      if(input$derivative == "bourdet"){
        dd <- as.numeric(input$points)
      }
      if(input$derivative == "spline"){
        dd <- as.numeric(input$sppoints)
      }
      if(input$derivative == "spane"){
        dd <- as.numeric(input$sppoints1)
      }
      if(input$derivative == "kernelreg"){
        dd <- as.numeric(input$krbw)
      }
      tmp.plt <- plot(current.ptest, type = "diagnostic",
                      dmethod = input$derivative, d = dd,
                      scale = input$scale,
                      slug = TRUE,
                      results = FALSE, y.intersp = 1)
    }
    tmp.plt
  })
  #
  # output$filter.slider <- renderUI({
  #   current.ptest <- server.env$current.ptest
  #   if(is.null(current.ptest))
  #     return(NULL)
  #   if(!input$filter.check)
  #     return(NULL)
  #   mnt <- min(current.ptest$t)
  #   mxt <- max(current.ptest$t)
  #   step <- (mxt-mnt)/100
  #   tmp <- sliderInput(inputId = "filter.slider1", label = "Filter Range", min = mnt,
  #                      max = mxt, step = step, value = c(1.5*mnt,0.8*mxt))
  # })
  ######################################################################################
  #                        Parameter estimation tab
  ######################################################################################
  output$parest_model <- renderUI({
    selectInput(inputId = "parest_model1", label = "Solution Type",
                choices = switch(input$type,
                               "Slug" = model.types[["Slug"]],
                               "Pumping" = model.types[["Pumping"]],
                               "Bailer" = mode.types[["Bailer"]]),
               selected = switch(input$type, "Slug"="None",
                                           "Pumping"="None",
                                           "Bailer"="None"))
  })
  #
  #
  output$parest_opt <- renderUI({
    if(input$parest_method == "automatic"){
      selectInput(inputId = "parest_opt1", label = "Optimization Method",
                  choices = c("nls", "l-bfgs-b", "SimulatedAnnealing",
                              "GeneticAlgorithms",
                              "ParticleSwarmOptimization",
                              "DifferentialEvolution"),
                  selected = "nls")
    }
  })
  #
  output$parest_objfn <- renderUI({
    if(input$parest_method == "automatic"){
      selectInput(inputId = "parest_objfn1", label = "Objetive Function",
                  choices = c("None", "Residual Sum of Squares", "Mean Absolute Deviation",
                              "Max. Absolute Deviation", "LogLikelihood"))
    }
  })
  #
  output$parest_opt_par_flag <- renderUI({
    if(input$parest_method == "automatic"){
      checkboxInput(inputId = "opt_par_flag1", label = "Optimization Parameters",
                    value = FALSE)
    }
  })
  #
  output$parest_logTr <- renderUI({
    if(input$parest_method == "visual"){
      sliderInput("LogTr","Log10(Transmissivity(m2/s))", min=-16.0,
                  max = -15.0, value = -15.5, step = 0.1)
    }
  })
  #
  output$parest_logSs <- renderUI({
    if(input$parest_method == "visual"){
      sliderInput("LogSs", "Log10(Storage Coefficient)", min=-16.0,
                  max = -15.0, value = -15.5, step = 0.1)
    }
  })
  #
 output$parest_omegad <- renderUI({
   if( input$parest_method =="visual" &&
       input$parest_model1 == "boulton"){
       sliderInput("omegad", "Log10(Drain. Porosity)", min = -16.0,
                   max = -15, value = -15.5, step = 0.1)
   }
 })
 #
 output$parest_phi <- renderUI({
   if(input$parest_method == "visual" &&
      input$parest_model1 == "boulton"){
       sliderInput("phi", "Log10(Phi)", min = -16.0,
                   max =-15, value = -15.5, step = 0.1)
     }
 })
#
 output$parest_hj_ka <- renderUI({
   if(input$parest_method == "visual" &&
      input$parest_model1 == "hantush_jacob"){
       sliderInput("Ka", "Log10(Ka)", min = -16.0,
                   max = -15.0, value = -15.5, step = 0.1)
   }
 })
 #
 output$parest_pc_cd <- renderUI({
   if(input$parest_method == "visual" &&
      input$parest_model1 == "papadopulos_cooper"){
     sliderInput("cd", "Log10(cd)", min = -16.0, max = -15.0, value = -15.5, step = 0.1)
   }
 })
 #
 output$parest_grf_n <- renderUI({
   if(input$parest_method == "visual" &&
      input$parest_model1 == "general_radial_flow"){
      sliderInput("n", "Flow Dimension", min= 1e-3, max = 3, value = 2, step = 0.1)
   }
 })
 #
 output$parest_wr_LogSm <- renderUI({
   if(input$parest_method == "visual" &&
      input$parest_model1 == "warren_root"){
     sliderInput(inputId = "LogSm", label = "Log10(Sm)", min = -16, max = -15,
                 value = -15.5, step = 0.1)
   }
 })
 #
 output$parest_wr_lambda <- renderUI({
   if(input$parest_method == "visual" &&
      input$parest_model1 == "warren_root"){
     sliderInput(inputId = "lambda", label = "Lambda", min = -16, max = -15,
                 value = -15.5, step = 0.1)
   }
 })
 #
 output$parest_derivative <- renderUI({
   selectInput(inputId = "parest_derivative1", label = "Derivative Type",
               choices = switch(input$type,
                                "Slug" = derivative.types,
                                "Pumping" = derivative.types,
                                "Bailer" = derivative.types),
               selected = switch(input$type, "Slug"="None",
                                 "Pumping"="None",
                                 "Bailer"="None"))
 })
 #
 output$parest_dpoints <- renderUI({
   tmp <- NULL
   if(!is.null(input$parest_derivative1)){
     if(input$parest_derivative1 == "bourdet"){
       tmp <- numericInput("parest_dpoints1", "Bourdet: Number of points", 2, min = 1,
                           max = 10, step = 1)
     }
     else if(input$parest_derivative1 == "spline"){
       tmp <- numericInput("parest_dpoints1", "Spline: Number of points", 20, min = 1,
                           max = 40, step = 1)
     }
     else if(input$parest_derivative1 == "spane"){
       tmp <- numericInput("parest_dpoints1", "Spane: Number of points", 2, min = 1,
                           max = 40, step = 1)
     }
   }
   tmp
 })
 #
 output$parest_scale <- renderUI({
   radioButtons("parest_scale1", "Plot Scale", c("loglog","slog"), selected = "loglog")
 })
 #
 get_objective_fn_value <- function(current.ptest){
   obj.fn <- 'residual_sum_squares'
   if(input$parest_method == "automatic"){
     obj.fn <- obj.fn1[[input$parest_objfn1]]
     #print('OBJFN')
     #print(c(input$parest_objfn1, obj.fn))
     if(is.null(obj.fn)) {
       obj.fn <- 'residual_sum_squares'
       updateSelectInput(session, inputId = "parest_objfn1",
                         selected = "Residual Sum of Squares")
     }
   }
   par <- as.numeric(current.ptest$parameters)
   if(obj.fn == "loglikelihood"){
     par <- as.numeric(current.ptest$parameters1)
   }
   args <- list(par = par, ptest = current.ptest, model = input$parest_model1)
   obj.fn_value <- do.call(obj.fn, args = args)
   return(obj.fn_value)
 }
 #
 output$estimation.plot <- renderPlot({
   current.table <- server.env$current.table
   current.ptest <- server.env$current.ptest
   if(is.null(current.table))
     return(NULL)
   if(is.null(current.ptest))
     return(NULL)
   if(input$parest_method == "None")
     return(NULL)
   if(input$parest_model1 == "None")
     return(NULL)
   if(input$parest_derivative1 == "None")
     return(NULL)
   if(input$parest_scale1 == "None")
     return(NULL)
   #
   tmp.plt <- NULL
   #
   if(input$parest_method == "automatic"){
     if(is.null(input$parest_opt1))
       return(NULL)
     if(input$parest_opt1 == "None")
       return(NULL)
     # Estimate hydraulic parameters
     opt.method <- NULL
     opt.method <- opt.methods[[input$parest_opt1]]
     #print(opt.method)
     #
     obj.fn <- NULL
     if(input$parest_objfn1 == "Residual Sum of Squares") obj.fn <- "rss"
     if(input$parest_objfn1 == "Mean Absolute Deviation") obj.fn <- "mnad"
     if(input$parest_objfn1 == "Max. Absolute Deviation") obj.fn <- "mxad"
     if(input$parest_objfn1 == "LogLikelihood") obj.fn <- "loglik"
     if(input$parest_opt1 == "nls" ) obj.fn <- "rss"
     #print(opt.method)
     #print(c(opt.method, input$parest_objfn1, obj.fn))
     current.fit <- fit.optimization(current.ptest, input$parest_model1,
                                     obj.fn = obj.fn,
                                     opt.method = opt.method, seed = 54321)
     hydraulic.parameters(current.ptest) <- current.fit$hydraulic_parameters
     fit.parameters(current.ptest) <- current.fit$parameters
     model(current.ptest) <- input$parest_model1
     estimated(current.ptest) <- TRUE
     server.env$confint_wald <- confint_wald(current.ptest)
     #print(confint_wald(current.ptest))
     if(obj.fn == "loglik"){
       current.ptest$parameters1 <- current.fit$parameters1
     }
     server.env$current.ptest <- current.ptest
   }
   else if(input$parest_method == "visual"){
     if(input$parest_model1 == "None")
       return(NULL)
     # Calculate initial solution if required
     initial.sol <- NULL
     if(!server.env$initial.solution.flag){
       server.env$initial.solution.flag <- TRUE
       fn <- paste0(input$parest_model1,"_solution_initial")
       #print(fn)
       initial.sol <- do.call(fn,list(ptest = server.env$current.ptest))
       print(initial.sol)
       server.env$initial.solution <- initial.sol
       #
       fn1 <- paste0(input$parest_model1, "_calculate_parameters")
       initial.par <- do.call(fn1, list(ptest = server.env$current.ptest,
                                        par = initial.sol,
                                        hydraulic = TRUE))
       #print(initial.par)
       # Update sliders
       Tr <- initial.par$Tr
       Trmn <- floor(log10(Tr))-3
       Trmx <- floor(log10(Tr))+3
       if(input$parest_model1 == "general_radial_flow"){
         Trmx <- Trmx + 5
       }
       Trstp <- (Trmx-Trmn)/100
       updateSliderInput(session, inputId = "LogTr", value = log10(Tr), min = Trmn,
                         max = Trmx, step = Trstp)
       Ss <- initial.par$Ss
       Ssmn <- floor(log10(Ss))-3
       Ssmx <- floor(log10(Ss))+3
       if(input$parest_model1 == "general_radial_flow"){
         Ssmx <- Ssmx + 5
       }
       Ssstp <-(Ssmx - Ssmn)/100
       #print(c(Ssmn, Ssmx, Ssstp))
       updateSliderInput(session, inputId = "LogSs", value = log10(Ss), min = Ssmn,
                         max = Ssmx, step = Ssstp)
       if(input$parest_model1 == "boulton"){
         Omegad <- initial.par$omegad
         Omegadmn <- floor(log10(Omegad))-3
         Omegadmx <- floor(log10(Omegad))+3
         Omegadstp <- (Omegadmx-Omegadmn)/100
         updateSliderInput(session, inputId = "omegad", value = log10(Omegad),
                           min = Omegadmn, max = Omegadmx, step = Omegadstp)
         #
         phi <- -4 #initial.par$phi
         phimn <- -6
         phimx <- -1
         phistp <- (phimx-phimn)/100
         updateSliderInput(session, inputId = "phi", value = phi, min = phimn,
                           max = phimx, step = phistp)
       }
       #
       if(input$parest_model1 == "hantush_jacob"){
         Ka <- initial.par$Ka
         Kamn <- floor(log10(Ka))-3
         Kamx <- floor(log10(Ka))+3
         Kastp <- (Kamx-Kamn)/100
         updateSliderInput(session, inputId = "Ka", value = log10(Ka), min = Kamn,
                           max = Kamx, step = Kastp)
       }
       #
       if(input$parest_model1 == "papadopulos_cooper"){
         cd <- initial.par$cd
         cdmn <- floor(log10(cd))-3
         cdmx <- floor(log10(cd))+3
         cdstp <- (cdmx - cdmn)/100
         updateSliderInput(session, inputId = "cd", value = log10(cd), min = cdmn,
                           max = cdmx, step = cdstp)
       }
       #
     }
     # Read model parameters from ui
     Tr <- 10^(base::as.numeric(input$LogTr))
     Ss <- 10^(base::as.numeric(input$LogSs))
     endp <- length(current.ptest$t)
     Ri <- 2.0*sqrt(Tr*current.ptest$t[endp]/Ss)
     par_rev <- list(Tr = Tr, Ss = Ss, radius_influence = Ri)
     if(input$parest_model1 == 'boulton'){
       omegad <- 10^(base::as.numeric(input$omegad))
       par_rev$omegad <- omegad
       par_rev$phi <- 10^(base::as.numeric(input$phi))
     }
     if(input$parest_model1 == 'hantush_jacob'){
       Ka <- 10^(base::as.numeric(input$Ka))
       par_rev$Ka <- Ka
     }
     if(input$parest_model1 == 'papadopulos_cooper'){
       cd <- 10^(base::as.numeric(input$pc_cd))
       par_rev$cd <- cd
       par_rev$rho <- ptest$r/as.numeric(input$pc_rw)
     }
     if(input$parest_model1 == 'general_radial_flow'){
       #print(par_rev)
       n <- base::as.numeric(input$n)
       par_rev$n <- n
     }
     fit_parameters_fn <- paste0(input$parest_model1,'_calculate_parameters')
     fit.par <- do.call(fit_parameters_fn, list(ptest = current.ptest, par = par_rev,
                                                hydraulic = FALSE))
     current.ptest$hydraulic_parameters <- par_rev
     fit.parameters(current.ptest) <- fit.par
     model(current.ptest) <- input$parest_model1
     estimated(current.ptest) <- TRUE
     server.env$current.ptest <- current.ptest
   }#end visual
   #
   output$parest_results <- renderUI({
     if(input$parest_method == "None")
       return(NULL)
     if(!server.env$current.ptest$estimated)
       return(NULL)
     obj_fn_value <- get_objective_fn_value(server.env$current.ptest)
     str1 <- "<h3>Results Parameter Estimation</h3><br>"
     str2 <- paste("<b>Analytical Model:</b>", input$parest_model1,"<br>")
     str3 <- NULL
     str4 <- NULL
     str5 <- "<br><br><h3>Parameter Confidence Intervals</h3>"
     str6 <- NULL
     if(input$parest_method == "visual"){
       str3 <- paste("<b>Objective Function Value</b>:", obj_fn_value,"<br>")
     }
     #
     if(input$parest_method == "automatic"){
       str3 <- paste("<b>Optimization Method</b>:", input$parest_opt1, "<br>",
                     "<b>Objective Function</b>:", input$parest_objfn1, "<br>",
                     "<b>Objective Function Value </b>:", obj_fn_value, "<br>")
     }
     #
     if(input$parest_model1 == "theis" | input$parest_model1 == "cooper_jacob"){
       Tr <- format(current.ptest$hydraulic_parameters$Tr, digits = 5, scientific = TRUE)
       Ss <- format(current.ptest$hydraulic_parameters$Ss, digits = 5, scientific = TRUE)
       ri <- format(current.ptest$hydraulic_parameters$radius_influence, digits = 3)
       str4 <- paste("<b>Transmisivity (m2/s):</b>", as.character(Tr), "<br>",
                     "<b>Storage Coefficient:</b>", as.character(Ss), "<br>",
                     "<b>Radius Influence (m):</b>", as.character(ri), "<br>")
     }
     #
     if(input$parest_model1 == "boulton"){
       Tr <- format(current.ptest$hydraulic_parameters$Tr, digits = 5, scientific = TRUE)
       Ss <- format(current.ptest$hydraulic_parameters$Ss, digits = 5, scientific = TRUE)
       omegad <- format(current.ptest$hydraulic_parameters$omegad, digits = 5, scientific = TRUE)
       phi <- format(current.ptest$hydraulic_parameters$phi, digits = 5, scientific = TRUE)
       ri <- format(current.ptest$hydraulic_parameters$radius_influence, digits = 3)
       str4 <- paste("<b>Transmisivity (m2/s):</b>", as.character(Tr), "<br>",
                     "<b>Storage Coefficient:</b>", as.character(Ss), "<br>",
                     "<b>Drainage Porosity:</b>", as.character(omegad), "<br>",
                     "<b>Radius influence (m)</b>:", as.character(ri), "<br>")
     }
     #
     if(input$parest_model1 == "hantush_jacob"){
       Tr <- format(current.ptest$hydraulic_parameters$Tr, digits = 5, scientific = TRUE)
       Ss <- format(current.ptest$hydraulic_parameters$Ss, digits = 5, scientific = TRUE)
       Ka <- format(current.ptest$hydraulic_parameters$Ka, digits = 5, scientific = TRUE)
       ri <- format(current.ptest$hydraulic_parameters$radius_influence, digits = 3)
       str4 <- paste("<b>Transmisivity (m2/s):</b>", as.character(Tr), "<br>",
                     "<b>Storage Coefficient:</b>", as.character(Ss), "<br>",
                     "<b>Aquitard Conductivity (m/s):</b>", as.character(Ka), "<br>",
                     "<b>Radius influence (m)</b>:", as.character(ri), "<br>")
     }
     #
     if(input$parest_model1 == "papadopulos_cooper"){
       Tr <- format(current.ptest$hydraulic_parameters$Tr, digits = 5, scientific = TRUE)
       Ss <- format(current.ptest$hydraulic_parameters$Ss, digits = 5, scientific = TRUE)
       cd <- format(current.ptest$hydraulic_parameters$cd, digits = 5, scientific = TRUE)
       str4 <- paste("<b>Transmisivity (m2/s):</b>", as.character(Tr), "<br>",
                     "<b>Storage Coefficient:</b>", as.character(Ss), "<br>",
                     "<b>Wellbore Storage:</b>", as.character(cd), "<br>")
     }
     #
     #print(server.env$confint_wald)
     # ci_wald <- server.env$confint_wald$hydraulic.parameters.ci
     # par_names <- row.names(server.env$confint_wald$hydraulic.parameters.ci)
     # level_names <- colnames(server.env$confint_wald$hydraulic.parameters.ci)
     # print(level_names)
     # str5 <- paste(str5, "<br>", level_names, "<br>")
     # str6 <- ""
     # for(ipar in 1:length(par_names)){
     #   str6 <- paste( str6,
     #                  "<b>", par_names[ipar], ":</b>",
     #                  as.character(ci_wald[ipar,1]), " ",
     #                  as.character(ci_wald[ipar,2]), "<br>")
     # }
     # #print(str6)
      HTML(paste(str1, str2, str3, str4, str5)) #, str6))
   })
   #
   #
   output$ci_parest_results <- renderDataTable({
     if(input$parest_method == "None")
       return(NULL)
     if(!server.env$current.ptest$estimated)
       return(NULL)
     ci_wald <- server.env$confint_wald$hydraulic.parameters.ci
     return(ci_wald)
   },
   extensions = c('Buttons'),
   options = list(
     pageLength = nrow(server.env$confint_wald$hydraulic.parameters.ci),
     dom = 'Bfrtip',
     buttons = c('copy', 'csv', 'excel', 'pdf', 'print'),
     text = 'Download',
     scrollY = 200,
     scroller = TRUE))

   # Now create the estimation plot
   dd <- 2
   if(!is.null(input$parest_dpoints1)){
     dd <- as.numeric(input$parest_dpoints1)
   }
   tmp.plt <- plot(current.ptest, type = "estimation", results = FALSE,
                   dmethod = input$parest_derivative1, d = dd,
                   scale = input$parest_scale1, y.intersp = 1)
   if(input$parest_method == "automatic"){
     server.env$initial.solution.flag <- FALSE
   }
   server.env$cint_jack <- NULL
   tmp.plt
 })#end estimation.plot
 #
 # Output the pumpingtest to be downloaded
 output$pumpingtest.download <- {
   fn <- paste0("pumpingtest",format(Sys.time(),format="-%Y-%m-%d_%H:%M"),".rda",sep="")
   downloadHandler(filename=function() {
     paste0("pumpingtest",format(Sys.time(),format="-%Y-%m-%d_%H:%M"),".rda",sep="")
   },
   content=function(file) {
     ptest.export <- server.env$current.ptest
     save(ptest.export, file=file)
   })
 }

 #######################################################################################
 #                    Diagnostic Model Tab
 #######################################################################################
 output$model.diagnostic <- renderPlot({
   current.ptest <- server.env$current.ptest
   validate(
     need(current.ptest, "The Pumping Test has not been defined.")
   )
   current.type <- input$diagnostic.type
   tmp.plt <- NULL
   if(is.null(current.ptest))
     return(NULL)
   if(!current.ptest$estimated)
     return(NULL)
   #
   if(current.type == "None")
     return(NULL)
   #
   if(current.type == "Model Diagnostic"){
     tmp.plt <- plot(current.ptest, type = "model.diagnostic")
   }
   else if(current.type == "Sample Influence"){
     cint_jack <- NULL
     if(is.null(server.env$cint_jack)){
       cint_jack <- confint_jackniffe(current.ptest)
       server.env$cint_jack <- cint_jack
     }
     else{
       cint_jack <- server.env$cint_jack
     }
     tmp.plt <- plot_sample_influence(cint_jack)
   }
   tmp.plt
 })#end model.diagnostic.plot
 #
 ##############################################################################
 #               Uncertainty Quantification Tab
 ##############################################################################
 # Bootstrapping options: level, d, nreal, seed
 output$uncert_level <- renderUI({
   tmp <- NULL
   if(input$uncertainty.method == "bootstrap"){
     tmp <- textInput("boot_level","Confidence Level",value = "0.025,0.5,0.975")
   }
   tmp
 })
 #
 output$uncert_d <- renderUI({
   if(is.null(server.env$current.ptest))
     return(NULL)
   if(!server.env$current.ptest$estimated)
     return(NULL)
   tmp <- NULL
   if(input$uncertainty.method == "bootstrap"){
     tmp <- textInput("boot_d","Number Points Derivative", value = "20")
   }
   tmp
 })
 #
 output$uncert_nreal <- renderUI({
   tmp <- NULL
   if(is.null(server.env$current.ptest))
     return(NULL)
   if(!server.env$current.ptest$estimated)
     return(NULL)
   tmp <- textInput("uncert_nreal1","Number Realizations",value = "100")
   return(tmp)
 })
 #
 output$uncert_burnin <- renderUI({
   tmp <- NULL
   if(is.null(server.env$current.ptest)){
     return(NULL)
   }
   if(!server.env$current.ptest$estimated)
     return(NULL)
   if(input$uncertainty.method == "bootstrap")
     return(NULL)
   if(!is.null(input$uncertainty.method) | input$uncertainty.method == "amcmc"){
     tmp <- textInput("uncert_burnin1","BurnIn",value = "500")
   }
   return(tmp)
 })
 #
 output$uncert_seed <- renderUI({
   if(is.null(server.env$current.ptest))
     return(NULL)
   if(!server.env$current.ptest$estimated)
     return(NULL)
   tmp <- textInput("uncert_seed1","Random Seed",value = "12345")
 })
 #
 output$uncert_cov.corr <- renderUI({
   tmp <- NULL
   if(is.null(server.env$current.ptest))
     return(NULL)
   if(!server.env$current.ptest$estimated)
     return(NULL)
   if(is.null(input$uncertainty.method))
     return(NULL)
   if(input$uncertainty.method == "amcmc"){
     tmp <- checkboxInput(inputId = "uncert_cov.corr1", label = "Update Covariance matrix", FALSE)
     }
   return(tmp)
 })
 #
 output$uncert_iter_update_cov <- renderUI({
   tmp <- NULL
   if(is.null(server.env$current.ptest))
     return(NULL)
   if(!server.env$current.ptest$estimated)
     return(NULL)
   if(is.null(input$uncertainty.method))
     return(NULL)
   if(input$uncertainty.method == "amcmc"){
     tmp <- textInput(inputId = "uncert_iter_update_cov1", label = "Covariance matrix update Period",
                      value = "20")
   }
   return(tmp)
 })
 #
output$uncert_prior.pdf1 <- renderUI({
  current.ptest <- server.env$current.ptest
  if(is.null(current.ptest))
    return(NULL)
  if(!current.ptest$estimated)
    return(NULL)
  if(input$uncertainty.method == "amcmc" | input$uncertainty.method == "twalk"){
    current.model <- current.ptest$model
    current.par <- hydr.par.mcmc[[current.model]]
    current.label <- paste0(current.par[1], " PDF type")
    selectInput("prior.pdf11", label = current.label,
                c("None", "Uniform", "Normal"),
                selected = "None")
  }
})
#
output$uncert_par1 <- renderUI({
  if(is.null(server.env$current.ptest))
    return(NULL)
  if(!server.env$current.ptest$estimated)
    return(NULL)
  current.par <- hydr.par.mcmc[[server.env$current.ptest$model]]
  uncert_method <- input$uncertainty.method
  tmp <- NULL
  if(uncert_method == "amcmc" | uncert_method == "twalk"){
    current.label = NULL
    if(is.null(input$prior.pdf11)){
      current.label <- paste0(current.par[1], " PDF Parameters")
    }
    else if(input$prior.pdf11 == "None"){
      current.label <- paste0(current.par[1], " PDF Parameters")
    }
    else if(input$prior.pdf11 == "Uniform"){
      current.label <- paste0(current.par[1], " PDF Parameters: min, max")
    }
    else if(input$prior.pdf11 == "Normal"){
      current.label <- paste0(current.par[1], " PDF Parameters: mean, sigma")
    }
    tmp <- textInput("uncert_par11", label = current.label, value = "")
  }
  return(tmp)
})
#
output$uncert_prior.pdf2 <- renderUI({
  current.ptest <- server.env$current.ptest
  if(is.null(current.ptest))
    return(NULL)
  if(!current.ptest$estimated)
    return(NULL)
  if(input$uncertainty.method == "amcmc" | input$uncertainty.method == "twalk"){
    current.model <- current.ptest$model
    current.par <- hydr.par.mcmc[[current.model]]
    current.label <- paste0(current.par[2], " PDF")
    selectInput("prior.pdf22", label = current.label,
                c("None", "Uniform", "Normal"),
                selected = "None")
  }
})
#
#
output$uncert_par2 <- renderUI({
  if(is.null(server.env$current.ptest))
    return(NULL)
  if(!server.env$current.ptest$estimated)
    return(NULL)
  current.par <- hydr.par.mcmc[[server.env$current.ptest$model]]
  uncert_method <- input$uncertainty.method
  tmp <- NULL
  if(uncert_method == "amcmc" | uncert_method == "twalk"){
    current.label = NULL
    if(is.null(input$prior.pdf22)){
      current.label <- paste0(current.par[2], " PDF Parameters")
    }
    else if(input$prior.pdf22 == "None"){
      current.label <- paste0(current.par[2], " PDF Parameters")
    }
    else if(input$prior.pdf22 == "Uniform"){
      current.label <- paste0(current.par[2], " PDF Parameters: min, max")
    }
    else if(input$prior.pdf22 == "Normal"){
      current.label <- paste0(current.par[2], " PDF Parameters: mean, sigma")
    }
    tmp <- textInput("uncert_par22", label = current.label, value = "")
  }
  return(tmp)
})
#
output$uncert_prior.pdf3 <- renderUI({
  tmp <- NULL
  current.ptest <- server.env$current.ptest
  if(is.null(current.ptest))
    return(NULL)
  if(!current.ptest$estimated)
    return(NULL)
  if(input$uncertainty.method == "amcmc" | input$uncertainty.method == "twalk"){
    current.model <- current.ptest$model
    current.par <- hydr.par.mcmc[[current.model]]
    npar <- length(current.par)
    if(npar > 2){
      current.label <- paste0(current.par[3], " PDF")
      tmp <- selectInput("prior.pdf33", label = current.label,
                         c("None", "Uniform", "Normal"),
                         selected = "None")
    }
  }
  return(tmp)
})
#
#
output$uncert_par3 <- renderUI({
  if(is.null(server.env$current.ptest))
    return(NULL)
  if(!server.env$current.ptest$estimated)
    return(NULL)
  current.par <- hydr.par.mcmc[[server.env$current.ptest$model]]
  npar <- length(current.par)
  if(npar <= 2){
    return(NULL)
  }
  uncert_method <- input$uncertainty.method
  tmp <- NULL
  if(uncert_method == "amcmc" | uncert_method == "twalk"){
    current.label = NULL
    if(is.null(input$prior.pdf33)){
      current.label <- paste0(current.par[3], " PDF Parameters")
    }
    else if(input$prior.pdf33 == "None"){
      current.label <- paste0(current.par[3], " PDF Parameters")
    }
    else if(input$prior.pdf33 == "Uniform"){
      current.label <- paste0(current.par[3], " PDF Parameters: min, max")
    }
    else if(input$prior.pdf33 == "Normal"){
      current.label <- paste0(current.par[3], " PDF Parameters: mean, sigma")
    }
    tmp <- textInput("uncert_par33", label = current.label, value = "")
  }
  return(tmp)
})
#
output$uncert_prior.pdf4 <- renderUI({
  current.ptest <- server.env$current.ptest
  if(is.null(current.ptest))
    return(NULL)
  if(!current.ptest$estimated)
    return(NULL)
  if(input$uncertainty.method == "amcmc" | input$uncertainty.method == "twalk"){
    current.model <- current.ptest$model
    current.par <- hydr.par.mcmc[[current.model]]
    npar <- length(current.par)
    if(npar == 4){
      current.label <- paste0(current.par[4], " PDF")
      selectInput("prior.pdf44", label = current.label,
                  c("None", "Uniform", "Normal"),
                  selected = "None")
    }
  }
})
#
output$uncert_par4 <- renderUI({
  if(is.null(server.env$current.ptest))
    return(NULL)
  if(!server.env$current.ptest$estimated)
    return(NULL)
  current.par <- hydr.par.mcmc[[server.env$current.ptest$model]]
  npar <- length(current.par)
  if(npar <= 3){
    return(NULL)
  }
  uncert_method <- input$uncertainty.method
  tmp <- NULL
  if(uncert_method == "amcmc" | uncert_method == "twalk"){
    current.label = NULL
    if(is.null(input$prior.pdf44)){
      current.label <- paste0(current.par[4], " PDF Parameters")
    }
    else if(input$prior.pdf44 == "None"){
      current.label <- paste0(current.par[4], " PDF Parameters")
    }
    else if(input$prior.pdf44 == "Uniform"){
      current.label <- paste0(current.par[4], " PDF Parameters: min, max")
    }
    else if(input$prior.pdf44 == "Normal"){
      current.label <- paste0(current.par[4], " PDF Parameters: mean, sigma")
    }
    tmp <- textInput("uncert_par44", label = current.label, value = "")
  }
  return(tmp)
})
#
output$uncert_sigma <- renderUI({
  if(is.null(server.env$current.ptest))
    return(NULL)
  if(!server.env$current.ptest$estimated)
    return(NULL)
  uncert_method <- input$uncertainty.method
  tmp <- NULL
  if(uncert_method == "amcmc"){
    tmp <- textInput(inputId = "uncert_sigma1", label = "Proposal sigma= ",
                     value = "0.1,0.1,0.1")
  }
  return(tmp)
})
#
output$uncert_run <- renderUI({
  if(is.null(server.env$current.ptest))
    return(NULL)
  if(!server.env$current.ptest$estimated)
    return(NULL)
  uncert_method <- input$uncertainty.method
#  tmp <- NULL
#  if(uncert_method == "amcmc" | uncert_method == "twalk"){
    tmp <- actionButton("uncert_run1","Run UQ, Run", icon = icon("bullseye"))
#  }
  return(tmp)
})
#
uncert.results <- reactive({
  if(input$uncert_run1 == 0){
    return(NULL)
  }
  else{
    print(input$uncert_run1)
  }
  input$uncert_run1
  #print('uncert.results INSIDE')
  tmp <- NULL
  current.ptest <- server.env$current.ptest
  #print(c(input$uncertainty.method,current.ptest$estimated))
  if(input$uncertainty.method == "None")
    return(NULL)
  if(is.null(current.ptest))
    return(NULL)
  if(!current.ptest$estimated)
    return(NULL)
  #if(is.null(input$boot_level)){
    #print('null.boot_level')
  #  return(NULL)
  #}

  if(input$uncertainty.method == "bootstrap"){
    current.level <- as.numeric(unlist(strsplit(input$boot_level,",")))
    #
    cat("UQ Method:", input$uncertainty.method, "\n")
    current.ptest.ci <- confint(current.ptest, level = current.level,
                                method = 'bootstrap',
                                d = as.numeric(input$boot_d),
                                neval = as.numeric(input$uncert_nreal1),
                                seed = as.integer(input$uncert_seed1))
    #
    #print(current.ptest.ci)
    #print(current.ptest$hydraulic_parameters)
    hydraulic.parameters(current.ptest) <- current.ptest.ci$hydraulic.parameters[]
    hydraulic.parameter.names(current.ptest) <- current.ptest.ci$hydraulic.parameters.names
    current.hydr.par <- as.data.frame(current.ptest.ci$hydraulic.parameters)
    current.hydr.par.ci <- current.ptest.ci$hydraulic.parameters.ci
    names(current.hydr.par) <- hydr.par1[[current.ptest$model]]
    server.env$current.hydr.par.ci <- current.hydr.par.ci
    server.env$current.ptest <- current.ptest
    tmp <- server.env$current.hydr.par.ci
    return(tmp)
  }
  else if(input$uncertainty.method == "amcmc"){
    cat("UQ Method:", input$uncertainty.method, "\n")
    nhydrpar <- length(hydr.par.mcmc[[current.ptest$model]])
    #print(c("amcmc", nhydrpar))
    # READ prior.pdf and prior.parameters!!!!!
    prior.pdf <- vector('character', length = nhydrpar)
    prior.pdf[1] <- pdf.models[[input$prior.pdf11]]
    prior.pdf[2] <- pdf.models[[input$prior.pdf22]]
    #
    prior.parameters <- matrix(0.0, nhydrpar, 2)
    prior.parameters[1,] <- as.numeric(unlist(strsplit(input$uncert_par11,",")))
    prior.parameters[2,] <- as.numeric(unlist(strsplit(input$uncert_par22,",")))
    #
    if(nhydrpar == 3 ){
      prior.pdf[3] <- pdf.models[[input$prior.pdf33]]
      prior.parameters[3,] <- as.numeric(unlist(strsplit(input$uncert_par33,",")))
    }
    if(nhydrpar == 4){
      prior.pdf[3] <- pdf.models[[input$prior.pdf33]]
      prior.parameters[3,] <- as.numeric(unlist(strsplit(input$uncert_par33,",")))
      prior.pdf[4] <- pdf.models[[input$prior.pdf44]]
      prior.parameters[4,] <- as.numeric(unlist(strsplit(input$uncert_par44,",")))
    }
    #
    iterations <- as.numeric(input$uncert_nreal1)
    proposal.sigma <- as.numeric(unlist(strsplit(input$uncert_sigma1,",")))
    iter_update_par <- as.numeric(input$uncert_iter_update_cov1)
    cov.corr <- input$uncert_cov.corr1
    #
    res.amcmc <- fit.sampling(ptest = current.ptest, model = current.ptest$model,
                               method = "amcmc", prior.pdf = prior.pdf,
                               prior.parameters = prior.parameters,
                               iterations = iterations,
                               proposal.sigma = proposal.sigma,
                               iter_update_par = iter_update_par,
                               cov.corr = cov.corr)
    #
    burnIn <- as.numeric(input$uncert_burnin1)
    #
    hydraulic.parameters(current.ptest) <- res.amcmc$hydraulic.parameters[burnIn:iterations,]
    #
    #
    current.ptest$hydraulic_parameters_names <- res.amcmc$hydraulic.parameters.names
    estimated(current.ptest) <- TRUE
    current.hydr.par.ci <- res.amcmc$hydraulic.parameters.ci
    server.env$current.hydr.par.ci <- current.hydr.par.ci
    server.env$current.ptest <- current.ptest
    #print(current.ptest$hydraulic_parameters)
    #print(names(current.ptest))
    tmp <- server.env$current.hydr.par.ci
  }
  else if(input$uncertainty.method == "twalk"){
    cat("UQ Method:", input$uncertainty.method, "\n")
    nhydrpar <- length(hydr.par.mcmc[[current.ptest$model]])
    #print(nhydrpar)
    # READ prior.pdf and prior.parameters!!!!!
    prior.pdf <- vector('character', length = nhydrpar)
    prior.pdf[1] <- pdf.models[[input$prior.pdf11]]
    prior.pdf[2] <- pdf.models[[input$prior.pdf22]]
    #
    prior.parameters <- matrix(0.0, nhydrpar, 2)
    prior.parameters[1,] <- as.numeric(unlist(strsplit(input$uncert_par11,",")))
    prior.parameters[2,] <- as.numeric(unlist(strsplit(input$uncert_par22,",")))
    #
    if(nhydrpar == 3 ){
      prior.pdf[3] <- pdf.models[[input$prior.pdf33]]
      prior.parameters[3,] <- as.numeric(unlist(strsplit(input$uncert_par33,",")))
    }
    if(nhydrpar == 4){
      prior.pdf[3] <- pdf.models[[input$prior.pdf33]]
      prior.parameters[3,] <- as.numeric(unlist(strsplit(input$uncert_par33,",")))
      prior.pdf[4] <- pdf.models[[input$prior.pdf44]]
      prior.parameters[4,] <- as.numeric(unlist(strsplit(input$uncert_par44,",")))
    }
    #
    iterations <- as.numeric(input$uncert_nreal1)
    res.twalk <- fit.sampling(ptest = current.ptest, model = current.ptest$model,
                              method = "twalk", prior.pdf = prior.pdf,
                              prior.parameters = prior.parameters,
                              iterations = iterations)
    #
    burnIn <- as.numeric(input$uncert_burnin1)
    hydraulic.parameters(current.ptest) <- res.twalk$hydraulic.parameters[burnIn:iterations,]
    current.ptest$hydraulic_parameters_names <- res.twalk$hydraulic.parameters.names
    estimated(current.ptest) <- TRUE
    current.hydr.par.ci <- res.twalk$hydraulic.parameters.ci
    server.env$current.hydr.par.ci <- current.hydr.par.ci
    server.env$current.ptest <- current.ptest
    tmp <- server.env$current.hydr.par.ci
  }
  return(tmp)
})
#
observeEvent(input$uncert_run1, {
  #print('INSIDE')
  tmp <- uncert.results()
  #print(tmp)
})
#
output$uncertainty_results <- renderPlot({
  input$uncert_run1
  #print(c('uncert_results',input$uncert_run1))
  tmp <- NULL
  current.ptest <- server.env$current.ptest
  validate(
    need(current.ptest, "The Pumping Test has not been defined")
    #need(current.ptest$estimated, "A Pumping Test with the estimated parameters is required")
  )
  if(is.null(input$uncert_run1) || input$uncert_run1 == 0)
    return(NULL)
  if(input$uncertainty.method == "None")
    return(NULL)
  if(is.null(current.ptest))
    return(NULL)
  if(!current.ptest$estimated)
    return(NULL)
  if(input$uncertainty.method == "bootstrap"){
    if(is.null(input$boot_level)){ return(NULL) }
  }
  if(class(current.ptest$hydraulic_parameters) != "matrix")
    return(NULL)
  if(input$uncertainty.method == "bootstrap"){
    tmp <- plot(current.ptest, type = "uncertainty")
  }
  else if(input$uncertainty.method == "amcmc"){
    tmp <- plot(current.ptest, type = "uncertainty")
  }
  else if(input$uncertainty.method == "twalk"){
    tmp <- plot(current.ptest, type = "uncertainty")
  }
  return(tmp)
})
#
output$uncertainty_table <- renderTable({
  #if(is.null(server.env$current.hydr.par.ci))
  #  return(NULL)
  if(is.null(input$uncert_run1) || input$uncert_run1 == "0")
    return(NULL)
  current.ci <- server.env$current.hydr.par.ci
  tmp <- current.ci
}, rownames = TRUE, digits = -5)
#

#########################################################################################
#                         Drawdown Calculator Tab
#########################################################################################
#
# Calculate u
#
calculate_u <- function(){
  if(input$calc_model == "None")
    return(NULL)
  Tr <- isolate(as.numeric(input$calc_transmissivity))
  Ss <- isolate(as.numeric(input$calc_storage))
  r <- isolate(as.numeric(input$calc_distance))
  tt <- isolate(as.numeric(input$calc_time))
  u <- (Ss*r**2)/(4*Tr*tt)
  return(u)
}
#
# Calculate well function
#
calculate_w <- function(){
  if(input$calc_model == "None")
    return(NULL)
  result <- NULL
  u <- calculate_u()
  td <- 1.0/u
  Q1 <- isolate(as.numeric(input$calc_pumprate))
  r1 <- isolate(as.numeric(input$calc_distance))
  ptest <- pumping_test("Test", Q = Q1, r = r1,
                        t = c(1e-10, 1e-10), s = c(1e-10, 1e-10))
  cmd1 <- paste("W<-", input$calc_model,"_well_function(td,par)", sep ="")
  par <- list(coeffs = ptest$coeffs)
  Tr <- isolate(as.numeric(input$calc_transmissivity))
  calc_model <- isolate(input$calc_model)
  if(calc_model == "boulton"){
    par$sigma <- isolate(as.numeric(input$calc_storage)/as.numeric(input$calc_omegad))
    par$phi <- isolate(as.numeric(input$cal_phi))
  }
  #
  if(calc_model == "hantush_jacob"){
    Ka <- isolate(as.numeric(input$calc_Ka))
    B <- isolate(as.numeric(input$calc_e))
    par$beta <- sqrt(Tr*B/Ka)
  }
  #
  if(calc_model == "general_radial_flow"){
    rw <- isolate(as.numeric(input$calc_grf_rw))
    rc <- isolate(as.numeric(input$calc_grf_rc))
    par$rd <- rw/rc
    par$n <- isolate(as.numeric(input$calc_grf_n))
  }
  #
  #print(par)
  #print(cmd1)
  eval(parse(text = cmd1))
  return(W)
}
#
observeEvent(input$importptest != 0, {
  #create_ptest()
  #print(input$type)
  if(input$type != "Pumping")
    return(NULL)
  # Fill all the textboxes with the information of the current pumpingtest
  if(!is.null(server.env$current.ptest)){
    current.ptest <- server.env$current.ptest
    if(is.null(current.ptest$model))
      return(NULL)
    if(!current.ptest$estimated)
      return(NULL)
    Q <- current.ptest$Q
    r <- current.ptest$r
    Tr <- current.ptest$hydraulic_parameters$Tr
    Ss <- current.ptest$hydraulic_parameters$Ss
    updateTextInput(session, inputId = "calc_pumprate", value = as.character(Q))
    updateTextInput(session, inputId = "calc_distance", value = as.character(current.ptest$r))
    updateTextInput(session, inputId = "calc_transmissivity", value = as.character(Tr))
    updateTextInput(session, inputId = "calc_storage", value = as.character(Ss))
    updateSelectInput(session, inputId = "calc_model", selected = current.ptest$model)
  }
})
#
calculate_drawdown <- function() {
  #input$calculatebutton
  #print(input$calculatebutton)
  if(input$calc_model == "None")
    return(NULL)
  Q1 <- isolate(as.numeric(input$calc_pumprate))
  r1 <- isolate(as.numeric(input$calc_distance))
  ptest <- pumping_test("Test", Q = Q1, r = r1,
                        t = c(1e-10, 1e-10), s = c(1e-10,1e-10))
  hydr.par <- list(Tr = isolate(as.numeric(input$calc_transmissivity)),
                   Ss = isolate(as.numeric(input$calc_storage)))
  calc_model <- isolate(input$calc_model)
  if(calc_model == "boulton"){
    hydr.par$omegad <- isolate(as.numeric(input$calc_omegad))
    hydr.par$phi <- isolate(as.numeric(input$calc_phi))
  }
  else if(calc_model == "hantush_jacob"){
    hydr.par$Ka <- isolate(as.numeric(input$calc_Ka))
    ptest$additional_parameters <- list(B = isolate(as.numeric(input$calc_e)))
  }
  else if(calc_model == "papadopulos_cooper"){
    ptest$additional_parameters <- list(rw = isolate(as.numeric(input$calc_pc_rw)),
                                        rc = isolate(as.numeric(input$calc_pc_rc)))
    hydr.par$cd <- isolate(as.numeric(input$calc_pc_cd))
  }
  else if(calc_model == "general_radial_flow"){
    ptest$additional_parameters <- list(rw = isolate(as.numeric(input$calc_grf_rw)),
                                        rc = isolate(as.numeric(input$calc_grf_rc)))
    hydr.par$n <- isolate(as.numeric(input$calc_grf_n))
  }
  cmd1 <- paste("statpar<-", calc_model,"_calculate_parameters(ptest,hydr.par,FALSE)", sep ="")
  #print(cmd1)
  eval(parse(text = cmd1))
  parameters <- model.par[[calc_model]]
  #print(parameters)
  spar <- parameters
  tt <- isolate(as.numeric(input$calc_time))
  current_cmd <- paste("drawdown<-", calc_model, "_solution(ptest,", sep = "")
  #print(current_cmd)
  for(i in spar){
    current_cmd <- paste(current_cmd, "statpar$",i,",", sep = "")
  }
  current_cmd <- paste(current_cmd,"tt)", sep = "")
  eval(parse(text = current_cmd))
}
#
# output$calc_u <- renderText({
#   validate(
#     need(input$calc_model != "None", "All the information is required for the calculation")
#   )
#   if(input$calc_model == "None")
#     return(NULL)
#   u <- calculate_u()
#   u <- format(u, digits = 5, scientific = TRUE)
#   paste0("<h3>Drawdown Results</h3><br><b>u (Dimensionless Variable)= </b>", u)
# })
# #
# output$calc_w <- renderText({
#   if(input$calc_model == "None") return(NULL)
#   w <- calculate_w()
#   w <- format(w, digits = 5, scientific = TRUE)
#   paste0("<b>W (Well Function)= </b>", w)
#
# })
# #
# output$calc_drawdown <- renderText({
#   if(input$calc_model == "None") return(NULL)
#   dr <- calculate_drawdown()
#   dr <- format(dr, digits = 5)
#   paste0("<b>Drawdown(m)= </b>", dr)
# })

#
observeEvent(input$calculatebutton, {
  #print("Button Pressed")
  output$calc_u <- renderText({
    u <- calculate_u()
    u <- format(u, digits = 5, scientific = TRUE)
    paste0("<h3>Drawdown Results</h3><br><b>u (Dimensionless Variable)= </b>", u)
  })
  #
  output$calc_w <- renderText({
    w <- calculate_w()
    w <- format(w, digits = 5, scientific = TRUE)
    paste0("<b>W (Well Function)= </b>", w)

  })
  #
  output$calc_drawdown <- renderText({
    dr <- calculate_drawdown()
    dr <- format(dr, digits = 5)
    paste0("<b>Drawdown(m)= </b>", dr)
  })
})
#
######################################################################################
#                        Report Generator Tab
######################################################################################
output$report.html.eng <- downloadHandler(
  filename = "report_eng.html",
  content = function(file) {
    tempReport <- file.path(tempdir(), "report_html_eng.Rmd")
    file.copy("report_html_eng.Rmd", tempReport, overwrite = TRUE)
    # Set up parameters to pass to Rmd document
    if(!is.null(server.env$current.ptest)){
      params <- list(current.ptest = server.env$current.ptest)
    }
    else{
      #params <- list(current.ves = NULL)
      return(NULL)
    }
    # Knit the document, passing in the `params` list, and eval it in a
    # child of the global environment (this isolates the code in the document
    # from the code in this app).
    rmarkdown::render(tempReport, output_file = file,
                      params = params,
                      envir = new.env(parent = globalenv())
    )
  }
)

})
