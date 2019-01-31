#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(rhandsontable)
library(shinyalert)
# Define UI for application that draws a histogram
# Define UI for dataset viewer application
shinyUI(pageWithSidebar(
  # Application title
  headerPanel("PumpingTest-GUI: Shiny Interface (v0.1)"),
  #### Panel 'About' (right hand side)
  ##############################################################################
  sidebarPanel(width = 3,
    useShinyalert(),
    imageOutput("uptc.logo", inline=TRUE),
    #selectInput('type', "Select the Test type:",
    #            c(Pumping="Pumping", Slug="Slug"#,
                  #ConstantHead="ConstantHead",
                  #VariableDischarge = "VariableDischarge"
    #              )),
    p(HTML("<h5>This is PumpingTest-GUI, the Shiny interface for analysis and
            evaluation of pumping and slug test data in <strong>R</strong>.</h5>
            This application can be used for the identification of the type of
            aquifer, estimation of hydraulic parameters, diagnosis of the estimated
            model, uncertainty quantification of the estimated parameters, and
            evaluation of drawdown in space and time using the  <strong>R</strong>
            package  <a href='http://www.github.com/khaors/pumpingtest'>pumpingtest</a>.
            The analysis of a pumping/slug test is achieved in six simple steps using the
            panels on the right.")),
    p(HTML('This GUI was developed by Oscar Garcia-Cabrejo, School of Geological
            Engineering, Universidad Pedagogica y Tecnologica de Colombia, Sogamoso,
            Boyaca, Colombia. It is also included in the
           <strong>R</strong> package
           <a href="http://www.github.com/khaors/pumpingtest">pumpingtest</a>. Its source code
           is freely available on github.')),
    br(),

    h3('References:'),
    p(HTML('<li> <span style="font-variant: small-caps;">G. P. Kruseman & de Ridder, N.</span>
           (1992).<I>Analysis and Evaluation of Pumping Test Data</I>.
           International Institute for Land Reclamation and Improvement, The Netherlands.</li>
           <li> <span style="font-variant: small-caps;">Lebbe, L.</span>(1999).
           <I>Hydraulic Parameter Identification: Generalized Interpretation
           Method for Single and Multiple Pumping Tests</I>. Berlin/Heidelberg:
           Springer-Verlag.</li>
           <li> <span style="font-variant: small-caps;">Cheng, A.</span> (2000).
           <I> Multilayered Aquifer Systems: Fundamentals and Applications</I>
           CRC Press, Boca Raton, Fl.</li>
           <li> <span style="font-variant: small-caps;">Walton, W.</span>
           (2006). <I>Aquifer Test Modeling</I>, CRC Press, Boca Raton, Fl.</li>
           <li> <span style="font-variant: small-caps;">Sindalovskiy, L.</span>(2016).
           <I>Aquifer Test Solutions: A Practitioner’s Guide With Algorithms Using ANSDIMAT</I>
           Springer-Verlag, Heidelberg</li>'))
    ),

    # Show a plot of the generated distribution
    mainPanel(
      tabsetPanel(
        #### Panel 'Import Data'
        #########################################################################
        tabPanel("Import Data",
                 icon = icon("file"),
                 h3("First step: import data"),

                 p(HTML("To stat using the application, you must import the data using the options
                         included in this tabset. There are different options included to give ample
                         flexibility in importing the information from a pumping test. In addition
                         you can see a sample of the imported data in the lower part of the tabset.
                         When this is done, proceed to the next step: specify pumping test information.")),
                 p(HTML('This interface can be tested using data files for
                        the <a href="http://github.com/khaors/pumpingtestdata/confined"target="_blank">confined</a>, <a href=
                        "http://github.com/khaors/pumpingtestdata/phreatic"target="_blank">phreatic</a> and <a href=
                        "http://github.com/khaors/pumpingtestdata/leaky"target="_blank">leaky </a> aquifers
                        (download these files on your computer, load them on the interface and proceed).')),

                 br(),
                 fluidRow(
                   column(4,
                          fileInput('file1', 'Choose CSV/TXT File'),
                          checkboxInput('header', ' Header?', TRUE),
                          checkboxInput('rownames', ' Row names?', FALSE),
                          selectInput('sep', 'Separator:',
                                      c("Comma","Semicolon","Tab","Space"), 'Comma'),
                          selectInput('quote', 'Quote:',
                                      c("None","Double Quote","Single Quote"),
                                      'Double Quote'),
                          selectInput('dec', 'Decimal mark', c("Period", "Comma"),
                                      'Period'),
                          numericInput('nrow.preview','Number of rows in the preview:',20),
                          numericInput('ncol.preview', 'Number of columns in the preview:',
                                       10)
                   ),
                   column(3,
                          h4('Data Table'),
                          br(),
                          helpText("Note: Even if the preview only shows a restricted
                          number of observations, the pumping_test object will be created
                          based on the full dataset."),
                          tableOutput("view")
                          ),
                   column(4,
                          h4("Summary"),
                          br(),
                          verbatimTextOutput("summary"))
                 )
        ),
        #
        # Test Information
        #
        tabPanel("Test Information",
                 icon = icon("newspaper-o"),
                 h3("Second Step: Specify Test information"),

                 p(HTML("Here the Pumping or Slug test information is specified,  including
                                 test name, pumping rate, distance to the observation well, and other
                                 additional parameters such as aquitard thickness, well radius, casing
                                 radius, etc. After this hit the 'Define Aquifer Test' button and then
                                 you can proceed to the next step: Identify the aquifer type.")),
                 fluidRow(
                   column(4,
                          #
                          selectInput('type', "Select the Test type:",
                                      c(Pumping = "Pumping", Slug = "Slug",
                                        VariableRate = "VariableRate",
                                        ConstantDischarge = "ConstantHead")),

                          #
                          textInput("testname", label = "Test Name", value = "Test0"),
                          #
                          uiOutput("prate"),
                          rHandsontableOutput("vrate"),
                          br(),
                          uiOutput("dist"),
                          #
                          br(),
                          h4("Additional Parameters"),
                          br(),
                          uiOutput('addpar1'),
                          uiOutput('setpar1'),
                          uiOutput('addpar2'),
                          uiOutput('setpar2'),
                          uiOutput('addpar3'),
                          uiOutput('setpar3'),
                          uiOutput('addpar4'),
                          uiOutput('setpar4'),
                          uiOutput('addpar5'),
                          uiOutput('setpar5'),
                          #
                          br(),
                          actionButton("defineptest", label = "Define Aquifer Test",
                                       icon = icon("bullseye"))),
                   column(6,  offset = 1,
                          h4("Plot"),
                          plotOutput('plot.simple'),
                          br(),
                          actionButton(inputId = "applyfilter", label = "Apply filter",
                                       icon = icon('filter')),
                          br(),
                          h5("Filter"),
                          br(),
                          uiOutput('filter.par1'),
                          uiOutput('filter.par2'))
                  )

        ),
        ######################################################################################
        #                            Smoothing Data
        ######################################################################################
        # tabPanel("Smoothing Data",
        #          icon = icon('paint-brush'),
        #          h3("Second and a half Step (Optional): Smooth drawdown information"),
        #          br(),
        #          p(HTML("Here the Pumping or Slug test information is specified,  including
        #                          test name, pumping rate, distance to the observation well, and other
        #                          additional parameters such as aquitard thickness, well radius, casing
        #                          radius, etc. After this hit the 'Define Aquifer Test' button and then
        #                          you can proceed to the next step: Identify the aquifer type.")),
        #          br(),
        #          sidebarLayout(
        #            sidebarPanel(
        #              uiOutput("smoothmethod0")
        #            ),
        #            mainPanel(
        #              plotOutput("smoothing.plot", height = 500, width = 500*1.5)
        #            )
        #          )
        #       ),
        ######################################################################################
        #                            Diagnostic Plot
        ######################################################################################
        tabPanel("Diagnostic Plot",
                 icon = icon("bar-chart-o"),
                 h3("Third step: Identify the aquifer type"),

                 p(HTML("The aquifer type can be identified from the flow conditions that
                         prevailed during the pumping test using the Diagnostic Plot. This
                         is a plot of the drawdown and its derivative with respect to the
                         logarithm of time. Here, there are six types of derivatives that
                         can be used to identify the flow conditions (aquifer type) using
                         the collection of diagnostic plots show in the lower part.")),
                 sidebarLayout(
                   sidebarPanel(
                     uiOutput("derivative0"),
                     #selectInput("derivative", "Derivative Type", c("None", "central", "horner",
                     #                                                "bourdet","spline","spane",
                     #                                                "smoothspline"),
                      #           selected = "None"),
                     conditionalPanel(
                       condition = "input.derivative == 'bourdet'",
                       sliderInput("points",
                                   "Bourder: Number of points",
                                   min = 2,
                                   max = 10,
                                   value = 5)
                     ),
                     conditionalPanel(
                       condition = "input.derivative == 'spline' ",
                       numericInput("sppoints", "Spline: Number of points", 20, min = 1, max = 40, step = 1)
                     ),
                     conditionalPanel(
                       condition = "input.derivative == 'spane' ",
                       numericInput("sppoints1", "LM: Number of points", 2, min = 1, max = 40, step = 1)
                     ),
                     conditionalPanel(
                       condition = "input.derivative == 'kernelreg'",
                       textInput(inputId = "krbw", label = "KernelReg: Bandwidth", value = "0.1")
                     ),
                     radioButtons("scale", "Plot Scale", c("loglog","slog"), selected = "loglog")
                   ),
                   mainPanel(
                     plotOutput("diagnostic.plot", height = 500, width = 500*1.5),
                     wellPanel(
                       imageOutput("diagnostic.table",  height = 500, width = 500*1.5,
                                   inline=FALSE)
                     )
                   )
                 )
                 ),
        ######################################################################################
        #                              Parameter estimation
        ######################################################################################
        tabPanel("Parameter Estimation",
                 icon = icon("puzzle-piece"),#signal
                 h3("Fourth Step: Estimate the hydraulic parameters"),

                 p(HTML("Once the aquifer model is identified, the hydraulic parameters can
                        be estimated. Pumpingtest includes three different methods for parameter
                        estimation: visual and automatic (via optimization).
                        ")),
                 br(),
                 br(),
                 downloadButton("pumpingtest.download", "Download the PumpingTest file"),

                 br(),
                 br(),
                 sidebarLayout(
                   sidebarPanel(width = 4,
                                selectInput(inputId = 'parest_method', "Select the Parameter Estimation Method:",
                                            c(None = "None", visual = "visual", automatic ="automatic")),
                                uiOutput(outputId = "parest_opt"),
                                uiOutput(outputId = "parest_objfn"),
                                uiOutput(outputId = "parest_opt_par_flag"),
                                uiOutput(outputId = "parest_opt_par"),
                                # Add ui Elements for hydraulic parameters
                                uiOutput(outputId = "parest_model"),
                                uiOutput(outputId = "parest_model_text"),
                                uiOutput(outputId = "parest_derivative"),
                                uiOutput(outputId = "parest_dpoints"),
                                uiOutput(outputId = "parest_scale"),
                                uiOutput(outputId = "parest_logTr"),
                                uiOutput(outputId = "parest_logSs"),
                                uiOutput(outputId = "parest_omegad"),
                                uiOutput(outputId = "parest_phi"),
                                uiOutput(outputId = "parest_pc_cd"),
                                uiOutput(outputId = "parest_hj_ka"),
                                uiOutput(outputId = "parest_grf_n"),
                                uiOutput(outputId = "parest_wr_LogSm"),
                                uiOutput(outputId = "parest_wr_lambda")
                   ),
                   mainPanel(
                     plotOutput("estimation.plot", height = 500, width = 500*1.5),
                     uiOutput(outputId = "parest_results"),
                     dataTableOutput(outputId = "ci_parest_results")
                   )
                 )),
        ######################################################################################
        #                             Model diagnostic
        ######################################################################################
        tabPanel("Model Diagnostic",
                 icon = icon("wrench"),
                 h3("Fifth Step: Evaluate the estimated model"),
                 #
                 p(HTML("Once the hydraulic parameters of the identified model are estimated then
                        it is critical to evaluate the estimated model to check if the assumptions
                        made during the parameter estimation are valid. In addition, it is possible
                        to identify the samples that are influential for each one of the estimated
                        parameters")),
                 sidebarLayout(
                   sidebarPanel(
                     width = 4,
                     selectInput(inputId = "diagnostic.type", label = "Select the Diagnostic Type",
                                 choices = c("None", "Model Diagnostic", "Sample Influence"), selected = "None")
                   ),
                   mainPanel(
                     # Add Plot
                     plotOutput("model.diagnostic", height = 500, width = 500*1.5)
                   )
                 )
               ),
        ######################################################################################
        #                        Uncertainty Quantification
        ######################################################################################
        tabPanel("Uncertainty Quantification",
                 icon = icon("random"),
                 h3("Sixth Step: Quantify the Uncertainty in the hydraulic parameters"),
                 br(),
                 p(HTML("The hydraulic parameters are estimated from the drawdown measurements and
                        therefore the noise present in the available information is propagated to the
                        estimates.")),
                 br(),
                 sidebarLayout(
                   sidebarPanel(
                     selectInput('uncertainty.method', "Select the UQ Method:",
                                 c(None = "None", bootstrap = "bootstrap", amcmc ="amcmc",
                                   twalk = "twalk")),
                     # Bootstrapping options: level, d, nreal, seed
                     # MCMC options: prior.pdf, prior.parameters, proposal.sigma,
                     #               iter, burIn, cov.corr, iter.update.cov, seed
                     uiOutput("uncert_level"),
                     uiOutput("uncert_d"),
                     uiOutput("uncert_prior.pdf1"),
                     uiOutput("uncert_par1"),
                     uiOutput("uncert_prior.pdf2"),
                     uiOutput("uncert_par2"),
                     uiOutput("uncert_prior.pdf3"),
                     uiOutput("uncert_par3"),
                     uiOutput("uncert_prior.pdf4"),
                     uiOutput("uncert_par4"),
                     #uiOuptut("uncert_prior.parameters"),
                     uiOutput("uncert_sigma"),
                     uiOutput("uncert_nreal"),
                     uiOutput("uncert_burnin"),
                     uiOutput("uncert_cov.corr"),
                     uiOutput("uncert_iter_update_cov"),
                     uiOutput("uncert_seed"),
                     uiOutput("uncert_run")
                   ),
                   mainPanel(
                     plotOutput("uncertainty_results", height = 500, width = 500*1.5),
                     br(),
                     h4("Hydraulic Parameters: Confidence Intervals"),
                     br(),
                     tableOutput("uncertainty_table")
                   ))
                 ),
        ######################################################################################
        #                          Drawdown Calculator
        ######################################################################################
        tabPanel("Drawdown Calculator",
                 icon = icon("calculator"),

                 h3("Optional Step: Calculate drawdown using the estimated Parameters"),
                 br(),
                 actionButton("importptest", "Import Current PumpingTest", icon = icon("file")),
                 br(),
                 br(),
                 textInput("calc_pumprate", label = "Pumping Rate(m3/s)", value = "0.0"),
                 textInput("calc_distance", label = "Distance(m)", value = "0.0"),
                 textInput("calc_time", label = "Time(s)", value = "0.0"),
                 textInput("calc_transmissivity", label = "Transm.(m2/s)", value = "1e-16"),
                 textInput("calc_storage", label = "Storage Coeff.", value = "1.0e-16"),
                 #
                 selectInput("calc_model", label = "Solution Type",
                             c("None", "theis", "cooper_jacob", "hantush_jacob",
                               "boulton", "general_radial_flow",
                               "papadopulos_cooper", "agarwal_skin",
                               "agarwal_recovery", "warren_root",
                               "gringarten"),
                             selected = "None"),
                 # UI elements of additional parameters
                 conditionalPanel(
                   condition = "input.calc_model == 'boulton'",
                   textInput("calc_omegad", label = "Drainage Porosity", value = "1e-3"),
                   textInput("calc_phi", label = "Phi", value = "1e-3")
                 ),
                 #
                 conditionalPanel(
                   condition = "input.calc_model == 'hantush_jacob'",
                   textInput("calc_e",label = "Aquitard Thickness(m)", value = "0.0"),
                   textInput("calc_Ka", label = "Aquitard Conductivity (m/s)", value = "1.0e-16")
                 ),
                 #
                 conditionalPanel(
                   condition = "input.calc_model == 'general_radial_flow'",
                   textInput("calc_grf_n", label = "Flow Dimension", value = "2.0"),
                   textInput("calc_grf_rw", label = "Well Radius(m)", value = "0.0"),
                   textInput("calc_grf_rc", label = "Well Casing(m)", value = "0.0")
                 ),
                 conditionalPanel(
                   condition = "input.calc_model == 'papadopulos_cooper'",
                   textInput("calc_pc_rw", label = "Well Radius(m)", value = "0.0"),
                   textInput("calc_pc_rc", label = "Well Casing(m)", value = "0.0")
                 ),
                 br(),
                 actionButton("calculatebutton","Calculate", icon = icon("bullseye")),
                 #br(),
                 uiOutput('calc_u'),
                 uiOutput('calc_w'),
                 uiOutput('calc_drawdown')
                 ),
        #################################################################################
        #                               Report Generator Panel
        #################################################################################
        tabPanel("Report Generator",
                 icon = icon("book"),
                 h2("Reports"),
                 br(),
                 "This tab allows the automatic generation and download of some preconfigured reports with
  tables and graphs in English.",
                 br(),
                 br(),br(),
                 fluidRow(
                   downloadButton(outputId = "report.html.eng", label = "Generate report: html - English"),
                   br(),
                   downloadButton(outputId = "report.word.eng", label = "Generate report: word - English")
                 )

        )
        ######################################################################################
        #                                 Help Panel
        ######################################################################################
#        tabPanel("Help",
#                 icon = icon("map-o"))

        )
    )
  )
)
