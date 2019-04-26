#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)

# Define UI for application that draws a histogram
shinyUI(pageWithSidebar(
  # Application title
  headerPanel("Analytical-GUI: Shiny Interface (v0.1)"),
  
  #### Panel 'About' (right hand side)
  ##############################################################################
  sidebarPanel(
    imageOutput("uptc.logo", inline=TRUE),
    selectInput('type', "Select the Test type:",
                c(Pumping="Pumping", Slug="Slug", Bailer = "Bailer")),
    p(HTML("<h5>This is Analytical-GUI, the Shiny interface for the evaluation of
           the analytical solutions of pumping and slug test data in <strong>R</strong>.</h5>
           This application can be used to understand several of the analytical
           solutions used in the description of the flow of water towards a well and their
           corresponding derivatives. This interface use the  <strong>R</strong>
           package  <a href='http://www.github.com/khaors/pumpingtest'>pumpingtest</a>.")),
    p(HTML('This package was developed by Oscar Garcia-Cabrejo, School of Geological
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
           <I> Multilayered Aquifier Systems: Fundamentals and Applications</I>
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
      tabPanel("Analytical Solution Explorer1",
               icon = icon("line-chart"),
               withMathJax(),
               h3("Explore analytical type curves and their derivatives"),
               br(),
               uiOutput(outputId = "curve0"),
               uiOutput(outputId = "curve1"),#Model
               uiOutput(outputId = "curve15"), #Scale
               uiOutput(outputId = "curve2"),
               uiOutput(outputId = "curve3"),
               uiOutput(outputId = "curve4"),
               uiOutput(outputId = "curve4.5"),
               uiOutput(outputId = "curve5"),
               plotOutput(outputId = "plot_curve"),
               uiOutput(outputId = "curve6"),
               uiOutput(outputId = "curve7")),
      #########################################################################################
      tabPanel("Compare Solutions",
               icon = icon("line-chart"),
               h3("Compare analytical type curves and their derivatives"),
               br(),
               p(HTML("<h5>On this tab it is possible to compare two drawdown curves</h5>")),
               br(),
               uiOutput(outputId = "compare_time"),
               uiOutput(outputId = "compare_scale"),
               uiOutput(outputId = "compare0"), #Dimensionless or dimension
               br(),
               h4("Curve 1"),
               uiOutput(outputId = "compare1a"),
               uiOutput(outputId = "compare1b"),
               uiOutput(outputId = "compare1c"),
               uiOutput(outputId = "compare1d"),
               br(),
               h4("Curve 2"),
               uiOutput(outputId = "compare2a"),
               uiOutput(outputId = "compare2b"),
               uiOutput(outputId = "compare2c"),
               uiOutput(outputId = "compare2d"),
               br(),
               plotOutput(outputId = "plot_comparison")
      )
      ##################################################################################
      #                     Drawdown Curve Calculator
      ##################################################################################
      # tabPanel("Drawdown Curve",
      #          icon = icon("calculator"),
      #          h3("Calculate and plot drawdown curves and their derivatives"),
      #          br(),
      #          br(),
      #          textInput("calc_pumprate", label = "Pumping Rate(m3/s)", value = "0.0"),
      #          textInput("calc_distance", label = "Distance(m)", value = "0.0"),
      #          #Time is a vector
      #          textInput("calc_time", label = "Time(s)", value = "0.0"),
      #          textInput("calc_transmissivity", label = "Transm.(m2/s)", value = "1e-16"),
      #          textInput("calc_storage", label = "Storage Coeff.", value = "1.0e-16"),
      #          #
      #          selectInput("calc_model", label = "Solution Type",
      #                      c("None", "theis", "cooper_jacob", "hantush_jacob",
      #                        "boulton", "general_radial_flow",
      #                        "papadopulos_cooper", "agarwal_skin",
      #                        "agarwal_recovery", "warren_root",
      #                        "gringarten"),
      #                      selected = "None"),
      #          # UI elements of additional parameters
      #          conditionalPanel(
      #            condition = "input.calc_model == 'boulton'",
      #            textInput("omegad", label = "Drainage Porosity", value = "1e-3"),
      #            textInput("phi", label = "Phi", value = "1e-3")
      #          ),
      #          #
      #          conditionalPanel(
      #            condition = "input.calc_model == 'hantush_jacob'",
      #            textInput("e",label = "Aquitard Thickness(m)", value = "0.0"),
      #            textInput("Ka", label = "Aquitard Conductivity (m/s)", value = "1.0e-16")
      #          ),
      #          #
      #          conditionalPanel(
      #            condition = "input.calc_model == 'general_radial_flow'",
      #            textInput("grf_n", label = "Flow Dimension", value = "2.0"),
      #            textInput("grf_rw", label = "Well Radius(m)", value = "0.0"),
      #            textInput("grf_rc", label = "Well Casing(m)", value = "0.0")
      #          ),
      #          conditionalPanel(
      #            condition = "input.calc_model == 'papadopulos_cooper'",
      #            textInput("pc_rw", label = "Well Radius(m)", value = "0.0"),
      #            textInput("pc_rc", label = "Well Casing(m)", value = "0.0")
      #          )),
      
      #######################################################################################
#      tabPanel("Sensitivity Analysis",
#               icon = icon("map"),
#               h3("Senstivity Analysis of analytical solution using different approaches"))
    )
    
  )
    ))
