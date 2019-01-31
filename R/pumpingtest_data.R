#' Pumping test data from a confined aquifer.
#'
#' A dataset containing the time and drawdowns of a pumping test in a confined aquifer.
#' The information for this pumping test is:
#' \itemize{
#' \item Pumping rate (in m3/s): Q <- 1.3888e-2
#' \item Distance to pumping well (in m): r <- 250
#' }
#'
#' @format A data frame with 22 rows and 2 variables:
#' \describe{
#'   \item{t}{time, in seconds}
#'   \item{s}{drawdown, in meters}
#' }
#' @source C.W. Fetter, 2001, Applied Hydrogeology, Fourth Edition. Prentice Hall,
#' Upper Saddle River, 598 pp. Table 5.1, page 172
"theis"
#' estimation of a slug test with the Cooper Solution
#'
#' The information for this slug test is:
#' \itemize{
#' \item Radius of the well (in m): rw <- 0.071
#' \item Radius of the casing (in m): rc <- 0.025
#' }
#'
#' Localisation: Test performed in a monitoring well in Lincoln County, Kansas
#' Geology. The well was screened in a deltaic sequence consisting of mudstone
#' interbedded with very fine sandstone.
#'
#' @format A data frame with 69 rows and 2 variables:
#' \describe{
#'   \item{t}{time, in seconds}
#'   \item{s}{drawdown, in meters}
#' }
#' @source Butler James, 1998, The design, performance and analysis of slug tests.
#' Lewis Publisher. Table 5.2 pp. 61-62.
"cooper_slug"
#' estimation of a pumping test in a Confined aquifer with leakage using the
#' Hantush-Jacob solution
#'
#' The information for this pumping test is:
#' \itemize{
#' \item Pumping rate (in m3/s):  Q <- 6.309e-3
#' \item Distance to the pumping well (in m): r <- 3.048
#' \item Thickness of the aquitard (in m): B <- 6.096
#' }
#'
#' @format A data frame with 43 rows and 2 variables:
#' \describe{
#'   \item{t}{time, in seconds}
#'   \item{s}{drawdown, in meters}
#' }
#' @source Hall, P. Water well and aquifer test analysis 1996, Water Resources
#' Publications, LLC, 412 pp. Data from Fig 11.14 p.267-268
"hantush_jacob"
#' estimation of a pumping test in an aquifer with well-bore storage
#'
#' The information for this pumping test is:
#' \itemize{
#'  \item Radius of the well (in m): rw <- 0.6096
#'  \item Radius of the casing (in m): rc <- 0.6096
#'  \item Distance to pumping well (in m): r <- 3.048
#'  \item Pumping Rate (in m3/s):  Q <- 0.0050472
#' }
#'
#' @format A data frame with 46 rows and 2 variables:
#' \describe{
#'   \item{t}{time, in seconds}
#'   \item{s}{drawdown, in meters}
#' }
#' @source Hall, P. Water well and aquifer test analysis 1996, Water Resources
#' Publications, LLC, 412 pp. Pages 171-174.
"papadopulos_cooper"
#' estimation of a pumping test in an unconfined aquifer with the Boulton's
#' Solution
#'
#' The information of this pumping test is:
#' \itemize{
#' \item Pumping rate (in m3/s): Q <- 0.030
#' \item Distance to pumping well (in m): r <- 20
#' }
#'
#' estimation with Boulton (1963) model for unconfined aquifer of a pumping
#' test in coastal aquifer
#'
#' @format A data frame with 132 rows and 2 variables:
#' \describe{
#'   \item{t}{time, in seconds}
#'   \item{s}{drawdown, in meters}
#' }
#' @source G. de Marsily, cours DEA Paris 6, "Aquifere cotier de Nefza, Tunisie".
#' Piezometer A3
"boulton"
#' estimation of the Shut-in slug test (pulse) using Neuzil (1982) solution
#'
#' The information for this slug test is:
#' \itemize{
#' \item Radius of the well (in m): rw <- 0.067
#' \item Effective compressibility (in Pa-1): Ceff <- 2.723e-09
#' \item Volume of the test section (in m3): Vs <- 1.59
#' }
#'
#' @format A data frame with 20 rows and 2 variables:
#' \describe{
#'   \item{t}{time, in seconds}
#'   \item{s}{drawdown, in meters}
#' }
#' @source Batu, V., Aquifer Hydraulics: A Comprehensive Guide to Hydrogeologic
#' Data Analysis, John Wiley, New York, 1998. Example 13-2, Pages 689-692
"neuzil"
#' Constant head test
#'
#' The information for this pumping test is:
#'\itemize{
#'\item Discharge (in m3/s): Q <- c(28.142,0.084)
#'\item well radius (in m): r <- 0.084
#'}
#' @format A data frame with 19 rows and 2 variables:
#' \describe{
#'   \item{t}{time, in seconds}
#'   \item{q}{discharge, in meters3/second}
#' }
#' @source Lohman (1965) Geology and artesian water supply of the Grand Junction area,
#' Colorado: U.S. Geol. Survey Prof. Paper 451, 149p. Tables 6 and 7, well 28
"jacob_lohman"
#' Cooper Slug Test - Initial solution
#'
#' @format A data frame with 100 rows and 2 variables:
#' \describe{
#'   \item{dsmax}{dsmax}
#'   \item{cdmax}{cdmax}
#' }
"cooper"
#' Cooper Slug Test - Initial solution
#'
#' @format A data frame with 100 rows and 1 variable:
#' \describe{
#'   \item{td05}{td05}
#' }
"td05"
#' General Radial Flow: Example 1
#'
#' Data from well well Sha1B. The information of this pumping tests is:
#' \itemize{
#' \item Distance to pumping well (in m): r <- 26.2
#' \item Radius of the well (in m): rw <- 0.1
#' \item Radius of the casing (in m): rc <- 0.1
#' \item Pumping rate in(m3/2): Q <- 0.02322
#'}
#'
#' @format A data frame with 78 rows and 2 variables
#' \describe{
#' \item{t}{Time in seconds}
#' \item{s}{Drawdown in m}
#' }
#'
#' @source G. Lods and P. Gouze (2004) WTFM, software for well test analysis in
#' fractured media combining fractional flow with double porosity approaches.
#' Computers and Geosciences. Vol 30. pp. 937-947 Example CS2 page 943-944
"general_radial_flow1"
#' General Radial Flow: Example 2
#'
#' Data from well well Sha1BT. The information of this pumping tests is:
#' \itemize{
#' \item Distance to pumping well (in m): r <- 26.2
#' \item Radius of the well (in m): rw <- 0.1
#' \item Radius of the casing (in m): rc <- 0.1
#' \item Pumping rate in(m3/2): Q <- 0.02322
#'}
#'
#' @format A data frame with 78 rows and 2 variables
#' \describe{
#' \item{t}{Time in seconds}
#' \item{s}{Drawdown in m}
#' }
#'
#' @source G. Lods and P. Gouze (2004) WTFM, software for well test analysis in
#' fractured media combining fractional flow with double porosity approaches.
#' Computers and Geosciences. Vol 30. pp. 937-947 Example CS3 page 943-944
"general_radial_flow2"
#' estimation of a recovery test after a constant rate or variable rate test.
#'
#' The information of this recovery test is the following:
#' \itemize{
#' \item Duration of the pumping test (in s): t <- 14400
#' \item Pumping rate (in m3/s): Q <- 0.02893519
#' \item Distance to the pumping well (in m): r <- 60
#'}
#' @format A data frame with 15 rows and 2 variables
#' \describe{
#' \item{t}{Time in seconds}
#' \item{s}{Drawdown in m}
#' }
#'
#'@source Batu, V., Aquifer Hydraulics: A Comprehensive Guide to Hydrogeologic
#'Data Analysis, John Wiley, New York, 1998. Example 4-12, Pages 183-186
"agarwal_recovery"
#' estimation of Step-drawdown test with Eden and Hazel (1973) method.
#'
#' The pumping rates in m3/s are given by:
#' \itemize{
#' \item pumping rates:  qp <- c(0.01511574,0.01959491,0.02804398,0.03774306,0.04738426,0.05809028)
#' \item final time for each step: tp <- c(10800, 21600, 32400, 43200, 54000, 64800)
#'}
#'
#' @format A data frame with 175 rows and 2 variables
#' \describe{
#' \item{t}{time in seconds}
#' \item{s}{Drawdown in m}
#' }
#'
#'@source Clark, L. (1977) The analysis and planning of step-drawdown tests.
#'Quarterly Journal of Engineering Geology and Hydrogeology 10(2):125-143
"hazel"
#' Pumping test in a double porosity aquifer. estimation with the Warren and Root
#' model.
#'
#' Fissured tertiary volcanic rocks in the vicinity of the Yucca Mountain
#' at Nevada Test Site. One pumped well UE-25b#1 and one observation well UE-25a#1.
#' Pumping duration was about three days. The information of this pumping test is:
#' \itemize{
#' \item pumping rate (m3/s): Q <- 3.58e-2
#' \item Well radius (m) rw <- 0.11
#' \item Casing radius (m): rc <- 0.11
#' \item Aquifer Thickness (m) B=400
#'}
#'
#' @format A data frame with 72 rows and 2 variables
#' \describe{
#' \item{t}{time in seconds}
#' \item{s}{Drawdown in m}
#' }
#'
#'@source Moench, A. 1984. Double porosity model for a fissured groundwater
#'reservoir with fracture skin. Water Resources Research, 20(7), 831-846.
#'Table 2, Page 838.
"warren_root"
#' Single fracture in a confined aquifer. estimation using the Gringarten Solution.
#'
#'The information of this pumping test is:
#' \itemize{
#' \item pumping rate (m3/s): Q <- 7.7e-4
#' \item Storage Coefficient S <- 1
#' \item Porosity : 0.12
#'}
#'
#' @format A data frame with 26 rows and 2 variables
#' \describe{
#' \item{t}{time in seconds}
#' \item{s}{Drawdown in m}
#' }
#'
#'@source Gringarten A. C., Ramey H. J., & Raghavan R., 1975.
#'Applied Pressure Analysis for Fractured Wells, Petroleum Transactions,
#'AIME follows page 784.
"gringarten"
#' Data from a pumping test with variable ppumping rates
#'
#' The information of this pumping test is:
#' \itemize{
#' \item Distance to observation well (r): 5 m
#' \item Pumping rate history:
#' - Step 1: 300s-1800s, 0.00578 m3/s
#' - Step 2: 2100s-4800s, 0.0081 m3/s
#' - Step 3: 5400s-7800s, 0.0069 m3/s
#' }
#' @format A data frame with 18 rows and 2 variables
#' \describe{
#' \item{t}{time in seconds}
#' \item{s}{Drawdown in m}
#' }
#'
#'@source Kruseman and de Ridder. 1994. Analysis and Evaluation of pumping test
#'data. Second Edition pp. 185.
"variable_rate"