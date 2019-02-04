#' @section Aquifer functions:
#' The aquifer functions include several functions designed to represent 
#' several types of aquifers including:
#'  
#' @docType package
#' @name pumpingtest
NULL
#' @title 
#' aquifer
#' @description 
#' Function to create an aquifer object
#' @param name Name of the aquifer
#' @param solution Analytical solution to be used
#' @param type Type of aquifer to be defined. Currently the supported types include:
#' \itemize{
#' \item Infinite
#' \item semi-infinite: with no-flow or constant head boundaries
#' \item strip:  with no-flow or constant head boundaries or a combination of both
#' \item wedge
#' \item U-shaped
#' \item rectangular
#' \item circular} 
#' @param xlim Numeric vector with limits of the aquifer in the x direction
#' @param ylim Numeric vector with limits of the aquifer in the y direction
#' @param nx Integer with the number of nodes in the x direction
#' @param ny Integer with the number of nodes in the y direction
#' @param hydraulic.pars A list with the hydraulic parameters of the aquifer
#' @return 
#' This function returns an aquifer object
#' @author 
#' Oscar Garcia-Cabrejo, \email{khaors@gmail.com}
#' @family aquifer functions
#' @export
aquifer <- function(name, solution, type, xlim, ylim, nx, ny, hydraulic.pars){
  res <- list(name = name,
              solution = solution,
              type = type, 
              xlim = xlim, 
              ylim = ylim, 
              nx = nx, 
              ny = ny, 
              grid = NULL, coords = NULL, 
              hydraulic.parameters = hydraulic.pars,
              trend = NULL, grid_sp = NULL, grid_sp1 = NULL, 
              number.wells = NULL,
              well.discharge = vector("numeric", length = 10),
              well.coords.x = vector("numeric", length = 10), 
              well.coords.y = vector("numeric", length = 10))
  coords <- calculate_aquifer_coordinates(nx, ny, xlim[1], xlim[2], 
                                          ylim[1], ylim[2])
  res$xcoords <- coords[,1]
  res$ycoords <- coords[,2]
  xsiz <- range(xlim)[2]/nx
  ysiz <- range(ylim)[2]/ny
  res$number.wells <- 0
  # boundaries: use points or lines?
  class(res) <- "aquifer"
  invisible(res)
  return(res)
}
#' @title
#' summary.aquifer
#' @description 
#' Function to print a summary of an aquifer object
#' @param object An aquifer object
#' @param ... Additional arguments
#' @return 
#' This function prints a summary of the aquifer object on the screen
#' @author 
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @family aquifer functions
#' @export 
summary.aquifer <- function(object, ...){
  if(class(object) != "aquifer"){
    stop("An aquifer object is required as input")
  }
}  
#' @title 
#' print.aquifer
#' @description 
#' Function to print the aquifer object on the console
#' @param x An aquifer object
#' @param ... Additional parameters
#' @return 
#' This function prints all the elements of the aquifer on the console
#' @author 
#' Oscar Garcia-Cabrejo, \email{khaors@gmail.com}
#' @family aquifer functions
#' @export
print.aquifer <- function(x, ...){
  if(class(x) != "aquifer"){
    stop("An aquifer object is required as input")
  }
  print("Aquifer")
  print(paste0("Name = ", x$name))
  print(paste0("Type= ", x$type))
  print(paste0("X Dim ", as.character(x$xmin), ',' ,as.character(x$xmax)))
  print(paste0("Y Dim ", as.character(x$xmin), ',', as.character(x$ymax)))
  print(paste0("nx,ny"), as.character(x$nx), ',',  as.character(x$ny))
}
#' @title
#' add.well<-
#' @description 
#' Function to add wells to an aquifer object
#' @param x An aquifer object
#' @param value A list with the values of the well discharges, x coordinates 
#' of the wells and the y coordinates of the wells
#' @return 
#' This function adds a well to an aquifer
#' @author
#' Oscar Garcia-Cabrejo, \email{khaors@gmail.com}
#' @family aquifer functions
#' @export 
`add.well<-` <- function(x, value) {
  if(class(x) != 'aquifer'){
    stop("An aquifer object is required as input")
  }
  x$well.discharge <- value$Q
  x$well.coords.x <- value$x0
  x$well.coords.y <- value$y0
  #
  return(x)
}
#' @title 
#' calculate_drawdown
#' @description 
#' Function to calculate the drawdown produced by well withdrawals at given times
#' @param aq An aquifer object
#' @param current.times A numeric vector with the current times
#' @return 
#' This function returns a matrix where the calculated drawdowns are assigned in 
#' the columns
#' @author 
#' Oscar Garcia-Cabrejo, \email{khaors@gmail.com} 
#' @family aquifer functions
calculate_drawdown <- function(aq, current.times){
  if(class(aq) != 'aquifer'){
    stop("An aquifer object is required as input")
  }
  #
  if(class(current.times) != "numeric"){
    stop("A numeric vector is required as input")
  }
  #
  res <- NULL
  if(aq$type == 'infinite'){
    res <- infinite_aquifer_calculate_drawdown_cpp(model = aq$type,
                                                   Q = aq$Q, 
                                                   x0 = aq$well_coordinates_x, 
                                                   y0 = aq$well_coordinates_y,
                                                   t = current.times, 
                                                   hydrpar = aq$hydraulic.parameters, 
                                                   nx = aq$nx, 
                                                   ny = aq$ny, 
                                                   xmin = aq$xmin, 
                                                   xmax = aq$xmax, 
                                                   ymin = aq$ymin, 
                                                   ymax = aq$ymax)
  }
  else{
    stop("Aquifer types is not currently supported")
  }
  #
  
  #
  return(res)
}