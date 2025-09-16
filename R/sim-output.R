#' @title SimOutput
#' @description A class to store simulation output information.
#' @field name Name of the simulation output
#' @field path Path to the simulation output
#' @field dimension Dimension of the simulation output
#' @field displayUnit Display unit of the simulation output
#' @export
SimOutput <- R6::R6Class(
  "SimOutput",
  public = list(
    
    #' @description Initialize the simulation output object
    #' @param name Name of the simulation output
    #' @param path Path to the simulation output
    #' @param dimension Dimension of the simulation output
    #' @param displayUnit Display unit of the simulation output
    initialize = function(name,
                          path,
                          dimension,
                          displayUnit) {
      private$.name <- name
      private$.path <- path
      private$.dimension <- dimension
      private$.displayUnit <- displayUnit
    }
  ),
  active = list(
    name = function(value) {
      if (missing(value)) {
        return(private$.name)
      }
    },
    path = function(value) {
      if (missing(value)) {
        return(private$.path)
      }
    },
    dimension = function(value) {
      if (missing(value)) {
        return(private$.dimension)
      }
    },
    displayUnit = function(value) {
      if (missing(value)) {
        return(private$.displayUnit)
      }
    }
  ),
  private = list(
    .name = NULL,
    .path = NULL,
    .dimension = NULL,
    .displayUnit = NULL
  )
)
