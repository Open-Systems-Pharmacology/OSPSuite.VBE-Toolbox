#' VBEParameter Class
#'
#' An R6 class for handling Virtual Bioequivalence (VBE) parameters.
#'
#' @export
VBEParameter <- R6::R6Class(
  classname = "VBEParameter",
  public = list(
    #' @description
    #' Initialize a VBEParameter object.
    #'
    #' @param pathInReferenceSimulation Character. Path to the parameter in the reference simulation.
    #' @param pathInTestSimulation Character. Path to the parameter in the test simulation.
    #' @param displayName Character. Display name of the parameter (defaults to `pathInReferenceSimulation`).
    #' @param dimension Character. Dimension of the parameter.
    #' @param unit Character. Unit of the parameter.
    #' @param lowerBound Numeric. Lower bound for the parameter (default: 1e-1).
    #' @param upperBound Numeric. Upper bound for the parameter (default: 1e1).
    #' @param logIncrement Logical. Whether the parameter is incremented logarithmically (default: TRUE).
    initialize = function(pathInReferenceSimulation,
                          pathInTestSimulation,
                          displayName = NULL,
                          dimension,
                          unit,
                          lowerBound = 1e-1,
                          upperBound = 1e1,
                          logIncrement = TRUE) {
      self$pathInReferenceSimulation <- pathInReferenceSimulation
      self$pathInTestSimulation <- pathInTestSimulation
      self$displayName <- displayName %||% pathInReferenceSimulation
      self$dimension <- dimension
      self$unit <- unit
      self$lowerBound <- lowerBound
      self$upperBound <- upperBound
      self$logIncrement <- logIncrement
    }
  ),
  private = list(
    .pathInReferenceSimulation = NULL,
    .pathInTestSimulation = NULL,
    .displayName = NULL,
    .dimension = NULL,
    .unit = NULL,
    .lowerBound = NULL,
    .upperBound = NULL,
    .logIncrement = NULL
  ),
  active = list(
    #' @field pathInReferenceSimulation Character. Path to the parameter in the reference simulation.
    pathInReferenceSimulation = function(value) {
      if (missing(value)) {
        return(private$.pathInReferenceSimulation)
      }
      error(
        condition = !is.character(value),
        errorMessage = "Parameter 'pathInReferenceSimulation' must be of type 'character'."
      )
      private$.pathInReferenceSimulation <- value
    },

    #' @field pathInTestSimulation Character. Path to the parameter in the test simulation.
    pathInTestSimulation = function(value) {
      if (missing(value)) {
        return(private$.pathInTestSimulation)
      }
      error(
        condition = !is.character(value),
        errorMessage = "Parameter 'pathInTestSimulation' must be of type 'character'."
      )
      private$.pathInTestSimulation <- value
    },

    #' @field displayName Character. Display name of the parameter.
    displayName = function(value) {
      if (missing(value)) {
        return(private$.displayName)
      }
      error(
        condition = !is.character(value),
        errorMessage = "Parameter 'displayName' must be of type 'character'."
      )
      private$.displayName <- value
    },

    #' @field dimension Character. Dimension of the parameter.
    dimension = function(value) {
      if (missing(value)) {
        return(private$.dimension)
      }
      error(
        condition = !is.character(value),
        errorMessage = "Parameter 'dimension' must be of type 'character'."
      )
      private$.dimension <- value
    },

    #' @field unit Character. Unit of the parameter.
    unit = function(value) {
      if (missing(value)) {
        return(private$.unit)
      }
      error(
        condition = !is.character(value),
        errorMessage = "Parameter 'unit' must be of type 'character'."
      )
      private$.unit <- value
    },

    #' @field lowerBound Numeric. Lower bound for the parameter.
    lowerBound = function(value) {
      if (missing(value)) {
        return(private$.lowerBound)
      }
      error(
        condition = !is.numeric(value),
        errorMessage = "Parameter 'lowerBound' must be of type 'numeric'."
      )
      private$.lowerBound <- value
    },

    #' @field upperBound Numeric. Upper bound for the parameter.
    upperBound = function(value) {
      if (missing(value)) {
        return(private$.upperBound)
      }
      error(
        condition = !is.numeric(value),
        errorMessage = "Parameter 'upperBound' must be of type 'numeric'."
      )
      private$.upperBound <- value
    },

    #' @field logIncrement Logical. Whether the parameter is incremented logarithmically.
    logIncrement = function(value) {
      if (missing(value)) {
        return(private$.logIncrement)
      }
      error(
        condition = !is.logical(value),
        errorMessage = "Parameter 'logIncrement' must be of type 'logical'."
      )
      private$.logIncrement <- value
    }
  )
)

#' Retrieve parameter bounds from a list of VBEParameter objects.
#'
#' @param parameterList List of `VBEParameter` objects.
#'
#' @return A list with two elements: lower bounds and upper bounds.
#' @export
getParameterBounds <- function(parameterList) {
  parameterBounds <- list()
  parameterBounds[[1]] <- sapply(parameterList, function(par) {
    par$lowerBound
  })
  parameterBounds[[2]] <- sapply(parameterList, function(par) {
    par$upperBound
  })
  return(parameterBounds)
}

#' Retrieve parameter paths in the reference simulation from a list of VBEParameter objects.
#'
#' @param parameterList List of `VBEParameter` objects.
#'
#' @return A character vector of parameter paths in the reference simulation.
#' @export
getParameterPathsInReferenceSimulation <- function(parameterList) {
  parameterPaths <- sapply(parameterList, function(par) {
    par$pathInReferenceSimulation
  })
  return(parameterPaths)
}

#' Retrieve parameter paths in the test simulation from a list of VBEParameter objects.
#'
#' @param parameterList List of `VBEParameter` objects.
#'
#' @return A character vector of parameter paths in the test simulation.
#' @export
getParameterPathsInTestSimulation <- function(parameterList) {
  parameterPaths <- sapply(parameterList, function(par) {
    par$pathInTestSimulation
  })
  return(parameterPaths)
}

#' @title vbeParameter
#' @description
#' Initialize a VBEParameter object.
#'
#' @param pathInReferenceSimulation Character. Path to the parameter in the reference simulation.
#' @param pathInTestSimulation Character. Path to the parameter in the test simulation.
#' @param displayName Character. Display name of the parameter (defaults to `pathInReferenceSimulation`).
#' @param dimension Character. Dimension of the parameter.
#' @param unit Character. Unit of the parameter.
#' @param lowerBound Numeric. Lower bound for the parameter (default: 1e-1).
#' @param upperBound Numeric. Upper bound for the parameter (default: 1e1).
#' @param logIncrement Logical. Whether the parameter is incremented log-arithmically (default: TRUE).
#'
#' @return A character vector of parameter paths in the test simulation.
#' @export
vbeParameter <- function(pathInReferenceSimulation,
            pathInTestSimulation,
            displayName = NULL,
            dimension,
            unit,
            lowerBound = 1e-1,
            upperBound = 1e1,
            logIncrement = TRUE){
  VBEParameter$new(
    pathInReferenceSimulation = pathInReferenceSimulation,
    pathInTestSimulation = pathInTestSimulation,
    displayName = displayName,
    dimension = dimension,
    unit = unit,
    lowerBound = lowerBound,
    upperBound = upperBound,
    logIncrement = logIncrement
  )
}
