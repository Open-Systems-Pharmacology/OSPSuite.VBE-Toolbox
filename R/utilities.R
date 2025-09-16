error <- function(condition, errorMessage = NULL) {
  if (condition) {
    stop(errorMessage)
  }
  return(NULL)
}

squareTheCircle <- function(x) {
  x <- sub(pattern = "[)]", x = x, replacement = "]")
  x <- sub(pattern = "[(]", x = x, replacement = "[")
  return(x)
}

minNull <- function(x) {
  if (is.null(x)) {
    return(NULL)
  }
  return(min(x))
}

maxNull <- function(x) {
  if (is.null(x)) {
    return(NULL)
  }
  return(max(x))
}
