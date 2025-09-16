#' Safe Logarithm Function
#'
#' Computes the logarithm of a numeric input, replacing non-positive values with a small positive threshold to avoid errors.
#'
#' @param x Numeric vector. The values for which the logarithm is computed.
#' @param base Numeric. The base of the logarithm. Default is the natural logarithm (`exp(1)`).
#' @param min_value Numeric. A small positive threshold to replace non-positive values. Default is `.Machine$double.xmin`.
#'
#' @return A numeric vector of logarithm values, where non-positive inputs are replaced with `min_value` before computing the log.
#'
#' @examples
#' logSafe(c(1, 10, 0, -5))
#' logSafe(c(0.1, 0.5, 0), base = 2)
#'
#' @export
logSafe <- function(x, base = exp(1), min_value = .Machine$double.xmin) {
  # Ensure x is numeric
  if (!is.numeric(x)) stop("Input must be numeric.")

  # Replace non-positive values with a small positive threshold
  x_safe <- ifelse(x > 0, x, min_value)

  # Compute log with the specified base
  log_result <- log(x_safe, base)

  # Handle NaN cases explicitly (if min_value is too small)
  log_result[is.nan(log_result)] <- -Inf

  return(log_result)
}

#' Safe Base-10 Logarithm Function
#'
#' Computes the base-10 logarithm of a numeric input using `logSafe()`.
#'
#' @param x Numeric vector. The values for which the base-10 logarithm is computed.
#'
#' @return A numeric vector of base-10 logarithm values, where non-positive inputs are replaced with a small positive threshold.
#'
#' @examples
#' log10Safe(c(1, 10, 0, -5))
#'
#' @export
log10Safe <- function(x) {
  logSafe(x, base = 10)
}


#' Compute Posterior Probabilities
#'
#' Computes posterior probabilities given estimated support points and weights.
#'
#' @param ans A list containing `PSI`, a matrix of support points (rows: subjects, columns: support points), and `w`, a vector of weights.
#'
#' @return A matrix of posterior probabilities for each subject and support point.
#'
#' @examples
#' ans <- list(PSI = matrix(runif(9), 3, 3), w = c(0.2, 0.5, 0.3))
#' posterior(ans)
#'
#' @export
posterior <- function(ans) {
  # psi[j,i] jth subject, ith support point
  psi <- ans$PSI
  nsub <- nrow(psi)
  nspp <- ncol(psi)
  w <- matrix(ans$w, nrow = nspp)
  post <- matrix(0, nsub, nspp)
  py <- psi %*% w
  for (j in 1:nsub) {
    for (i in 1:nspp) {
      post[j, i] <- psi[j, i] * w[i] / py[j]
    }
  }
  return(post)
}

#' Compute Weighted Support Points
#'
#' Computes weighted support points using posterior probabilities and demographic parameters.
#'
#' @param res A list containing `theta` (a matrix of support points), `w` (weights), and `demographicData` (a data frame of demographic information).
#' @param params A character vector of parameter names to extract from `demographicData`. Default is `c('Age', 'Weight', 'Height')`.
#'
#' @return A list with `points`, a data frame of weighted support points, and `weights`, a vector of posterior probabilities.
#'
#' @examples
#' res <- list(
#' theta = matrix(c(0.25,0.625,0.125,0.9375),ncol = 2),
#' w = c(0.7, 0.3),
#' PSI = matrix(runif(6),ncol=2),
#' demographicData = data.frame(Age = c(25, 30, 35),
#' Weight = c(70, 80, 90), Height = c(170, 175, 180))
#' )
#' getWeightedPoints(res)
#'
#' @export
getWeightedPoints <- function(res, params = c("Age", "Weight", "Height")) {
  posteriorProbabilities <- posterior(res)
  numberTheta <- nrow(res$theta)
  numberSupportPoints <- ncol(res$theta)
  parameterNames <- c(paste0("theta", seq(numberTheta)), params)
  columnNames <- c("posteriorProbability", parameterNames)
  numberOfIndividuals <- nrow(res$demographicData)
  weightedPointsDf <- data.frame(matrix(ncol = length(columnNames), nrow = 0))
  for (thetaIndex in seq_len(numberSupportPoints)) {
    thetaValues <- res$theta[, thetaIndex]
    for (ind in seq_len(numberOfIndividuals)) {
      posteriorProbability <- posteriorProbabilities[ind, thetaIndex]
      demParams <- sapply(params, function(par) {
        res$demographicData[ind, ][[par]]
      })
      rowValues <- c(posteriorProbability, thetaValues, demParams)
      weightedPointsDf <- rbind.data.frame(weightedPointsDf, rowValues)
    }
    colnames(weightedPointsDf) <- columnNames
  }
  return(list(points = weightedPointsDf[, parameterNames], weights = weightedPointsDf$posteriorProbability))
}
