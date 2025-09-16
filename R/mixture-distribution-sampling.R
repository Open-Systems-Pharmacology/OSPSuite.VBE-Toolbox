getMuList <- function(mclstResults) {
  means <- mclstResults$parameters$mean
  return(lapply(1:ncol(means), function(n) {
    means[, n]
  }))
}

getSigmaList <- function(mclstResults) {
  sigmas <- mclstResults$parameters$variance$sigma
  return(lapply(1:dim(sigmas)[3], function(n) {
    sigmas[, , n]
  }))
}

getSubMus <- function(mclstResults, elementNumbers) {
  muList <- getMuList(mclstResults)
  return(lapply(muList, function(muVec) {
    muVec[elementNumbers]
  }))
}

getSubSigmas <- function(mclstResults, elementNumbers) {
  sigmaList <- getSigmaList(mclstResults)
  return(lapply(sigmaList, function(sigmaMatrix) {
    sigmaMatrix[elementNumbers, elementNumbers]
  }))
}

getProbabilityOfPoint <- function(point, mu, sigma) {
  if (length(point) == 1) {
    return(dnorm(x = point, mean = mu, sd = sigma))
  }

  return(mclust::dmvnorm(matrix(point, ncol = length(mu)), mean = mu, sigma = sigma))
}

getClusterConditionalProbabilities <- function(mclstResults, givenPoints) {
  givenElementNumbers <- which(!is.na(givenPoints))
  proportions <- mclstResults$parameters$pro
  subMus <- getSubMus(mclstResults = mclstResults, elementNumbers = givenElementNumbers)
  subSigmas <- getSubSigmas(mclstResults = mclstResults, elementNumbers = givenElementNumbers)
  pointProbabilityDensitiesInEachCluster <- lapply(1:mclstResults$G, function(n) {
    getProbabilityOfPoint(
      point = givenPoints[givenElementNumbers],
      mu = subMus[[n]],
      sigma = subSigmas[[n]]
    )
  })

  scaledProbabilityOfEachCluster <- sapply(1:mclstResults$G, function(n) {
    pointProbabilityDensitiesInEachCluster[[n]] * proportions[n]
  })

  if (all(scaledProbabilityOfEachCluster == 0)) {
    scaledProbabilityOfEachCluster <- rep(1 / mclstResults$G, mclstResults$G) * proportions
  }

  conditionalProbabilityOfEachClusterGivenPoints <- sapply(1:mclstResults$G, function(n) {
    scaledProbabilityOfEachCluster[n] / sum(scaledProbabilityOfEachCluster)
  })

  return(conditionalProbabilityOfEachClusterGivenPoints)
}


pickCluster <- function(mclstResults, givenPoints, numberOfSamples) {
  conditionalProbabilityOfEachClusterGivenPoints <- getClusterConditionalProbabilities(mclstResults = mclstResults, givenPoints = givenPoints)
  sample(x = 1:mclstResults$G, size = numberOfSamples, replace = TRUE, prob = conditionalProbabilityOfEachClusterGivenPoints)
}


samplesFromConditionalGaussian <- function(mu, sigma, givenPoints, lowerBounds, upperBounds) {
  variableElementNumbers <- which(is.na(givenPoints))
  givenElementNumbers <- which(!is.na(givenPoints))

  a <- givenPoints[givenElementNumbers]

  mu1 <- mu[variableElementNumbers]
  mu2 <- mu[givenElementNumbers]

  S11 <- sigma[variableElementNumbers, variableElementNumbers]
  S12 <- sigma[variableElementNumbers, givenElementNumbers]
  S21 <- sigma[givenElementNumbers, variableElementNumbers]
  S22 <- sigma[givenElementNumbers, givenElementNumbers]

  muConditional <- mu1 + (S12 %*% solve(S22) %*% (a - mu2))
  sigmaConditional <- S11 - (S12 %*% solve(S22) %*% S21)
  numberOfGibbsSamples <- 1000
  sampleVector <- tmvtnorm::rtmvnorm(numberOfGibbsSamples, mean = muConditional[, 1], sigma = sigmaConditional, lower = lowerBounds, upper = upperBounds, algorithm = "gibbs")[numberOfGibbsSamples, ]
  if (NaN %in% sampleVector) {
    sampleVector <- mu1
  }
  return(sampleVector)
}

#' @export
samplesFromMixture <- function(mclstResults, givenPoints, numberOfSamples, lowerBounds, upperBounds) {
  clusterSamples <- pickCluster(mclstResults, givenPoints, numberOfSamples)
  muList <- getMuList(mclstResults)
  sigmaList <- getSigmaList(mclstResults)
  samplesDf <- NULL
  for (cs in clusterSamples) {
    samplesDf <- rbind.data.frame(
      samplesDf,
      samplesFromConditionalGaussian(muList[[cs]], sigmaList[[cs]], givenPoints, lowerBounds, upperBounds)
    )
  }
  return(unname(samplesDf))
}
