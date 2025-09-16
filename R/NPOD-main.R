#' @title NPODObject
#' @description An NPOD object
#' @export
#' @importFrom R6 R6Class
NPODObject <- R6::R6Class(
  classname = "NPODObject",
  public = list(
    #' @field simulation An `ospsuite` Simulation object.
    simulation = NULL,
    #' @field optimizationParameterList A list of parameters for optimization
    optimizationParameterList = NULL,
    #' @field outputPath The path to the output of the simulation.
    outputPath = NULL,
    #' @field pkData A data.frame of PK parameters
    pkData = NULL,
    #' @field demographicData A data.frame of demographic parameters
    demographicData = NULL,
    #' @field populationData A data.frame of population parameters
    populationData = NULL,
    #' @field npodRunSettings A `NPODRunSettings` object
    npodRunSettings = NULL,
    #' @field cached_mu Numeric value of cached means
    cached_mu = NULL,
    #' @field n_err A data.frame of demographic parameters
    n_err = 0,
    #' @field theta_0 A data.frame of demographic parameters
    theta_0 = NULL,
    #' @field pyl A data.frame of demographic parameters
    pyl = NULL,
    #' @field err_log A data.frame of demographic parameters
    err_log = rep(list(), 5),
    #' @field sigma A data.frame of demographic parameters
    sigma = list(),
    #' @field useLogNormalLikelihood A logical to define if using log normal likelihood
    useLogNormalLikelihood = FALSE,
    
    #' @description Initialize a new NPODObject
    #' @param simulationFilePath A string path to the simulation file
    #' @param optimizationParameterList A list of `VBEParameter` objects
    #' @param outputPath A string path to the output of the simulation
    #' @param pkDataFilePath A string path to the PK data file
    #' @param studyPopulationDataFilePath A string path to the study population data file
    #' @param cofactorPaths A list of cofactor paths
    #' @param initialGridSize An integer for the initial grid size
    #' @param useLogNormalLikelihood A logical to define if using log normal likelihood
    #' @param npodRunSettings A `NPODRunSettings` object
    initialize = function(simulationFilePath,
                          optimizationParameterList,
                          outputPath,
                          pkDataFilePath,
                          studyPopulationDataFilePath = NULL,
                          cofactorPaths = NULL,
                          initialGridSize = NULL,
                          useLogNormalLikelihood = NULL,
                          npodRunSettings = NULL) {
      # if there is a populationData check that the individuals in it are the same as in the pkData
      self$npodRunSettings <- npodRunSettings %||% NPODRunSettings$new()
      self$outputPath <- outputPath
      self$simulation <- ospsuite::loadSimulation(simulationFilePath)
      ospsuite::addOutputs(
        quantitiesOrPaths = self$outputPath,
        simulation = self$simulation
      )

      # check that this is a list of VBEParameter class:
      checkOptimizationParameterList <- all(sapply(optimizationParameterList, inherits, "VBEParameter"))
      if (!checkOptimizationParameterList) {
        stop("Not all optimization parameters are of class VBEParameter.")
      }
      self$optimizationParameterList <- optimizationParameterList

      self$parsePKData(pkDataFilePath)
      self$parsePopulationData(studyPopulationDataFilePath, cofactorPaths)

      self$simulation$outputSchema$addTimePoints(sort(as.numeric(unique(unlist(self$pkData$timeVectorList))))) # add times

      self$cached_mu <- self$multi_mu
      if (!is.null(self$npodRunSettings$npod_cache)) {
        self$cached_mu <- memoise::memoise(self$multi_mu, cache = self$npodRunSettings$npod_cache)
      }

      if (!is.null(useLogNormalLikelihood)) {
        self$useLogNormalLikelihood <- useLogNormalLikelihood
      }

      self$setGrid(initialGridSize)
    },

    #' @description Parse PK data and convert to base unit
    #' @param pkDataFilePath A string path to the PK data file
    parsePKData = function(pkDataFilePath) {
      pkDataDf <- read.csv(pkDataFilePath)
      # Time
      if ("timeUnits" %in% colnames(pkDataDf)) {
        pkDataDf$time <- sapply(
          1:nrow(pkDataDf),
          function(index) {
            ospsuite::toBaseUnit(
              quantityOrDimension = "Time",
              values = pkDataDf$time[index],
              unit = pkDataDf$timeUnits[index]
            )
          }
        )
      }
      pkDataDf$timeUnits <- ospsuite::getBaseUnit("Time")
      # Output
      outputQuantity <- ospsuite::getQuantity(path = self$outputPath, container = self$simulation)
      if ("outputUnits" %in% colnames(pkDataDf)) {
        pkDataDf$outputValues <- sapply(
          1:nrow(pkDataDf), function(index) {
            ospsuite::toBaseUnit(
              quantityOrDimension = outputQuantity$dimension,
              values = pkDataDf$outputValues[index],
              unit = pkDataDf$outputUnits[index]
            )
          }
        )
      }
      pkDataDf$outputUnits <- ospsuite::getBaseUnit(outputQuantity$dimension)

      uniqueIndividualIds <- unique(pkDataDf$id)
      numberOfIndividuals <- length(uniqueIndividualIds)
      timeVectorList <- vector(mode = "list", length = numberOfIndividuals)
      plasmaConcentrationVectorList <- vector(mode = "list", length = numberOfIndividuals)
      idList <- vector(mode = "list", length = numberOfIndividuals)
      min_y <- min(na.omit(pkDataDf$outputValues))

      for (idNo in seq_along(uniqueIndividualIds)) {
        individualId <- uniqueIndividualIds[idNo]
        individualTimeVector <- pkDataDf$time[pkDataDf$id == individualId & pkDataDf$outputValues != 999 & !is.na(pkDataDf$outputValues)]
        individualPlasmaConcentrationProfile <- pkDataDf$outputValues[pkDataDf$id == individualId & pkDataDf$outputValues != 999 & !is.na(pkDataDf$outputValues)]

        if (length(individualPlasmaConcentrationProfile) == 0) {
          next
        }

        timeVectorList[[idNo]] <- individualTimeVector
        plasmaConcentrationVectorList[[idNo]] <- individualPlasmaConcentrationProfile
        idList[[idNo]] <- individualId
      }

      self$pkData <- list(
        timeVectorList = timeVectorList,
        plasmaConcentrationVectorList = plasmaConcentrationVectorList,
        idList = idList
      )
    },

    #' @description parse Population data and convert to base unit
    #' @param studyPopulationDataFilePath A string path to the study population data file
    #' @param cofactorPaths A list of cofactor paths
    parsePopulationData = function(studyPopulationDataFilePath, cofactorPaths = NULL) {
      if (is.null(studyPopulationDataFilePath)) {
        self$populationData <- NULL
        return()
      }

      self$demographicData <- read.csv(studyPopulationDataFilePath)

      if (!is.null(cofactorPaths)) {
        populationDataFrame <- data.frame(IndividualId = self$demographicData$id)
        for (cofactor in names(cofactorPaths)) {
          populationDataFrame[[cofactorPaths[[cofactor]]]] <- self$demographicData[[cofactor]]
        }
        self$populationData <- populationDataFrame
        return()
      }

      individualCharacteristics <- getIndividualCharacteristics(studyPopulationData = self$demographicData)

      individualList <- list()
      populationDataFrame <- NULL
      for (indNo in seq_along(individualCharacteristics)) {
        id <- individualCharacteristics[[indNo]]$id
        individual <- createIndividual(individualCharacteristics[[indNo]]$characteristics)
        individualDf <- cbind(data.frame(t(individual$derivedParameters$values)), data.frame(t(individual$distributedParameters$values)))
        names(individualDf) <- c(individual$derivedParameters$paths, individual$distributedParameters$paths)
        individualDf <- cbind(data.frame(IndividualId = id), individualDf)
        populationDataFrame <- rbind.data.frame(populationDataFrame, individualDf)
      }

      self$populationData <- populationDataFrame
    },

    #' @description setGrid Set Grid for optimization
    #' @param initialGridSize An integer for the initial grid size
    setGrid = function(initialGridSize) {
      numberOfParameters <- length(self$optimizationParameterList)

      size_theta0 <- initialGridSize
      if (is.null(size_theta0)) {
        size_theta0 <- length(self$pkData$idList)
      }

      if (numberOfParameters == 1) {
        uniformMatrix <- t(DiceDesign::runif.faure(size_theta0, 2)$design)
        uniformMatrix <- matrix(uniformMatrix[1, ], ncol = length(uniformMatrix[1, ]))
      } else {
        uniformMatrix <- t(DiceDesign::runif.faure(size_theta0, numberOfParameters)$design)
      }

      parameterMatrix <- uniformMatrix

      for (parNo in seq_along(self$optimizationParameterList)) {
        a <- self$optimizationParameterList[[parNo]]$lowerBound
        b <- self$optimizationParameterList[[parNo]]$upperBound
        logInc <- self$optimizationParameterList[[parNo]]$logIncrement

        if (logInc) {
          log_a <- logSafe(a)
          log_b <- logSafe(b)
          parameterMatrix[parNo, ] <- exp(log_a + uniformMatrix[parNo, ] * (log_b - log_a))
        } else {
          parameterMatrix[parNo, ] <- a + uniformMatrix[parNo, ] * (b - a)
        }
      }

      self$theta_0 <- parameterMatrix
    },

    #' @description Run NPOD
    runNPOD = function() {
      npodResults <- self$Dopt()
      return(npodResults)
    },

    #' @description Dopt Optimization algorithm
    Dopt = function() {
      # browser()
      old_theta <- self$theta_0
      counter <- 1
      F0 <- -10^(30)
      F1 <- 2 * F0
      lam <- rep(1, ncol(old_theta)) / ncol(old_theta)
      P1 <- self$getPSI(theta = old_theta, optimizeSigma = TRUE, lambda = lam)
      options <- neldermead::optimset(MaxFunEvals = 2000000000, TolX = 1e-14, MaxIter = self$npodRunSettings$MaxIter, TolFun = 1e-14)
      objective_function_values <- c() # Initialize a vector to store the objective function values at each iteration
      objective_function_values[counter] <- F0
      objective_function_values[counter + 1] <- F1
      
      cli::cli_h2("Optimization")
      while (abs(objective_function_values[counter + 1] - objective_function_values[counter]) > self$npodRunSettings$theta_F) {
        cli::cli_alert_info("Run {counter}:")
        cli::cli_alert("theta:")
        print(old_theta)
        cli::cli_alert("lambda: {lam}")
        IPM_results_old_theta <- IPM(P1)
        lam <- IPM_results_old_theta$lambda # Dual variable obtained from IPM optimization with theta = old_theta

        # CONDENSE - remove low probability support points
        cli::cli_h3("Condense - removing low probability support points")
        L1 <- 0.00000001
        L2 <- max(lam) / 1000
        ind1 <- (lam > L1) & (lam > L2) # Find elements of lambda that are greater than the lower bounds L1 and L2
        lam <- lam[ind1]
        inb_theta <- old_theta[, ind1, drop = FALSE] # Select the columns of theta corresponding to values of lambda that fall within the constraints L_lower and L_upper
        cli::cli_alert("theta:")
        print(inb_theta)
        cli::cli_alert("lambda: {lam}")
        
        # Calculate the matrix of log likelihoods P2 for the set of support points inb_theta
        P2 <- self$getPSI(theta = inb_theta, optimizeSigma = TRUE, lambda = lam)

        IPM_results_inb_theta <- IPM(P2)
        lam <- IPM_results_inb_theta$lambda # Dual variable obtain from IPM optimization with theta = inb_theta
        objective_function_values[counter + 2] <- IPM_results_inb_theta$fobj # Add new objective value to objective_function_values

        # CONDENSE - remove low probability support points
        cli::cli_h3("Condense - removing low probability support points")
        L3 <- max(lam) / 1000
        ind2 <- lam > L3 # Find elements of lambda that are greater than the lower bound L3
        lam <- lam[ind2]
        new_weights <- lam / sum(lam)
        new_theta <- inb_theta[, ind2, drop = FALSE]
        cli::cli_alert("theta:")
        print(new_theta)
        cli::cli_alert("lambda: {lam}")

        if (abs(objective_function_values[counter + 2] - objective_function_values[counter + 1]) <= self$npodRunSettings$theta_F) {
          cli::cli_alert_success("{.strong Stopping}: Absolute change in objective function is smaller than {.strong theta_F} ({self$npodRunSettings$theta_F})")
          cli::cli_alert_info("theta:")
          print(new_theta)
          break
        }

        if (counter >= self$npodRunSettings$ncycles) {
          cli::cli_alert_warning("{.strong Stopping}: Maximum number of cycles reached ({counter})")
          cli::cli_alert_info("theta:")
          print(new_theta)
          break
        }

        K <- length(new_theta[1, ])
        self$pyl <- P2[, ind2, drop = FALSE] %*% new_weights
        cli::cli_h3("fminsearch and prune")
        for (l in 1:K) {
          cli::cli_alert_info("Step {l}/{K}:")
          cli::cli_process_start("fminsearch")
          cand_theta <- tryCatch(
            {
              neldermead::fminsearch(self$multi_D, new_theta[, l], options)
            },
            error = function(e) {
              NULL
            }
          )
          cli::cli_process_done()
          cli::cli_process_start("prune")
          new_theta <- self$prune(theta = new_theta, theta_plus = cand_theta$optbase$xopt)
          cli::cli_process_done()
          cli::cli_alert_info("theta:")
          print(new_theta)
        }
        old_theta <- new_theta
        counter <- counter + 1
        # Calculate the matrix of log likelihoods P1 for the set of support points old_theta
        P1 <- self$getPSI(theta = old_theta)
      }

      return(list("count" = counter, "theta" = new_theta, "w" = new_weights, "LogLikelihood" = objective_function_values[length(objective_function_values)], "PSI" = P2, demographicData = self$demographicData))
    },

    #' @description The `prune` function checks if a candidate support point (`theta_plus`)
    #' is sufficiently different from the existing set of support points (`theta`) and within
    #' specified bounds (`a`, `b`). If these conditions are met, `theta_plus` is
    #' added as a new column to `theta`.
    #' @param theta A matrix of support points
    #' @param theta_plus A candidate support point
    prune = function(theta, theta_plus) {
      if (is.null(theta_plus)) {
        return(theta)
      }
      # Ensure theta_plus is treated as a column vector
      theta_plus <- as.numeric(theta_plus)

      # Determine lower and upper bounds with log scaling where applicable
      a <- sapply(self$optimizationParameterList, function(par) {
        ifelse(test = isTRUE(par$logIncrement), log10Safe(par$lowerBound), par$lowerBound)
      })
      b <- sapply(self$optimizationParameterList, function(par) {
        ifelse(test = isTRUE(par$logIncrement), log10Safe(par$upperBound), par$upperBound)
      })

      # Scale theta_plus based on logIncrement flag
      theta_plus_scaled <- sapply(seq_along(theta_plus), function(n) {
        ifelse(test = isTRUE(self$optimizationParameterList[[n]]$logIncrement), log10Safe(theta_plus[n]), theta_plus[n])
      })

      # Scale theta matrix based on logIncrement flag
      theta_scaled <- theta
      for (n in seq_len(nrow(theta_scaled))) {
        if (self$optimizationParameterList[[n]]$logIncrement) {
          theta_scaled[n, ] <- log10Safe(theta_scaled[n, ])
        }
      }

      if (ncol(theta_scaled) > 0) {
        # Broadcast theta_plus_scaled to the same dimensions as theta_scaled
        distances <- colSums(abs(theta_scaled - matrix(theta_plus_scaled, nrow = nrow(theta_scaled), ncol = ncol(theta_scaled), byrow = FALSE)) / (b - a))
        min_dist <- min(distances)
      } else {
        min_dist <- Inf
      }

      # Boundary checks
      within_bounds <- all(theta_plus_scaled > a & theta_plus_scaled < b)
      # Add theta_plus_scaled if conditions are met
      if (min_dist > self$npodRunSettings$theta_d && within_bounds) {
        cli::cli_alert_success("Adding new support point: ({theta_plus})")
        theta <- cbind(theta, theta_plus)
        return(theta)
      }

      # Otherwise reject candidate support point and provide reason.
      cli::cli_alert_danger("Rejecting candidate support point: ({theta_plus})")
      if (!within_bounds) {
        cli::cli_alert("{.strong Reason}: Outside of bounds.")
      } else {
        cli::cli_alert("{.strong Reason}: Too close to existing support point.")
      }
      return(theta)
    },

    #' @description multi_D computes a likelihood-based metric, `D_comp`, by summing the
    #' probabilities from `getPSI()`, normalized by `self$pyl`.
    #' It starts with a penalty term of `-N`, calls `getPSI()`, and updates `D_comp`.
    #' @param theta_parameter A matrix of support points
    multi_D = function(theta_parameter) {
      N <- length(self$pkData$plasmaConcentrationVectorList) # nsub
      D_comp <- -N
      D_comp <- D_comp + sum(self$getPSI(theta = theta_parameter) / self$pyl)
      return(D_comp)
    },

    #' @description getPIS The `getPSI` function calculates the probability of the observed data
    #' `mprob` is an N x K matrix where each entry `[i, l]` represents the probability
    #' of observing `y[[i]]` given the model prediction `ySimList[[i, l]]`.
    #' A small lower bound, `1e-100`, is enforced for numerical stability.
    #' @param theta A matrix of support points
    #' @param ySimList A list of model predictions
    #' @param useLog A logical to define if using log normal likelihood
    #' @param optimizeSigma A logical to define if optimizing sigma
    #' @param lambda A list of lambda values
    getPSI = function(theta, ySimList = NULL, useLog = self$useLogNormalLikelihood, optimizeSigma = FALSE, lambda = NULL) {
      y <- self$pkData$plasmaConcentrationVectorList
      t <- self$pkData$timeVectorList

      # Determine K efficiently
      K <- if (is.matrix(theta)) ncol(theta) else length(theta)

      # Compute ySimList if not provided
      if (is.null(ySimList)) {
        t1 <- system.time({
          ySimList <- self$cached_mu(theta)
        })
      }

      if (optimizeSigma) {
        # Optimize sigma
        for (ind in seq_len(nrow(ySimList))) {
          if (useLog) {
            ySimWeighted <- rowSums(sapply(seq_along(lambda), function(n) {
              lambda[n] * logSafe(ySimList[[ind, n]])
            }))
            sigma <- sqrt(sum((logSafe(ySimWeighted) - logSafe(y[[ind]]))^2))
          } else {
            ySimWeighted <- rowSums(sapply(seq_along(lambda), function(n) {
              lambda[n] * ySimList[[ind, n]]
            }))
            sigma <- sqrt(sum((ySimWeighted - y[[ind]])^2))
          }
          self$sigma[[ind]] <- rep(sigma, length(y[[ind]]))
        }
      }

      # Vectorized probability computation
      t1 <- system.time({
        mprob <- sapply(1:K, function(l) {
          sapply(1:length(y), function(i) {
            max(1e-100, self$getLikelihood(y = y[[i]], sigma = self$sigma[[i]], ySim = ySimList[[i, l]], useLog = useLog))
          })
        })
      })
      return(mprob)
    },

    #' @description Calculates the likelihood of the observed data, `y`,
    #' given the model predictions, `ySim`, and measurement error, `sigma`.
    #' The default case uses arithmetic error, but if `use_log = TRUE`,
    #' it calculates the likelihood using the logarithmic (ratio) error.
    #' @param y A vector of observed data
    #' @param t A vector of time points
    #' @param theta A matrix of support points
    #' @param sigma A vector of measurement errors
    #' @param ySim A vector of model predictions
    #' @param useLog A logical to define if using log normal likelihood
    getLikelihood = function(y, t, theta, sigma, ySim, useLog = FALSE) {
      j <- length(y)
      z <- matrix(0, 1, j)

      # Ensure sizes match
      if (length(y) != length(ySim)) {
        warning("Observed and predicted sizes don't match")
        stop()
      }

      if (useLog) {
        # Log-transformed error (ratio error)
        z <- (logSafe(y) - logSafe(ySim))^2
        logProb <- -logSafe((sqrt(2 * pi) * sigma * y)) - (z / (2 * (sigma^2)))
        # Calculate the final probability by exponentiating the logSafe likelihood and summing over all data points
        prob <- exp(sum(logProb[sigma != 0 & y != 0]))
      } else {
        # Arithmetic error (default case)
        z <- (y - ySim)^2
        logProb <- -logSafe((sqrt(2 * pi) * sigma)) - (z / (2 * (sigma^2)))
        # Calculate the final probability by exponentiating the logSafe likelihood and summing over all data points
        prob <- exp(sum(logProb[sigma != 0]))
      }
      return(prob)
    },

    #' @description computes the model predictions for each individual in the study population
    #' given a set of support points, `theta`.
    #' @param theta A matrix of support points
    multi_mu = function(theta) {
      timePointsList <- self$pkData$timeVectorList

      n_ind <- length(timePointsList)

      inferenceParameterPaths <- sapply(self$optimizationParameterList, function(x) {
        x$pathInReferenceSimulation
      })

      popDf <- data.frame(matrix(ncol = 0, nrow = n_ind))

      if (!is.null(self$populationData)) {
        populationParameterPaths <- setdiff(colnames(self$populationData), c("IndividualId", inferenceParameterPaths))
        for (path in populationParameterPaths) {
          popDf[[path]] <- self$populationData[[path]]
        }
      }

      for (path in inferenceParameterPaths) {
        popDf[[path]] <- 0
      }

      ySimList <- matrix(rep(list(), n_ind * length(theta[1, ])), nrow = n_ind, ncol = length(theta[1, ]))

      bigPopDf <- NULL
      for (sp in seq_len(ncol(theta))) {
        df <- popDf
        for (pathNo in seq_along(inferenceParameterPaths)) {
          df[[inferenceParameterPaths[pathNo]]] <- theta[pathNo, sp]
        }
        bigPopDf <- rbind.data.frame(bigPopDf, df)
      }
      bigPopDf <- cbind(data.frame(IndividualId = seq(nrow(bigPopDf))), bigPopDf)

      write.csv(x = bigPopDf, file = "population.csv", row.names = FALSE)

      population <- ospsuite::loadPopulation("population.csv")
      file.remove("population.csv")

      res <- ospsuite::runSimulation(simulation = self$simulation, population = population)
      resData <- ospsuite::getOutputValues(res, quantitiesOrPaths = self$outputPath)$data

      for (ind in unique(resData$IndividualId)) {
        df <- resData[resData$IndividualId == ind, ]
        individualNumber <- ((ind - 1) %% n_ind) + 1
        ptNumber <- ((ind - 1) %/% n_ind) + 1
        ySimList[[individualNumber, ptNumber]] <- df[[self$outputPath]][df$Time %in% timePointsList[[individualNumber]]]
      }

      return(ySimList)
    }
  )
)


#' @title NPODRunSettings
#' @description An NPODRunSettings object
#' @export
NPODRunSettings <- R6::R6Class(
  classname = "NPODRunSettings",
  public = list(
    #' @field noise A logical to define if noise is present
    noise = NULL,
    #' @field model A model object
    model = NULL,
    #' @field ncycles The number of cycles
    ncycles = NULL,
    #' @field MaxIter The maximum number of iterations
    MaxIter = NULL,
    #' @field theta_F The theta_F value
    theta_F = NULL,
    #' @field theta_d The theta_d value
    theta_d = NULL,
    #' @field size_theta0 The size of theta0
    size_theta0 = NULL,
    #' @field cache_folder_name The name of the cache folder
    cache_folder_name = NULL,
    #' @field npod_cache The cache object
    npod_cache = NULL,
    
    #' @description Initialize a new NPODRunSettings
    #' @param noise A logical to define if noise is present
    #' @param model A model object
    #' @param ncycles The number of cycles
    #' @param MaxIter The maximum number of iterations
    #' @param theta_F The theta_F value
    #' @param theta_d The theta_d value
    #' @param size_theta0 The size of theta0
    #' @param cache_folder_name The name of the cache folder
    initialize = function(noise = FALSE,
                          model = NULL,
                          ncycles = Inf,
                          MaxIter = 10,
                          theta_F = 10e-2,
                          theta_d = 10e-4,
                          size_theta0 = NULL,
                          cache_folder_name = NULL) {
      self$noise <- noise
      self$model <- model
      self$ncycles <- ncycles
      self$MaxIter <- MaxIter
      self$theta_F <- theta_F
      self$theta_d <- theta_d
      self$size_theta0 <- size_theta0

      self$npod_cache <- cachem::cache_mem()
      if (!is.null(cache_folder_name)) {
        self$cache_folder_name <- cache_folder_name
        self$npod_cache <- cachem::cache_disk(cache_folder_name)
      }
    }
  )
)

#' @title getIndividualCharacteristics
#' @description A function to create a list of individuals with `ospsuite` characteristics
#' @param studyPopulationData A data.frame
#' @return A list of individuals with characteristics
#' @importFrom ospsuite createIndividualCharacteristics
getIndividualCharacteristics <- function(studyPopulationData) {
  individualCharacteristics <- list()
  for (indNo in seq_along(studyPopulationData[[standardColumnNames$idColumn]])) {
    individualCharacteristics[[indNo]] <- list()
    id <- studyPopulationData[[standardColumnNames$idColumn]][indNo]
    individualCharacteristics[[indNo]]$id <- id
    individualCharacteristics[[indNo]]$characteristics <- ospsuite::createIndividualCharacteristics(
      species = studyPopulationData[[standardColumnNames$speciesColumn]][studyPopulationData$id == id],
      population = studyPopulationData[[standardColumnNames$populationColumn]][studyPopulationData$id == id],
      gender = studyPopulationData[[standardColumnNames$genderColumn]][studyPopulationData$id == id],
      weight = studyPopulationData[[standardColumnNames$weightColumn]][studyPopulationData$id == id],
      weightUnit = studyPopulationData[[standardColumnNames$weightUnitsColumn]][studyPopulationData$id == id] %||% ospsuite::ospUnits$Mass$kg,
      height = studyPopulationData[[standardColumnNames$heightColumn]][studyPopulationData$id == id],
      heightUnit = studyPopulationData[[standardColumnNames$heightUnitsColumn]][studyPopulationData$id == id] %||% ospsuite::ospUnits$Length$cm,
      age = studyPopulationData[[standardColumnNames$ageColumn]][studyPopulationData$id == id],
      ageUnit = studyPopulationData[[standardColumnNames$ageUnitsColumn]][studyPopulationData$id == id] %||% ospsuite::ospUnits$`Age in years`$`year(s)`,
      gestationalAge = studyPopulationData[[standardColumnNames$gestationalAgeColumn]][studyPopulationData$id == id] %||% 40,
      gestationalAgeUnit = studyPopulationData[[standardColumnNames$gestationalAgeUnitsColumn]][studyPopulationData$id == id] %||% ospsuite::ospUnits$`Age in weeks`$`week(s)`
    )
  }
  return(individualCharacteristics)
}
