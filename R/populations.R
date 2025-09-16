unifyUnits <- function(values, units, dimension) {
  if (is.null(values)) {
    return(list(values = NULL, units = NULL))
  }

  allUnits <- unique(units)
  if (length(allUnits) == 1) {
    return(list(values = values, units = allUnits))
  }

  values <- sapply(seq_along(values), function(n) {
    ospsuite::toBaseUnit(
      quantityOrDimension = dimension,
      values = values[n],
      unit = units[n]
    )
  })
  units <- ospsuite::getBaseUnit(dimension)
  return(list(values = values, units = units))
}

unifyPopulationDataUnits <- function(demographicData) {
  cofactorUnits <- NULL

  weightsDimension <- cofactorDimensions[[standardColumnNames$weightColumn]]
  weights <- unifyUnits(
    values = demographicData[[standardColumnNames$weightColumn]],
    units = demographicData[[standardColumnNames$weightUnitsColumn]],
    dimension = weightsDimension
  )
  demographicData[[standardColumnNames$weightColumn]] <- weights$values
  demographicData[[standardColumnNames$weightUnitsColumn]] <- weights$units
  cofactorUnits[[standardColumnNames$weightColumn]] <- weights$units

  heightsDimension <- cofactorDimensions[[standardColumnNames$heightColumn]]
  heights <- unifyUnits(
    values = demographicData[[standardColumnNames$heightColumn]],
    units = demographicData[[standardColumnNames$heightUnitsColumn]],
    dimension = heightsDimension
  )
  demographicData[[standardColumnNames$heightColumn]] <- heights$values
  demographicData[[standardColumnNames$heightUnitsColumn]] <- heights$units
  cofactorUnits[[standardColumnNames$heightColumn]] <- heights$units

  agesDimension <- cofactorDimensions[[standardColumnNames$ageColumn]]
  ages <- unifyUnits(
    values = demographicData[[standardColumnNames$ageColumn]],
    units = demographicData[[standardColumnNames$ageUnitsColumn]],
    dimension = agesDimension
  )
  demographicData[[standardColumnNames$ageColumn]] <- ages$values
  demographicData[[standardColumnNames$ageUnitsColumn]] <- ages$units
  cofactorUnits[[standardColumnNames$ageColumn]] <- ages$units

  gestationalAgesDimension <- cofactorDimensions[[standardColumnNames$gestationalAgeColumn]]
  gestationalAges <- unifyUnits(
    values = demographicData[[standardColumnNames$gestationalAgeColumn]],
    units = demographicData[[standardColumnNames$gestationalAgeUnitsColumn]],
    dimension = gestationalAgesDimension
  )
  demographicData[[standardColumnNames$gestationalAgeColumn]] <- gestationalAges$values
  demographicData[[standardColumnNames$gestationalAgeUnitsColumn]] <- gestationalAges$units
  cofactorUnits[[standardColumnNames$gestationalAgeColumn]] <- gestationalAges$units

  return(list(demographicData = demographicData, cofactorUnits = cofactorUnits))
}

getClusteringNPOD <- function(inferredDistribution, parameterNames, cofactorNames, numberOfClusters) {
  cli::cli_progress_step("Building cluster model based on NPOD results.")
  weighted_points <- getWeightedPoints(
    res = inferredDistribution,
    params = cofactorNames
  )
  sampledPoints <- weighted_points$points[sample(1:length(weighted_points$weights),
                                                 size = 100,
                                                 replace = TRUE,
                                                 prob = weighted_points$weights
  ), ]
  mclustBIC <- mclust::mclustBIC
  clusters <- mclust::Mclust(data = sampledPoints, G = numberOfClusters)
  plot(clusters, what = "classification")
  cli::cli_progress_done()
  return(clusters)
}

getClustersFunctions <- list(NPOD = getClusteringNPOD)

getParameterListNPOD <- function(inferredDistribution) {
  return(inferredDistribution$parameters)
}

getParameterListFunctions <- list(NPOD = getParameterListNPOD)


#' @title getClusters
#' @description Get clusters.
#' @param inferredDistribution The inferred distribution.
#' @param numberOfClusters The number of clusters.
#' @export
getClusters <- function(inferredDistribution, numberOfClusters) {
  parameterList <- getParameterListFunctions[[inferredDistribution$method]](inferredDistribution)
  parameterNames <- sapply(parameterList, function(par) {
    par$displayName
  })
  cofactorNames <- intersect(allCofactorNames, names(inferredDistribution$demographicData))
  clusters <- getClustersFunctions[[inferredDistribution$method]](inferredDistribution, parameterNames, cofactorNames, numberOfClusters)
  return(clusters)
}

#' VirtualPopulation R6 Class
#'
#' A class for generating and simulating virtual populations with support for inter-occasion variability (IOV).
#'
#' @description
#' The `VirtualPopulation` class provides functionality to:
#' - Create virtual populations based on an inferred multivariate distribution and demographic constraints.
#' - Add additive or multiplicative inter-occasion variability (IOV) to model parameters.
#' - Simulate model outputs for a virtual population using reference and test simulation files.
#'
#' @section Public Methods:
#' - `initialize(inferredDistribution, demographyRanges, numberOfVirtualIndividuals, proportionOfFemales, numberOfClusters)`: Constructor for the class.
#' - `addAdditiveIOV(parameterPath, SD)`: Add additive IOV to a specified parameter.
#' - `addMultiplicativeIOV(parameterPath, CV)`: Add multiplicative IOV to a specified parameter.
#' - `simulateVirtualPopulation(referenceSimulationFilePath, testSimulationFilePath, numberOfReplicates, outputPath, startTime, endTime, resolutionPtsMin)`: Simulate the virtual population using specified input files and time settings.
#'
#' @section Active Bindings:
#' - `inferredDistribution`: Returns the distribution used to generate the virtual population.
#' - `referencePopulationDataframe`: Returns the reference population data frame.
#' - `testPopulationDataframe`: Returns the test population data frame.
#' - `iovListSD`: Returns the list of additive IOV values.
#' - `iovListCV`: Returns the list of multiplicative IOV values.
#'
#' @return An object of class `VirtualPopulation`.
#'
#' @examples
#' \dontrun{
#' vpop <- VirtualPopulation$new(
#'   inferredDistribution = inferredDistribution,
#'   demographyRanges = list(age = c(18, 65)),
#'   numberOfVirtualIndividuals = 200,
#'   proportionOfFemales = 60,
#'   numberOfClusters = 5
#' )
#' vpop$addAdditiveIOV("PATH|TO|PARAMETER1", SD = 0.1)
#' vpop$addMultiplicativeIOV("PATH|TO|PARAMETER2", CV = 0.2)
#' vpop$simulateVirtualPopulation(
#'   referenceSimulationFilePath = "path/to/ref.sim",
#'   testSimulationFilePath = "path/to/test.sim",
#'   outputPath = "sim_outputs",
#'   startTime = 0,
#'   endTime = 24,
#'   resolutionPtsMin = 10
#' )
#' }
#'
#' @export
VirtualPopulation <- R6::R6Class(classname = "VirtualPopulation",
                                 public = list(
                                   #' @description
                                   #' Initialize a new `VirtualPopulation` object.
                                   #'
                                   #' @param inferredDistribution A distribution object for generating virtual individuals.
                                   #' @param demographyRanges A list specifying demographic constraints (e.g., age range, weight limits).
                                   #' @param numberOfVirtualIndividuals Integer. Number of virtual individuals to simulate. Default is 100.
                                   #' @param proportionOfFemales Numeric. Proportion (%) of females in the virtual population. Default is 50.
                                   #' @param numberOfClusters Integer. Number of clusters to use for stratifying the distribution. Default is 4.
                                   #' @return A new `VirtualPopulation` object.
                                   initialize = function(inferredDistribution,
                                                         demographyRanges,
                                                         numberOfVirtualIndividuals = 100,
                                                         proportionOfFemales = 50,
                                                         numberOfClusters = 4){
                                     private$.inferredDistribution <- inferredDistribution
                                     private$.createVirtualPopulation(demographyRanges = demographyRanges,
                                                                      numberOfVirtualIndividuals = numberOfVirtualIndividuals,
                                                                      proportionOfFemales = proportionOfFemales,
                                                                      numberOfClusters = numberOfClusters)

                                   },
                                   #' @description
                                   #' Add additive inter-occasion variability (IOV) to a specific parameter in the virtual population.
                                   #'
                                   #' @param parameterPath Character. The path to the parameter in the model to which additive IOV is applied.
                                   #' @param SD Numeric. Standard deviation of the additive IOV. Must be a single positive number.
                                   #' @return None. Adds to internal IOV list.
                                   #' @examples
                                   #' \dontrun{
                                   #' vpop$addAdditiveIOV("PATH|TO|PARAMETER1", SD = 0.1)
                                   #' }
                                   addAdditiveIOV = function(parameterPath,SD){
                                     if (!is.numeric(SD) || length(SD) != 1 || SD <= 0){
                                       stop("SD must be a single positive numeric value.")
                                     }
                                     private$.checkIOVParameterReplacement(parameterPath = parameterPath,variabilityType = "additive")
                                     private$.iovListSD[[parameterPath]] <- SD
                                   },
                                   #' @description
                                   #' Add multiplicative inter-occasion variability (IOV) to a specific parameter.
                                   #' @param parameterPath Character. The path to the parameter in the model.
                                   #' @param CV Numeric. Coefficient of variation (CV) for the multiplicative IOV. Must be a single positive number.
                                   #'
                                   #' @return None. Adds to internal IOV list.
                                   #'
                                   #' @examples
                                   #' \dontrun{
                                   #' vpop$addMultiplicativeIOV("PATH|TO|PARAMETER2", CV = 0.2)
                                   #' }

                                   addMultiplicativeIOV = function(parameterPath,CV){
                                     if (!is.numeric(CV) || length(CV) != 1 || CV <= 0){
                                       stop("CV must be a single positive numeric value.")
                                     }
                                     private$.checkIOVParameterReplacement(parameterPath = parameterPath,variabilityType = "multiplicative")
                                     private$.iovListCV[[parameterPath]] <- CV
                                   },

                                   #' @description
                                   #' Run simulations on the virtual population using provided simulation templates.
                                   #' @param referenceSimulationFilePath Character. File path to the reference model simulation template.
                                   #' @param testSimulationFilePath Character. File path to the test model simulation template.
                                   #' @param numberOfReplicates Integer. Number of replicate simulations per virtual subject. Default is 1.
                                   #' @param outputPath Character. Directory to store simulation results.
                                   #' @param startTime Numeric. Start time of the simulation window.
                                   #' @param endTime Numeric. End time of the simulation window.
                                   #' @param resolutionPtsMin Numeric. Output resolution in points per minute.
                                   #'
                                   #' @return A list containing results from reference and test population simulations.
                                   #'
                                   #' @examples
                                   #' \dontrun{
                                   #' vpop$simulateVirtualPopulation(
                                   #'   referenceSimulationFilePath = "ref.sim",
                                   #'   testSimulationFilePath = "test.sim",
                                   #'   outputPath = "PATH|TO|OUTPUT",
                                   #'   startTime = 0,
                                   #'   endTime = 1440,
                                   #'   resolutionPtsMin = 10
                                   #' )
                                   #' }
                                   simulateVirtualPopulation = function(referenceSimulationFilePath,
                                                                        testSimulationFilePath,
                                                                        numberOfReplicates = 1,
                                                                        outputPath,
                                                                        startTime,
                                                                        endTime,
                                                                        resolutionPtsMin){

                                     refAndTestSimulationsInVirtualPopulation<- simulateVirtualPopulation(referenceSimulationFilePath = referenceSimulationFilePath,
                                                                                                          testSimulationFilePath = testSimulationFilePath,
                                                                                                          referencePopulationDataframe = private$.referencePopulationDataframe,
                                                                                                          testPopulationDataframe = private$.testPopulationDataframe,
                                                                                                          numberOfReplicates = numberOfReplicates,
                                                                                                          outputPath = outputPath,
                                                                                                          startTime = startTime,
                                                                                                          endTime = endTime,
                                                                                                          resolutionPtsMin = resolutionPtsMin,
                                                                                                          iovListSD = private$.iovListSD,
                                                                                                          iovListCV = private$.iovListCV)

                                     return(refAndTestSimulationsInVirtualPopulation)

                                   }

                                 ),
                                 active = list(
                                   #' @field inferredDistribution Get the distribution object used to generate the virtual population.
                                   inferredDistribution = function(value) {
                                     if (missing(value)) {
                                       return(private$.inferredDistribution)
                                     }
                                   },

                                   #' @field referencePopulationDataframe Get the data frame of the reference virtual population.
                                   referencePopulationDataframe = function(value) {
                                     if (missing(value)) {
                                       return(private$.referencePopulationDataframe)
                                     }
                                   },

                                   #' @field testPopulationDataframe Get the data frame of the test virtual population.
                                   testPopulationDataframe = function(value) {
                                     if (missing(value)) {
                                       return(private$.testPopulationDataframe)
                                     }
                                   },

                                   #' @field iovListSD Get the list of additive IOV standard deviations by parameter.
                                   iovListSD = function(value) {
                                     if (missing(value)) {
                                       return(private$.iovListSD)
                                     }
                                   },

                                   #' @field iovListCV Get the list of multiplicative IOV coefficients of variation by parameter.
                                   iovListCV = function(value) {
                                     if (missing(value)) {
                                       return(private$.iovListCV)
                                     }
                                   }

                                 ),

                                 private = list(

                                   .inferredDistribution = NULL,
                                   .createVirtualPopulation = function(demographyRanges,
                                                                       numberOfVirtualIndividuals,
                                                                       proportionOfFemales,
                                                                       numberOfClusters){
                                     clusters <- getClusters(private$.inferredDistribution,
                                                             numberOfClusters)
                                     populationDataFrames <- createVirtualPopulation(private$.inferredDistribution,
                                                                                     proportionOfFemales,
                                                                                     numberOfVirtualIndividuals,
                                                                                     clusters,
                                                                                     demographyRanges)
                                     private$.referencePopulationDataframe <- populationDataFrames$referencePopulationDataframe
                                     private$.testPopulationDataframe <- populationDataFrames$testPopulationDataframe
                                   },
                                   .referencePopulationDataframe = NULL,
                                   .testPopulationDataframe = NULL,
                                   .iovListSD = list(),
                                   .iovListCV = list(),
                                   .checkIOVParameterReplacement = function(parameterPath,variabilityType){
                                     if(parameterPath %in% names(private$.iovListSD)){
                                       warning(paste0("Replacing additive variability for parameter path ",parameterPath," with new ",variabilityType," variability."))
                                       private$.iovListSD[[parameterPath]] <- NULL
                                     }
                                     if(parameterPath %in% names(private$.iovListCV)){
                                       warning(paste0("Replacing multiplicative variability for parameter path ",parameterPath," with new ",variabilityType," variability."))
                                       private$.iovListCV[[parameterPath]] <- NULL
                                     }
                                   }
                                 ))

# rdf <- data.frame(A = c(1,2), B = c(1,2))
# tdf <- data.frame(A = c(10,20), B = c(10,20))
# ii <- IOV$new(referencePopulationDataframe =rdf , testPopulationDataframe = tdf)
# ii$addAdditiveIOV(parameterPath = "A",SD = 2)


# demographyRanges <- list(weight = list(range = c(50, 110),
#                                        units = "kg"),
#                          height = list(range = c(16,20),
#                                        units = "dm"))

#' @title createVirtualPopulation
#' @description Create a virtual population.
#' @param inferredDistribution The inferred distribution.
#' @param proportionOfFemales Proportion of females in the virtual population in percent
#' @param numberOfVirtualIndividuals The number of virtual individuals to create.
#' @param clusters List of clusters.
#' @param demographyRanges The demography ranges.
#' @importFrom stats rlnorm
createVirtualPopulation <- function(inferredDistribution,
                                    proportionOfFemales,
                                    numberOfVirtualIndividuals,
                                    clusters,
                                    demographyRanges) {
  # create virtual population
  cli::cli_h1("Virtual Population")
  cli::cli_progress_step("Create virtual population characteristics")

  if ("cofactorPaths" %in% names(inferredDistribution)) {
    populationDataframe <- data.frame(IndividualId = seq(0, (numberOfVirtualIndividuals - 1))) # data.frame(matrix(ncol = 0, nrow = numberOfVirtualIndividuals))

    for (cofactor in names(inferredDistribution$cofactorPaths)) {
      geomean <- mean(c(log(min(demographyRanges[[cofactor]]$range)), log(max(demographyRanges[[cofactor]]$range))))
      geosd <- (log(max(demographyRanges[[cofactor]]$range)) - log(min(demographyRanges[[cofactor]]$range))) / (2 * 1.96)
      values <- stats::rlnorm(n = numberOfVirtualIndividuals, meanlog = geomean, sdlog = geosd)
      populationDataframe[[cofactorPaths[[cofactor]]]] <- toBaseUnit(
        quantityOrDimension = cofactorDimensions[[cofactor]],
        values = values,
        unit = demographyRanges[[cofactor]]$units
      )
    }

    cofactorPathDictionary <- cofactorPaths
  } else {
    weightRange <- demographyRanges$weight$range
    weightUnits <- demographyRanges$weight$units
    heightRange <- demographyRanges$height$range
    heightUnits <- demographyRanges$height$units
    ageRange <- demographyRanges$age$range
    ageUnits <- demographyRanges$age$units
    gestationalAgeRange <- demographyRanges$gestationalAge$range
    gestationalAgeUnits <- demographyRanges$gestationalAge$units
    BMIRange <- demographyRanges$BMI$range
    BMIUnits <- demographyRanges$BMI$units

    virtualPopnChars <- createPopulationCharacteristics(
      species = unique(inferredDistribution$demographicData$species),
      population = unique(inferredDistribution$demographicData$population),
      numberOfIndividuals = numberOfVirtualIndividuals,
      proportionOfFemales = proportionOfFemales,
      weightMin = weightRange[1] %||% minNull(inferredDistribution$demographicData[[standardColumnNames$weightColumn]]),
      weightMax = weightRange[2] %||% maxNull(inferredDistribution$demographicData[[standardColumnNames$weightColumn]]),
      weightUnit = weightUnits %||% inferredDistribution$cofactorUnits[[standardColumnNames$weightColumn]] %||% ospsuite::ospUnits$Mass$kg,
      heightMin = heightRange[1] %||% minNull(inferredDistribution$demographicData[[standardColumnNames$heightColumn]]),
      heightMax = heightRange[2] %||% maxNull(inferredDistribution$demographicData[[standardColumnNames$heightColumn]]),
      heightUnit = heightUnits %||% inferredDistribution$cofactorUnits[[standardColumnNames$heightColumn]] %||% ospsuite::ospUnits$Length$cm,
      ageMin = ageRange[1] %||% minNull(inferredDistribution$demographicData[[standardColumnNames$ageColumn]]),
      ageMax = ageRange[2] %||% maxNull(inferredDistribution$demographicData[[standardColumnNames$ageColumn]]),
      ageUnit = ageUnits %||% inferredDistribution$cofactorUnits[[standardColumnNames$ageColumn]] %||% ospsuite::ospUnits$`Age in years`$`year(s)`,
      gestationalAgeMin = gestationalAgeRange[1] %||% minNull(inferredDistribution$demographicData[[standardColumnNames$gestationalAgeColumn]]),
      gestationalAgeMax = gestationalAgeRange[2] %||% maxNull(inferredDistribution$demographicData[[standardColumnNames$gestationalAgeColumn]]),
      gestationalAgeUnit = gestationalAgeUnits %||% inferredDistribution$cofactorUnits[[standardColumnNames$gestationalAgeColumn]] %||% ospsuite::ospUnits$`Age in weeks`$`week(s)`,
      BMIMin = BMIRange[1],
      BMIMax = BMIRange[2],
      BMIUnit = BMIUnits %||% ospsuite::ospUnits$BMI$`kg/m²`
    )

    virtualPopn <- createPopulation(populationCharacteristics = virtualPopnChars)
    populationDataframe <- populationToDataFrame(population = virtualPopn$population)
    cofactorPathDictionary <- individualParameterPaths
  }
  referencePopulationDataframe <- populationDataframe
  testPopulationDataframe <- populationDataframe

  # update the parameter paths to match new simulation
  parameterList <- getParameterListFunctions[[inferredDistribution$method]](inferredDistribution)
  referenceSimulationParameterPaths <- getParameterPathsInReferenceSimulation(parameterList)
  testSimulationParameterPaths <- getParameterPathsInTestSimulation(parameterList)

  parameterPaths <- names(testPopulationDataframe)
  for (j in seq_along(parameterList)) {
    parameterPaths[parameterPaths == referenceSimulationParameterPaths[j]] <- testSimulationParameterPaths[j]
  }
  colnames(testPopulationDataframe) <- parameterPaths

  cofactorNames <- intersect(allCofactorNames, names(inferredDistribution$demographicData))

  # find theta values of closest individual in reference population
  i <- 1
  cli::cli_progress_step(
    "Sampling from conditional distribution: individual {i}/{numberOfVirtualIndividuals}",
    spinner = TRUE
  )
  for (i in 1:numberOfVirtualIndividuals) {
    # read in ith person's weight and height
    cofactorValues <- sapply(cofactorNames, function(cof) {
      dimension <- cofactorDimensions[[cof]]
      ospsuite::toUnit(
        quantityOrDimension = dimension,
        values = populationDataframe[[cofactorPathDictionary[[cof]]]][i],
        sourceUnit = ospsuite::getBaseUnit(dimension),
        targetUnit = inferredDistribution$cofactorUnits[[cof]]
      )
    })

    givenPoints <- c(rep(NA, length(parameterList)), cofactorValues)

    # Sample from conditional distribution and ensure that sample is within parameter bounds
    thetaSamples <- samplesFromMixture(
      mclstResults = clusters,
      givenPoints = givenPoints,
      numberOfSamples = 1,
      lowerBounds = sapply(parameterList, function(x) {
        x$lowerBound
      }),
      upperBounds = sapply(parameterList, function(x) {
        x$upperBound
      })
    )
    cli::cli_progress_update()

    for (j in seq_along(parameterList)) {
      referencePopulationDataframe[[referenceSimulationParameterPaths[j]]][i] <- thetaSamples[, j]
      testPopulationDataframe[[testSimulationParameterPaths[j]]][i] <- thetaSamples[, j]
    }
  }
  cli::cli_progress_done()

  return(list(referencePopulationDataframe = referencePopulationDataframe, testPopulationDataframe = testPopulationDataframe))
}

scaleLognormally <- function(params, CV) {
  SD <- sqrt(log(1 + CV^2))
  Z <- rnorm(length(params),sd = SD)
  scaledParams <- params * exp(Z)
  return(scaledParams)
}

perturbNormally <- function(params, SD) {
  Z <- rnorm(length(params),sd = SD)
  perturbedParams <- params + Z
  return(perturbedParams)
}

#' @title simulateVirtualPopulation
#' @description Simulate a virtual population.
#' @param referenceSimulationFilePath The path to the reference simulation file.
#' @param testSimulationFilePath The path to the test simulation file.
#' @param referencePopulationDataframe The reference population dataframe.
#' @param testPopulationDataframe The test population dataframe.
#' @param outputPath The output path.
#' @param startTime The start time of the simulation.
#' @param endTime The end time of the simulation.
#' @param resolutionPtsMin The resolution of the simulation.
#' @import ospsuite
simulateVirtualPopulation <- function(referenceSimulationFilePath,
                                      testSimulationFilePath,
                                      referencePopulationDataframe,
                                      testPopulationDataframe,
                                      numberOfReplicates,
                                      outputPath,
                                      startTime,
                                      endTime,
                                      resolutionPtsMin,
                                      iovListSD,
                                      iovListCV) {
  #load simulations
  cli::cli_progress_step("Simulating virtual population", spinner = TRUE)
  referenceSimulation <- ospsuite::loadSimulation(referenceSimulationFilePath)
  referenceSimulation$outputSelections$clear()
  ospsuite::addOutputs(quantitiesOrPaths = outputPath, simulation = referenceSimulation)
  cli::cli_progress_update()

  testSimulation <- ospsuite::loadSimulation(testSimulationFilePath)
  testSimulation$outputSelections$clear()
  ospsuite::addOutputs(quantitiesOrPaths = outputPath, simulation = testSimulation)
  cli::cli_progress_update()

  # set output interval
  ospsuite::setOutputInterval(simulation = referenceSimulation, startTime = startTime, endTime = endTime, resolution = resolutionPtsMin)
  ospsuite::setOutputInterval(simulation = testSimulation, startTime = startTime, endTime = endTime, resolution = resolutionPtsMin)
  cli::cli_progress_update()

  # add iov parameters to population dataframe if missing
  for (iovParameterPath in unique(c(names(iovListCV),names(iovListSD)))){

    if(!(iovParameterPath %in% names(referencePopulationDataframe))){
      referencePopulationDataframe[[iovParameterPath]] <- getParameter(path = iovParameterPath,container = referenceSimulation)$value
    }

    if(!(iovParameterPath %in% names(testPopulationDataframe))){
      testPopulationDataframe[[iovParameterPath]] <- getParameter(path = iovParameterPath,container = testSimulation)$value
    }

  }

  combinedResultsData <- NULL

  for (rep in seq_len(numberOfReplicates)){

    referencePopulationDataframeWithIOV <- referencePopulationDataframe
    testPopulationDataframeWithIOV <- testPopulationDataframe

    #perturb multiplicative IOV parameters in the dataframe
    for (iovParameterPath in names(iovListCV)){

      referencePopulationDataframeWithIOV[[iovParameterPath]] <- scaleLognormally(params = referencePopulationDataframeWithIOV[[iovParameterPath]],
                                                                                  CV = iovListCV[[iovParameterPath]])
      testPopulationDataframeWithIOV[[iovParameterPath]] <- scaleLognormally(params = testPopulationDataframeWithIOV[[iovParameterPath]],
                                                                             CV = iovListCV[[iovParameterPath]])

    }

    #perturb additive IOV parameters in the dataframe
    for (iovParameterPath in names(iovListSD)){
      referencePopulationDataframeWithIOV[[iovParameterPath]] <- perturbNormally(params = referencePopulationDataframeWithIOV[[iovParameterPath]],
                                                                                 SD = iovListSD[[iovParameterPath]])
      testPopulationDataframeWithIOV[[iovParameterPath]] <- perturbNormally(params = testPopulationDataframeWithIOV[[iovParameterPath]],
                                                                            SD = iovListSD[[iovParameterPath]])
    }

    write.csv(referencePopulationDataframeWithIOV, file = "referenceVirtualPopulation.csv", row.names = FALSE)
    write.csv(testPopulationDataframeWithIOV, file = "testVirtualPopulation.csv", row.names = FALSE)

    referenceVirtualPopulation <- ospsuite::loadPopulation("referenceVirtualPopulation.csv")
    testVirtualPopulation <- ospsuite::loadPopulation("testVirtualPopulation.csv")

    # generate plasma concentration time profiles
    referenceSimulationResults <- ospsuite::runSimulation(simulation = referenceSimulation, population = referenceVirtualPopulation)
    testSimulationResults <- ospsuite::runSimulation(simulation = testSimulation, population = testVirtualPopulation)
    cli::cli_progress_update()

    referenceResultsData <- ospsuite::getOutputValues(referenceSimulationResults, quantitiesOrPaths = outputPath)$data
    referenceResultsData$drug <- "R"
    referenceResultsData$replicate <- rep
    cli::cli_progress_update()

    testResultsData <- ospsuite::getOutputValues(testSimulationResults, quantitiesOrPaths = outputPath)$data
    testResultsData$drug <- "T1"
    testResultsData$replicate <- rep
    cli::cli_progress_update()

    combinedResultsData <- rbind(combinedResultsData,referenceResultsData, testResultsData)


    resultsData <- data.frame(
      id = 1 + combinedResultsData[["IndividualId"]],
      time = combinedResultsData[["Time"]],
      conc = combinedResultsData[[outputPath]],
      rep = combinedResultsData[["replicate"]],
      drug = combinedResultsData[["drug"]]
    )

    resultsData$id <- resultsData$id
    resultsData$per <- resultsData$per
    resultsData$drug <- as.factor(resultsData$drug)

    file.remove("referenceVirtualPopulation.csv")
    file.remove("testVirtualPopulation.csv")

  }
  cli::cli_progress_done()
  return(resultsData)
}

#' @title runClinicalTrialSimulation
#' @description Run a clinical trial simulation.
#' @param virtualPopulationSimulationResults The results of the virtual population simulation.
#' @param n_trials The number of trials to run.
#' @param subj_min The minimum number of subjects to include in the trial.
#' @param subj_max The maximum number of subjects to include in the trial.
#' @param subj_step The step size between `subj_min` and `subj_max`.
#' @export
runClinicalTrialSimulation <- function(virtualPopulationSimulationResults,
                                       n_trials = 50,
                                       subj_min = 5,
                                       subj_max = 50,
                                       subj_step = 5,
                                       aucICV = 0,
                                       cmaxICV = 0) {
  cli::cli_progress_step("Running clinical trial simulation", spinner = TRUE)
  res <- calc_be(virtualPopulationSimulationResults, n_trials = n_trials, subj_min = subj_min, subj_max = subj_max, subj_step = subj_step, aucICV = aucICV, cmaxICV = cmaxICV)
  cli::cli_progress_update()
  extRes <- extract_be(res)
  cli::cli_progress_done()
  # make_report(sim=virtualPopulationSimulationResults,
  #             be=extRes)
  return(extRes)
}
