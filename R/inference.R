#' @title runDistributionInference
#' @description Run distribution inference
#' @param method The method to use for distribution inference. Currently only `"NPOD"` is supported.
#' @param referenceSimulationFilePath The path to the reference simulation file.
#' @param outputPath The path to the output directory.
#' @param pkDataFilePath The path to the PK data file.
#' @param studyPopulationDataFilePath The path to the study population data file.
#' @param cofactorPaths The paths to the cofactor files.
#' @param inferenceParameters The parameters to infer.
#' @param numberOfIterations The number of iterations to run.
#' @param useLogNormalLikelihood Whether to use a log-normal likelihood.
#' @param initialGridSize The initial grid size.
#' @param saveResultsPath The path to save the results.
#' @param cacheFolder The folder to use for caching.
#' @export
runDistributionInference <- function(method = "NPOD",
                                     referenceSimulationFilePath,
                                     outputPath,
                                     pkDataFilePath,
                                     studyPopulationDataFilePath = NULL,
                                     cofactorPaths = NULL,
                                     inferenceParameters = NULL,
                                     numberOfIterations = 10,
                                     useLogNormalLikelihood = NULL,
                                     initialGridSize = NULL,
                                     saveResultsPath = NULL,
                                     cacheFolder = NULL) {
  return(runDistributionInferenceFunctionsList[[method]](referenceSimulationFilePath = referenceSimulationFilePath,
    outputPath = outputPath,
    pkDataFilePath = pkDataFilePath,
    studyPopulationDataFilePath = studyPopulationDataFilePath,
    cofactorPaths = cofactorPaths,
    inferenceParameters = inferenceParameters,
    numberOfIterations = numberOfIterations,
    useLogNormalLikelihood = useLogNormalLikelihood,
    initialGridSize = initialGridSize,
    saveResultsPath = saveResultsPath,
    cacheFolder = cacheFolder
  ))
}

runDistributionInferenceNPOD <- function(referenceSimulationFilePath,
                                         outputPath,
                                         pkDataFilePath,
                                         studyPopulationDataFilePath,
                                         cofactorPaths,
                                         inferenceParameters,
                                         numberOfIterations,
                                         useLogNormalLikelihood,
                                         initialGridSize,
                                         saveResultsPath,
                                         cacheFolder,
                                         ...) {
  cli::cli_h1("Non-Parametric Optimal Design (NPOD)")
  
  individualCharacteristics <- NULL
  cofactorUnits <- NULL
  if (!is.null(studyPopulationDataFilePath)) {
    demographicData <- read.csv(studyPopulationDataFilePath)
    unifiedCofactorsUnits <- unifyPopulationDataUnits(demographicData)
    cofactorUnits <- unifiedCofactorsUnits$cofactorUnits
  }

  npodSettings <- NPODRunSettings$new(
    MaxIter = numberOfIterations,
    cache_folder_name = cacheFolder
  )

  npod <- NPODObject$new(
    simulationFilePath = referenceSimulationFilePath,
    optimizationParameterList = inferenceParameters,
    outputPath = outputPath,
    pkDataFilePath = pkDataFilePath,
    studyPopulationDataFilePath = studyPopulationDataFilePath,
    cofactorPaths = cofactorPaths,
    initialGridSize = initialGridSize,
    useLogNormalLikelihood = useLogNormalLikelihood,
    npodRunSettings = npodSettings
  )

  npodResults <- npod$runNPOD()
  npodResults$parameters <- inferenceParameters
  npodResults$cofactorUnits <- cofactorUnits
  npodResults$cofactorPaths <- cofactorPaths
  npodResults$method <- "NPOD"

  if (!is.null(saveResultsPath)) {
    saveRDS(object = npodResults, file = saveResultsPath)
  }
  return(npodResults)
}


runDistributionInferenceFunctionsList <- list(NPOD = runDistributionInferenceNPOD)
