rm(list=ls())
graphics.off()
library(ospsuite.openVBE)
dateTime <- paste0(format(Sys.Date(), "%Y%m%d"), "_", format(Sys.time(), "%H%M%S"))
exampleName <- "testosterone-dermal"
subfolder <- file.path("examples",exampleName)

#Data set D1
pkDataFilePath <- file.path(subfolder,"testosterone_pk_data_reference.csv")

#Data set D2
studyPopulationDataFilePath <- file.path(subfolder,"testosterone_study_population_data.csv")

########>   STEP S1   <########
referenceSimulationFilePath <- file.path(subfolder,"testosterone-dermal-invivo-reference-petrolatum-2ug-per-cm2.pkml")
testSimulationFilePath <- file.path(subfolder,"testosterone-dermal-invivo-test-ethylene-glycol-2ug-per-cm2.pkml")


########>   STEP S2   <########
outputPath <- "DERMAL_APPLICATION_AREA|in_vivo_sink|permeant|whole_body_concentration"

uncertainty_logktrans <- VBEParameter$new(pathInReferenceSimulation = "permeant|uncertainty_logktrans",
                                          pathInTestSimulation = "permeant|uncertainty_logktrans",
                                          displayName = "uncertainty_logktrans",
                                          dimension = ospDimensions$Dimensionless,
                                          unit = ospUnits$Dimensionless$Unitless,
                                          lowerBound = -1,
                                          upperBound = 1,
                                          logIncrement = FALSE)

scaling_total_clearance <- VBEParameter$new(pathInReferenceSimulation = "in_vivo_sink_clearance|scaling_total_clearance",
                                            pathInTestSimulation = "in_vivo_sink_clearance|scaling_total_clearance",
                                            displayName = "scaling_total_clearance",
                                            dimension = ospDimensions$Dimensionless,
                                            unit = ospUnits$Dimensionless$Unitless,
                                            lowerBound = 0.5,
                                            upperBound = 1.5,
                                            logIncrement = FALSE)

inferenceParameters <- list(uncertainty_logktrans,
                            scaling_total_clearance)

cofactorPaths <- list(weight = "DERMAL_APPLICATION_AREA|in_vivo_sink|body_weight")


inferredDistribution <- runDistributionInference(method = "NPOD",
                                                 referenceSimulationFilePath = referenceSimulationFilePath,
                                                 outputPath = outputPath,
                                                 pkDataFilePath = pkDataFilePath,
                                                 studyPopulationDataFilePath = studyPopulationDataFilePath,
                                                 cofactorPaths = cofactorPaths,
                                                 inferenceParameters = inferenceParameters,
                                                 initialGridSize = 40,
                                                 numberOfIterations = 10,
                                                 useLogNormalLikelihood = TRUE,
                                                 saveResultsPath = file.path(subfolder,paste0(exampleName,"-npod-results-",dateTime,".rds")),
                                                 cacheFolder = file.path("examples/cacheFolder"))


########>   STEP 3   <########

## Parameters to describe demographics (body wieght) of study population to create a virtual population
# Given values
popBWMedian <- 85 #Median body weight
popBWCV <- 0.148 #CV body weight
# Compute lognormal parameters
mu <- log(popBWMedian)
sigma <- sqrt(log(1 + popBWCV^2))
# Compute confidence intervals
CI_95 <- qlnorm(c(0.025, 0.975), meanlog = mu, sdlog = sigma)

demographyRanges <- list(weight = list(range = c(CI_95[1], CI_95[2]),
                                       units = "kg"))

virtualPopulation <- createVirtualPopulation(inferredDistribution = inferredDistribution,
                                             numberOfVirtualIndividuals = 100,
                                             clusters = getClusters(inferredDistribution = inferredDistribution,numberOfClusters = 4),
                                             demographyRanges = demographyRanges)

refAndTestSimulationsInVirtualPopulation <- simulateVirtualPopulation(referenceSimulationFilePath = referenceSimulationFilePath,
                                                                      testSimulationFilePath = testSimulationFilePath,
                                                                      referencePopulationDataframe = virtualPopulation$referencePopulationDataframe,
                                                                      testPopulationDataframe = virtualPopulation$testPopulationDataframe,
                                                                      outputPath = outputPath,
                                                                      startTime = 0,
                                                                      endTime = 6000,
                                                                      resolutionPtsMin =1/15)  # pts/min


########>   PLOTTING   <########
source(file.path(subfolder,"plotting-script-vbe-testosterone-dermal.R"), encoding = 'UTF-8')


########>   STEP S4   <########
ctsResults <- runClinicalTrialSimulation(virtualPopulationSimulationResults = refAndTestSimulationsInVirtualPopulation,n_trials = 1000)

pltList <- plotClinicalTrialSimulationResults(ctsResults)
show(pltList)

bandsPlot <- plotVBEBands(ctsResults)
show(bandsPlot)

probabilityPlot <- plotVBEProbability(ctsResults)
show(probabilityPlot)



