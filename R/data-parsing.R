#' @export
standardColumnNames <- list(
  idColumn = "id",
  occasionColumn = "occasion",
  timeColumn = "time",
  timeUnitsColumn = "timeUnits",
  outputNamesColumn = "outputNames",
  outputValuesColumn = "outputValues",
  outputUnitsColumn = "outputUnits",
  speciesColumn = "species",
  populationColumn = "population",
  genderColumn = "gender",
  weightColumn = "weight",
  weightUnitsColumn = "weightUnits",
  heightColumn = "height",
  heightUnitsColumn = "heightUnits",
  ageColumn = "age",
  ageUnitsColumn = "ageUnits",
  gestationalAgeColumn = "gestationalAge",
  gestationalAgeUnitsColumn = "gestationalAgeUnits"
)


#' @export
experimentFileInputs <- c("idColumn", "occasionColumn", "timeColumn", "timeUnitsColumn", "outputNamesColumn", "outputValuesColumn", "outputUnitsColumn")


#' @export
populationFileInputs <- c("idColumn", "speciesColumn", "populationColumn", "genderColumn", "weightColumn", "weightUnitsColumn", "heightColumn", "heightUnitsColumn", "ageColumn", "ageUnitsColumn", "gestationalAgeColumn", "gestationalAgeUnitsColumn")

individualParameterPaths <- list()
individualParameterPaths[[standardColumnNames$weightColumn]] <- "Organism|Weight"
individualParameterPaths[[standardColumnNames$heightColumn]] <- "Organism|Height"
individualParameterPaths[[standardColumnNames$ageColumn]] <- "Organism|Age"
individualParameterPaths[[standardColumnNames$gestationalAgeColumn]] <- "Organism|Gestational age"

allCofactorNames <- c(
  standardColumnNames$weightColumn,
  standardColumnNames$heightColumn,
  standardColumnNames$ageColumn,
  standardColumnNames$gestationalAgeColumn
)

cofactorDimensions <- list()
cofactorDimensions[[standardColumnNames$weightColumn]] <- ospDimensions$Mass
cofactorDimensions[[standardColumnNames$heightColumn]] <- ospDimensions$Length
cofactorDimensions[[standardColumnNames$ageColumn]] <- ospDimensions$`Age in years`
cofactorDimensions[[standardColumnNames$gestationalAgeColumn]] <- ospDimensions$`Age in weeks`

#' @export
processPopulationDataFile <- function(populationFilePathOrDataframe,
                                      idColumn,
                                      speciesColumn,
                                      populationColumn,
                                      genderColumn = NULL,
                                      weightColumn = NULL,
                                      weightUnitsColumn = NULL,
                                      heightColumn = NULL,
                                      heightUnitsColumn = NULL,
                                      ageColumn = NULL,
                                      ageUnitsColumn = NULL,
                                      gestationalAgeColumn = NULL,
                                      gestationalAgeUnitsColumn = NULL) {
  dataDf <- populationFilePathOrDataframe
  if (is.character(populationFilePathOrDataframe)) {
    dataDf <- read.csv(populationFilePathOrDataframe,
      check.names = FALSE,
      fileEncoding = "UTF-8-BOM"
    )
  }
  colnames(dataDf)[colnames(dataDf) == idColumn] <- standardColumnNames$idColumn
  colnames(dataDf)[colnames(dataDf) == speciesColumn] <- standardColumnNames$speciesColumn
  colnames(dataDf)[colnames(dataDf) == populationColumn] <- standardColumnNames$populationColumn
  colnames(dataDf)[colnames(dataDf) == genderColumn] <- standardColumnNames$genderColumn
  colnames(dataDf)[colnames(dataDf) == weightColumn] <- standardColumnNames$weightColumn
  colnames(dataDf)[colnames(dataDf) == weightUnitsColumn] <- standardColumnNames$weightUnitsColumn
  colnames(dataDf)[colnames(dataDf) == heightColumn] <- standardColumnNames$heightColumn
  colnames(dataDf)[colnames(dataDf) == heightUnitsColumn] <- standardColumnNames$heightUnitsColumn
  colnames(dataDf)[colnames(dataDf) == ageColumn] <- standardColumnNames$ageColumn
  colnames(dataDf)[colnames(dataDf) == ageUnitsColumn] <- standardColumnNames$ageUnitsColumn
  colnames(dataDf)[colnames(dataDf) == gestationalAgeColumn] <- standardColumnNames$gestationalAgeColumn
  colnames(dataDf)[colnames(dataDf) == gestationalAgeUnitsColumn] <- standardColumnNames$gestationalAgeUnitsColumn

  validColumnNames <- sapply(populationFileInputs, function(x) {
    standardColumnNames[[x]]
  })
  dataDf <- unique(dataDf[, intersect(colnames(dataDf), validColumnNames)])

  ## Validate no more than one row per individual id

  return(dataDf)
}



#' @export
processExperimentDataFile <- function(experimentDataFilePathOrDataframe,
                                      idColumn,
                                      occasionColumn,
                                      timeColumn,
                                      timeUnitsColumn,
                                      outputNamesColumn,
                                      outputValuesColumn,
                                      outputUnitsColumn) {
  dataDf <- experimentDataFilePathOrDataframe
  if (is.character(experimentDataFilePathOrDataframe)) {
    dataDf <- read.csv(experimentDataFilePathOrDataframe,
      check.names = FALSE,
      fileEncoding = "UTF-8-BOM"
    )
  }

  colnames(dataDf)[colnames(dataDf) == idColumn] <- standardColumnNames$idColumn
  colnames(dataDf)[colnames(dataDf) == occasionColumn] <- standardColumnNames$occasionColumn
  colnames(dataDf)[colnames(dataDf) == timeColumn] <- standardColumnNames$timeColumn
  colnames(dataDf)[colnames(dataDf) == timeUnitsColumn] <- standardColumnNames$timeUnitsColumn
  colnames(dataDf)[colnames(dataDf) == outputNamesColumn] <- standardColumnNames$outputNamesColumn
  colnames(dataDf)[colnames(dataDf) == outputValuesColumn] <- standardColumnNames$outputValuesColumn
  colnames(dataDf)[colnames(dataDf) == outputUnitsColumn] <- standardColumnNames$outputUnitsColumn
  validColumnNames <- sapply(experimentFileInputs, function(x) {
    standardColumnNames[[x]]
  })
  dataDf <- dataDf[, intersect(colnames(dataDf), validColumnNames)]

  return(dataDf)
}
