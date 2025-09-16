#' @title getReferenceGoodnessOfFitPlot
#' @description Generates plots of the goodness of fit for the reference drug.
#' @param refAndTestSimulationsInVirtualPopulation A data frame containing the simulated plasma concentrations for the reference and test drugs.
#' @param pkDataFilePath The path to the observed plasma concentrations data.
#' @param referenceDrugName The name of the reference drug.
#' @param yDimension The dimension of the y-axis.
#' @param yUnit The unit of the y-axis.
#' @param timeUnit The unit of the time axis.
#' @export
#' @import ggplot2
getReferenceGoodnessOfFitPlot <- function(refAndTestSimulationsInVirtualPopulation,
                                          pkDataFilePath,
                                          referenceDrugName = "Reference",
                                          yDimension = NULL,
                                          yUnit = NULL,
                                          timeUnit = ospsuite::ospUnits$Time$min) {
  plasmaDf <- refAndTestSimulationsInVirtualPopulation

  if (!is.null(yDimension) & !is.null(yUnit)) {
    plasmaDf$conc <- ospsuite::toUnit(
      quantityOrDimension = yDimension,
      values = plasmaDf$conc,
      targetUnit = yUnit
    )
  }
  plasmaDf$time <- ospsuite::toUnit(
    quantityOrDimension = ospsuite::ospDimensions$Time,
    values = plasmaDf$time,
    targetUnit = timeUnit
  )

  mmdf <- NULL

  q05 <- function(x) {
    quantile(x = x, 0.05)
  }
  q95 <- function(x) {
    quantile(x = x, 0.95)
  }
  fnVec <- c("median" = median, "min" = min, "max" = max)
  for (fn in names(fnVec)) {
    df <- aggregate(
      x = plasmaDf$conc,
      by = list(
        plasmaDf$time,
        plasmaDf$drug
      ),
      fnVec[[fn]]
    )
    names(df) <- c("time", "drug", fn)
    mmdf[["time"]] <- df[["time"]]
    mmdf[["drug"]] <- df[["drug"]]
    mmdf[[fn]] <- df[[fn]]
  }
  mmdf <- as.data.frame(mmdf)
  mmdf$drug <- as.character(mmdf$drug)
  mmdf$drug[mmdf$drug == "R"] <- referenceDrugName

  plt <- ggplot2::ggplot() +
    ggplot2::geom_ribbon(
      data = mmdf[mmdf$drug == referenceDrugName, ],
      mapping = ggplot2::aes(
        x = time,
        y = (median),
        ymin = (min),
        ymax = (max),
        fill = drug
      ), alpha = 0.5
    ) +
    ggplot2::theme_minimal()

  pkData <- read.csv(pkDataFilePath)
  pkData$drug <- referenceDrugName
  pkData$id <- as.factor(pkData$id)
  plt <- plt + ggplot2::geom_point(
    data = pkData,
    mapping = ggplot2::aes(
      x = time,
      y = outputValues,
      color = id,
      shape = "Individual"
    )
  )

  plt <- plt +
    ggplot2::xlab(paste0("Time (", timeUnit, ")")) +
    ggplot2::ylab(paste0("Reference (", yUnit, ")"))
  plt <- plt + ggplot2::guides(color = "none")
  plt <- plt + ggplot2::theme(
    legend.position = "bottom", legend.direction = "horizontal",
    axis.text = ggplot2::element_text(size = 12), axis.title = ggplot2::element_text(size = 14),
    legend.text = ggplot2::element_text(size = 12), legend.title = ggplot2::element_text(size = 14),
    strip.text = ggplot2::element_text(size = 14)
  )

  plt <- plt + scale_fill_manual(name = referenceDrugName, labels = "Simulated", values = "darkgrey")
  plt <- plt + scale_shape_manual(name = "", labels = "Observed", values = 16)

  methods::show(plt)

  return(plt)
}


#' @title getReferenceVBEPlot
#' @description Generates plots of the virtual bioequivalence for the reference drug.
#' @param refAndTestSimulationsInVirtualPopulation A data frame containing the simulated plasma concentrations for the reference and test drugs.
#' @param referenceDrugName The name of the reference drug.
#' @param testDrugName The name of the test drug.
#' @param yDimension The dimension of the y-axis.
#' @param yUnit The unit of the y-axis.
#' @param timeUnit The unit of the time axis.
#' @export
getReferenceVBEPlot <- function(refAndTestSimulationsInVirtualPopulation,
                                referenceDrugName = "Reference",
                                testDrugName = "Test",
                                yDimension = NULL,
                                yUnit = NULL,
                                timeUnit = ospsuite::ospUnits$Time$min) {
  plasmaDf <- refAndTestSimulationsInVirtualPopulation

  if (!is.null(yDimension) & !is.null(yUnit)) {
    plasmaDf$conc <- ospsuite::toUnit(
      quantityOrDimension = yDimension,
      values = plasmaDf$conc,
      targetUnit = yUnit
    )
  }
  plasmaDf$time <- ospsuite::toUnit(
    quantityOrDimension = ospsuite::ospDimensions$Time,
    values = plasmaDf$time,
    targetUnit = timeUnit
  )

  vbemmdf <- NULL

  q05 <- function(x) {
    quantile(x = x, 0.05)
  }
  q95 <- function(x) {
    quantile(x = x, 0.95)
  }
  fnVec <- c("median" = median, "min" = min, "max" = max)
  for (fn in names(fnVec)) {
    df <- aggregate(
      x = plasmaDf$conc,
      by = list(
        plasmaDf$time,
        plasmaDf$drug
      ),
      fnVec[[fn]]
    )
    names(df) <- c("time", "drug", fn)
    vbemmdf[["time"]] <- df[["time"]]
    vbemmdf[["drug"]] <- df[["drug"]]
    vbemmdf[[fn]] <- df[[fn]]
  }
  vbemmdf <- as.data.frame(vbemmdf)
  vbemmdf$drug <- as.character(vbemmdf$drug)
  vbemmdf$drug[vbemmdf$drug == "R"] <- referenceDrugName
  vbemmdf$drug[vbemmdf$drug == "T1"] <- testDrugName

  vbemmdf$drug <- factor(vbemmdf$drug, levels = c(referenceDrugName, testDrugName))
  vbeplt <- ggplot2::ggplot() +
    ggplot2::geom_ribbon(
      data = vbemmdf,
      mapping = ggplot2::aes(
        x = time,
        y = (median),
        ymin = (min),
        ymax = (max),
        fill = drug
      ), alpha = 0.33
    )

  vbeplt <- vbeplt + ggplot2::geom_line(
    data = vbemmdf,
    mapping = ggplot2::aes(
      x = time,
      y = (median),
      color = drug
    ),
    size = 2
  ) + ggplot2::theme_minimal()

  vbeplt <- vbeplt + ggplot2::theme(
    legend.position = "bottom", legend.direction = "horizontal",
    axis.text = ggplot2::element_text(size = 12), axis.title = ggplot2::element_text(size = 14),
    legend.text = ggplot2::element_text(size = 12), legend.title = ggplot2::element_text(size = 14),
    strip.text = ggplot2::element_text(size = 14)
  )

  vbeplt <- vbeplt +
    ggplot2::xlab(paste0("Time (", timeUnit, ")")) +
    ggplot2::ylab(paste0(yDimension, " (", yUnit, ")"))
  vbeplt <- vbeplt + ggplot2::scale_fill_manual(name = "Drug", labels = c(referenceDrugName, testDrugName), values = c("#ff0000", "#1e90ff"))
  vbeplt <- vbeplt + ggplot2::scale_color_manual(name = "Drug", labels = c(referenceDrugName, testDrugName), values = c("#ff0000", "#1e90ff"))

  methods::show(vbeplt)
  return(vbeplt)
}


#' @title plotVBEProbability
#' @description Generates plots of the probability of bioequivalence.
#' @param ctsResults The results of the clinical trial simulation.
#' @export
plotVBEProbability <- function(ctsResults) {

  ctsResultsSummary <- ctsResults$summary

  trialNames <- c("auc_par", "cmax_par", "auc_cross", "cmax_cross","auc_rep","cmax_rep")

  pkList <- c(sapply(c("auc_par", "auc_cross","auc_rep"),function(x){"AUC"}),
              sapply(c("cmax_par", "cmax_cross","cmax_rep"),function(x){"Cmax"}))

  designList <- c(sapply(c("auc_par", "cmax_par"),function(x){"Parallel"}),
                  sapply(c("cmax_cross","auc_cross"),function(x){"Crossover"}),
                  sapply(c("auc_rep","cmax_rep"),function(x){"Replicate"}))

  probSuccessDfColNames <- c("n_subj","pBE_rep")

  probSuccessDf <- NULL

  for (trial in trialNames) {

    if (!is.data.frame(ctsResultsSummary[[trial]])) {
      next
    }

    df <- as.data.frame(ctsResultsSummary[[trial]])

    df <- df[,c(probSuccessDfColNames)]

    df$trial <- trial
    probSuccessDf <- rbind.data.frame(probSuccessDf, df)
  }

  probSuccessDf$pk <- sapply(probSuccessDf$trial, function(x) {
    pkList[[x]]
  })

  probSuccessDf$design <- sapply(probSuccessDf$trial, function(x) {
    designList[[x]]
  })

  proPlt <- ggplot2::ggplot(data = probSuccessDf, mapping = ggplot2::aes(x = n_subj, y = pBE_rep)) +
    ggplot2::geom_line() +
    ggplot2::geom_point() +
    ggplot2::facet_grid(pk ~ design) +
    ggplot2::scale_x_continuous(name = "Number of subjects") +
    ggplot2::scale_y_continuous(name = "Probability of BE")
  proPlt <- proPlt + ggplot2::geom_hline(yintercept = 0.9, linetype = "dashed")
  proPlt <- proPlt + ggplot2::annotate("text", x = 15, y = 0.82, label = "0.9")
  proPlt <- proPlt + ggplot2::theme(
    axis.title = ggplot2::element_text(size = 14),
    axis.text = ggplot2::element_text(size = 12),
    strip.text = ggplot2::element_text(size = 16)
  )
  return(proPlt)
}

getMeanMinMaxDataframe <- function(dataframe, valuesColumnName, byColumnNames,
                                   returnGeometricMean = FALSE,
                                   yminFunction = min, ymaxFunction = max) {
  # Define mean function
  meanFn <- if (returnGeometricMean) {
    function(x) exp(mean(log(x), na.rm = TRUE))
  } else {
    function(x) mean(x, na.rm = TRUE)
  }

  # Create list of aggregation functions
  aggFunctions <- list(ymean = meanFn, ymin = yminFunction, ymax = ymaxFunction)

  # Aggregate data
  aggregated <- lapply(aggFunctions, function(fn) {
    aggregate(dataframe[[valuesColumnName]],
              by = dataframe[byColumnNames],
              FUN = fn
    )
  })

  # Combine results into a single dataframe
  result <- aggregated[[1]][, byColumnNames, drop = FALSE] # Base structure
  for (name in names(aggFunctions)) {
    result[[name]] <- aggregated[[name]]$x
  }

  return(result)
}


#' @title plotVBEBands
#' @description Generates plots of the virtual bioequivalence bands.
#' @param ctsResults The results of the clinical trial simulation.
#' @export
plotVBEBands <- function(ctsResults) {
  loEnd <- 0.05

  yminFunction <- function(x) {
    quantile(x = x, loEnd)
  }
  ymaxFunction <- function(x) {
    quantile(x = x, 1 - loEnd)
  }

  dfAll <- ctsResults$trials
  if ("trials" %in% names(dfAll)) {
    dfAll <- dfAll[["trials"]]
  }

  trialLabels <- c(
    "auc_par_all" = paste0("AUC (Parallel)"),
    "auc_cross_all" = paste0("AUC (Crossover)"),
    "auc_rep_all" = paste0("AUC (Replicate)"),
    "cmax_par_all" = paste0("Cmax (Parallel)"),
    "cmax_cross_all" = paste0("Cmax (Crossover)"),
    "cmax_rep_all" = paste0("Cmax (Replicate)")
  )

  trialNames <- names(trialLabels)

  pltList <- list()

  for (trial in trialNames) {
    df <- dfAll[[trial]]

    if(!is.data.frame(df)){
      next
    }

    df <- reshape2::melt(data = df, measure.vars = c("de_ratio", "ci_lo", "ci_up"))

    df <- getMeanMinMaxDataframe(
      dataframe = df,
      valuesColumnName = "value",
      byColumnNames = c("n_subj", "variable"),
      returnGeometricMean = TRUE,
      yminFunction = yminFunction,
      ymaxFunction = ymaxFunction
    )

    plt <- ggplot2::ggplot() +
      ggplot2::geom_ribbon(data = df, mapping = ggplot2::aes(x = n_subj, y = ymean, ymin = ymin, ymax = ymax, fill = variable), alpha = 0.33)

    plt <- plt + ggplot2::geom_line(data = df, mapping = aes(x = n_subj, y = ymean, color = variable), size = 0.5)
    plt <- plt + ggplot2::scale_fill_manual(name = "X", values = c(ci_up = "#1e90ff", ci_lo = "#1e90ff", de_ratio = "#ff0000"))
    plt <- plt + ggplot2::scale_color_manual(name = "X", values = c(ci_up = "#1e90ff", ci_lo = "#1e90ff", de_ratio = "#ff0000"))
    plt <- plt + ggplot2::geom_hline(
      yintercept = c(0.8, 1 / 0.8),
      color = "#000000", linetype = "dashed", size = 0.5
    )
    plt <- plt + ggplot2::xlab("Sample size") + ggplot2::ylab("Ratio")
    plt <- plt + ggplot2::theme_minimal()
    plt <- plt + ggplot2::theme(
      axis.text = ggplot2::element_text(size = 12), legend.position = "none",
      panel.background = ggplot2::element_rect(fill = "white", colour = "white"),
      plot.background = ggplot2::element_rect(fill = "white", colour = "white"),
      axis.title = ggplot2::element_text(size = 14), strip.text = ggplot2::element_text(size = 14)
    )

    plt <- plt + ggplot2::ggtitle(label = trialLabels[[trial]])

    pltList[[trial]] <- plt
  }

  return(pltList)
}

#' @title plotClinicalTrialSimulationResults
#' @description Generates plots of the clinical trial simulation results.
#' @param ctsResults The results of the clinical trial simulation.
#' @param addPlotHeaders Whether to add plot headers.
#' @export
plotClinicalTrialSimulationResults <- function(ctsResults, addPlotHeaders = TRUE) {
  pkDictionary <- list(
    auc_par = "AUC (Parallel)",
    cmax_par = "Cmax (Parallel)",
    auc_cross = "AUC (Crossover)",
    cmax_cross = "Cmax (Crossover)",
    auc_rep = "AUC (Replicate)",
    cmax_rep = "Cmax (Replicate)"
  )
  dataDf <- NULL
  pltList <- list()
  for (trial in names(pkDictionary)) {
    df <- ctsResults$summary[[trial]]
    if (!is.data.frame(df)) {
      next
    }
    df$trial <- pkDictionary[[trial]]
    plt <- ggplot2::ggplot(data = df)
    plt <- plt + ggplot2::geom_ribbon(mapping = ggplot2::aes(x = n_subj, ymin = ci_lo, ymax = ci_up), fill = "#1e90ff", color = "#1e90ff", alpha = 0.33)
    plt <- plt + ggplot2::geom_line(mapping = ggplot2::aes(x = n_subj, y = de_ratio))
    plt <- plt + ggplot2::geom_point(mapping = ggplot2::aes(x = n_subj, y = de_ratio), size = 3)
    plt <- plt + ggplot2::geom_hline(yintercept = c(0.8, 1 / 0.8), color = "#ff0000", linetype = "dashed", size = 1)
    plt <- plt + ggplot2::xlab("Sample size") + ggplot2::ylab("Ratio")
    plt <- plt + ggplot2::theme_minimal()
    plt <- plt + ggplot2::theme(
      panel.background = ggplot2::element_rect(fill = "white", colour = "white"),
      plot.background = ggplot2::element_rect(fill = "white", colour = "white"),
      axis.text = ggplot2::element_text(size = 12), axis.title = ggplot2::element_text(size = 14),
      strip.text = ggplot2::element_text(size = 14)
    )


    if (addPlotHeaders) {
      plt <- plt + ggplot2::facet_grid(. ~ trial)
    }
    pltList[[trial]] <- plt
  }
  return(pltList)
}
