#' @title globalVariableNames
#' @description
#' Remove check warning for the variables listed below
#' Variables called withing dplyr and tidyr functions
#' @keywords internal
globalVariableNames <- c(
  "ci_lo", "ci_up", "cofactorPaths", "de_ratio", "drug",
  "id", "n_subj", "outputValues", "pBE_rep", "per", 
  "variable", "ymax", "ymean", "ymin"
)

utils::globalVariables(c(globalVariableNames))
