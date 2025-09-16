getLmerParallel <- function(this_lmer, summary_settings) {
  summary_settings <- list(ci = 0.9, ci_lcut = 0.8, ci_ucut = 1.25)

  ci <- summary_settings$ci
  ci_lcut <- summary_settings$ci_lcut
  ci_ucut <- summary_settings$ci_ucut

  sum_lm <- summary(this_lmer)
  n_subj <- length(this_lmer$residuals) / length(this_lmer$xlevels$drug)  ##WHAT IF NOT INTEGER?
  dd_eff <- sum_lm$coefficients[2,1] #drug effect
  se <- sum_lm$coefficients[2,2] #standard error
  df <- sum_lm$df[2] #degrees of freedom

  de_ratio <- exp(dd_eff)
  deltaCI <- stats::qt((1 - ci) / 2, df, lower.tail = FALSE) * se
  ci_lo <- de_ratio * exp(-deltaCI)
  ci_up <- de_ratio * exp(deltaCI)
  if ((ci_lo >= ci_lcut) & (ci_up <= ci_ucut)) {
    pBE_rep <- 1
  } else {
    pBE_rep <- 0
  }

  return(data.frame(
    de_ratio = de_ratio, ci_lo = ci_lo, ci_up = ci_up,
    beta = dd_eff, se = se, df = df, deltaCI = deltaCI,
    pBE_rep = pBE_rep, n_subj = n_subj
  ))
}

getLmerCrossover <- function(this_lmer, summary_settings) {
  summary_settings <- list(ci = 0.9, ci_lcut = 0.8, ci_ucut = 1.25)

  ci <- summary_settings$ci
  ci_lcut <- summary_settings$ci_lcut
  ci_ucut <- summary_settings$ci_ucut

  sum_lm <- summary(this_lmer)
  n_subj <-  length(unique(this_lmer$groups$id))

  dd_eff <- sum_lm$tTable[2, 1] # drug effect
  se <- sum_lm$tTable[2, 2] # standard error
  df <- sum_lm$fixDF$X[2] # degrees of freedom

  de_ratio <- exp(dd_eff)
  deltaCI <- stats::qt((1 - ci) / 2, df, lower.tail = FALSE) * se
  ci_lo <- de_ratio * exp(-deltaCI)
  ci_up <- de_ratio * exp(deltaCI)
  if ((ci_lo >= ci_lcut) & (ci_up <= ci_ucut)) {
    pBE_rep <- 1
  } else {
    pBE_rep <- 0
  }

  mse <- sum_lm$sigma^2
  cv_intra <- sqrt(exp(mse) - 1)

  return(data.frame(
    de_ratio = de_ratio, ci_lo = ci_lo, ci_up = ci_up,
    beta = dd_eff, se = se, df = df, deltaCI = deltaCI,
    cv_intra = cv_intra, pBE_rep = pBE_rep, n_subj = n_subj
  ))
}

getLmerReplicate <- getLmerCrossover

getLmerFunctionsList <- list(parallel = getLmerParallel, crossover = getLmerCrossover, replicate = getLmerReplicate)

get_this_lmer <- function(this_lmer, design, summary_settings) {
  results_summary <- getLmerFunctionsList[[design]](this_lmer, summary_settings)
  return(results_summary)
}

#' @title get_summary
#' @param lmer_res The output of (`calc_BE()`).
#' @param summary_settings A list of settings for the summary function.
#' @param design The design of the study, either `"parallel"`, `"crossover"`, or `"replicate"`.
#' @param index The index of the design, either 1 or 2.
#' @return A tibble with the summary statistics.
#' @author Michael Neely
#' @author Abdullah Hamadeh
#' @export
#' @importFrom purrr map map_dfr
#' @importFrom dplyr group_by as_tibble
get_summary <- function(lmer_res, summary_settings, design, index) {
  designIndex <- list(parallel = 1, crossover = 2, replicate = 3)
  this_res <- purrr::flatten(lmer_res)
  this_res <- purrr::map(this_res, c(designIndex[[design]], index))
  this_res <- purrr::map(this_res, get_this_lmer, designIndex[[design]], summary_settings)
  this_res <- purrr::map_dfr(this_res, ~ dplyr::as_tibble(.))
  this_res <- dplyr::group_by(this_res, n_subj)
  return(this_res)
}


#' Extracts BE statistics from results of (`calc_BE()`).
#'
#' @title Extract  BE statistics
#' @param lmer_res The output of (`calc_BE()`).
#' @param ci Confidence interval for geometric mean ratio, default is 0.9 (90%)
#' per FDA guidance.
#' @param ci_lcut Lower confidence interval cutoff for geometric mean AUC or Cmax to be
#' considered BE, default 0.8 per FDA guidance.
#' @param ci_ucut Upper confidence interval cutoff for geometric mean AUC or Cmax to be
#' considered BE, default `1.25` per FDA guidance.
#' @return A list of lists: `[[1:n_samp]][[1:n_trials]][["par", "cross", "rep"]][[1,2]]`,
#' where `n_samp` is the number of sample sizes tested, e.g., (`subj_max` - `subj-min`) / `subj_step`.
#' The first item in the design is auc, and the second is cmax.
#' For example the comparison of auc using the third trial of the second sample size
#' for the replicated design would be `x[[2]][[3]]$rep[[1]]`.
#' @author Michael Neely
#' @author Abdullah Hamadeh
#' @export
#' @importFrom dplyr summarize across matches
extract_be <- function(lmer_res, ci = 0.9, ci_lcut = 0.8, ci_ucut = 1.25) {
  summary_settings <- list(ci = 0.9, ci_lcut = 0.8, ci_ucut = 1.25)

  n_samp <- length(lmer_res)
  n_trials <- length(lmer_res[[1]])
  rep_design <- !is.null(lmer_res[[1]][[1]]$rep[[1]]) #inherits(lmer_res[[1]][[1]]$rep[[1]], "lmerModLmerTest")

  auc_par_data <- get_summary(lmer_res = lmer_res, summary_settings = summary_settings, design = "parallel", 1)
  auc_par <- dplyr::summarize(
    auc_par_data,
    dplyr::across(.cols = !dplyr::matches("pBE_rep"), mean),
    pBE_rep = sum(pBE_rep) / n_trials
  )

  cmax_par_data <- get_summary(lmer_res = lmer_res, summary_settings = summary_settings, design = "parallel", 2)
  cmax_par <- dplyr::summarize(
    cmax_par_data,
    dplyr::across(.cols = !dplyr::matches("pBE_rep"), mean),
    pBE_rep = sum(pBE_rep) / n_trials
  )

  auc_cross_data <- get_summary(lmer_res = lmer_res, summary_settings = summary_settings, design = "crossover", 1)
  auc_cross <- dplyr::summarize(
    auc_cross_data,
    dplyr::across(.cols = !dplyr::matches("pBE_rep"), mean),
    pBE_rep = sum(pBE_rep) / n_trials
  )

  cmax_cross_data <- get_summary(lmer_res = lmer_res, summary_settings = summary_settings, design = "crossover", 2)
  cmax_cross <- dplyr::summarize(
    cmax_cross_data,
    dplyr::across(.cols = !dplyr::matches("pBE_rep"), mean),
    pBE_rep = sum(pBE_rep) / n_trials
  )


  auc_par_all <- get_summary(lmer_res = lmer_res, summary_settings = summary_settings, design = "parallel", 1)
  cmax_par_all <- get_summary(lmer_res = lmer_res, summary_settings = summary_settings, design = "parallel", 2)
  auc_cross_all <- get_summary(lmer_res = lmer_res, summary_settings = summary_settings, design = "crossover", 1)
  cmax_cross_all <- get_summary(lmer_res = lmer_res, summary_settings = summary_settings, design = "crossover", 2)

  if (rep_design) {
    auc_rep_all <- get_summary(lmer_res = lmer_res, summary_settings = summary_settings, design = "replicate", 1)
    auc_rep <- dplyr::summarize(
      auc_rep_all,
      dplyr::across(.cols = !dplyr::matches("pBE_rep"), mean),
      pBE_rep = sum(pBE_rep) / n_trials
    )

    cmax_rep_all <- get_summary(lmer_res = lmer_res, summary_settings = summary_settings, design = "replicate", 2)
    cmax_rep <- dplyr::summarize(
      cmax_rep_all,
      dplyr::across(.cols = !dplyr::matches("pBE_rep"), mean),
      pBE_rep = sum(pBE_rep) / n_trials
    )
  } else {
    auc_rep <- NA
    cmax_rep <- NA
    auc_rep_all <- NA
    cmax_rep_all <- NA
  }

  return(list(
    summary = list(
      auc_par = auc_par, cmax_par = cmax_par,
      auc_cross = auc_cross, cmax_cross = cmax_cross,
      auc_rep = auc_rep, cmax_rep = cmax_rep
    ),
    trials = list(
      auc_par_all = auc_par_all, cmax_par_all = cmax_par_all,
      auc_cross_all = auc_cross_all, cmax_cross_all = cmax_cross_all,
      auc_rep_all = auc_rep_all, cmax_rep_all = cmax_rep_all
    )
  ))
}
