#' @title compute_auc
#' @description Compute AUC and Cmax
#' @param data A data.frame or tibble with these columns: time, conc.
#' @import stats
#' @import utils
compute_auc <- function(data) {
  # Ensure the dataframe has the required columns
  if (!all(c("time", "conc") %in% names(data))) {
    stop("Data must contain 'time' and 'conc' columns.")
  }

  # Sort data by time to ensure correct integration
  data <- data[order(data$time), ]

  # Compute AUC using the trapezoidal rule
  auc <- sum(diff(data$time) * (utils::head(data$conc, -1) + utils::tail(data$conc, -1)) / 2)

  return(auc)
}

compute_cmax <- function(data) {
  # Ensure the dataframe has the required columns
  if (!all(c("time", "conc") %in% names(data))) {
    stop("Data must contain 'time' and 'conc' columns.")
  }

  # Compute Cmax
  cmax <- max(data$conc)

  return(cmax)
}

getAucCmaxDf <- function(df) {
  categoryNames <- c("id","drug","rep")
  #categoryNames <- c("id","drug","rep","seq")
  #  results <- data.frame(id = numeric(), drug = character(), rep = numeric(), seq = character(), auc = numeric(), cmax = numeric())
  results <- data.frame(id = numeric(), drug = character(), rep = numeric(), auc = numeric(), cmax = numeric())
  uniqueTimeProfiles <- unique(df[,categoryNames])

  for (nr in seq_len(nrow(uniqueTimeProfiles))){
    drug <- uniqueTimeProfiles$drug[nr]
    id <- uniqueTimeProfiles$id[nr]
    rep <- uniqueTimeProfiles$rep[nr]
    #seq <- uniqueTimeProfiles$seq[nr]

    #subset_data <- df[df$id == id & df$drug == drug & df$rep == rep & df$seq == seq, ]
    subset_data <- df[df$id == id & df$drug == drug & df$rep == rep, ]

    auc_value <- compute_auc(subset_data)
    cmax_value <- compute_cmax(subset_data)
    results <- rbind.data.frame(results,
                                data.frame(id = id,
                                           drug = drug,
                                           rep = rep,
                                           #seq = seq,
                                           auc = auc_value,
                                           cmax = cmax_value))
  }



  return(results)
}


get_lmer <- function(nSubjects, pkParameterData, n_trials, aucICV, cmaxICV) {
  res <- vector("list", n_trials)  # Preallocate list

  #loop through each of n_trials
  for (i in seq_len(n_trials)) {
    res[[i]] <- this_samp(pkParameterData, nSubjects, aucICV, cmaxICV)
  }
  return(res)
}

get_parallel_design_data <- function(pkDf,nSubjects){
  par_df <- NULL
  drugNames <- unique(pkDf$drug)
  for (drugName in drugNames) {
    drugPkDf <- pkDf[pkDf$drug == drugName & pkDf$rep == "1", ]
    idsThisDrug <- sample(x = unique(drugPkDf$id), size = nSubjects, replace = FALSE)
    par_df <- rbind.data.frame(par_df, drugPkDf[drugPkDf$id %in% idsThisDrug, ])
    pkDf <- pkDf[!(pkDf$id %in% idsThisDrug), ]
  }
  return(par_df)
}

addICV <- function(pkParameterData,aucICV,cmaxICV){
  sigmaAUC <- sqrt( log(1 + aucICV^2 ))
  sigmaCmax <- sqrt( log(1 + cmaxICV^2) )
  pkParameterData$auc <- exp( log(pkParameterData$auc) + rnorm( n = length(pkParameterData$auc) , mean = 0 , sd =  sigmaAUC) )
  pkParameterData$cmax <- exp( log(pkParameterData$cmax) + rnorm( n = length(pkParameterData$cmax) , mean = 0 , sd =  sigmaCmax) )
  return(pkParameterData)
}

checkIdAppearsOncePerReplicate <- function(pkDf){
  pkDfSummary <- as.data.frame(table(pkDf$id, pkDf$drug, pkDf$rep))
  colnames(pkDfSummary) <- c("id", "drug", "rep", "N")
  if(!all(pkDfSummary$N == 1)){
    repeatedIds <- pkDfSummary$id[pkDfSummary$N != 1]
    stop(paste("Subject IDs", repeatedIds ,"appear more than once."))
  }
}

# summarizePerSeqDrug <- function(pkDf){
#   pkDfSummary <- as.data.frame(table(pkDf$per, pkDf$drug, pkDf$seq))
#   colnames(pkDfSummary) <- c("per", "drug" , "seq", "N")
#   pkDfSummary <- pkDfSummary[pkDfSummary$N != 0, ]
#   return(pkDfSummary)
# }
#
# checkIdsSeqDisjoint <- function(pkDf){
#   pkDfSummary <- summarizePerSeqDrug(pkDf)
#   idsDf <- NULL
#   for (seq in unique(pkDfSummary$seq)){
#     ids <- unique(pkDf$id[pkDf$seq == seq])
#     idsDf <- rbind.data.frame(idsDf,
#                               data.frame(id = ids,
#                                          seq = seq))
#   }
#   duplicatedIds <- duplicated(idsDf$id)
#   if( any(duplicatedIds)){
#     stop( paste("Ids",paste0(idsDf$id[duplicatedIds],collapse = ", "),"present in multiple sequences."))
#   }
# }


# summarizePerDrug <- function(pkDf){
#   pkDfSummary <- as.data.frame(table(pkDf$per, pkDf$drug))
#   colnames(pkDfSummary) <- c("per", "drug" , "N")
#   pkDfSummary <- pkDfSummary[pkDfSummary$N != 0, ]
#   return(pkDfSummary)
# }

# checkIds <- function(pkDf){
#   #want to say that the ids appearing in per 1 for drug R are the same as the ids appearing in per 2 for drug T1 and in per 3 for drug T4
#   pkDfSummary <- summarizePerSeqDrug(pkDf)
#   browser()
#   for (seq in unique(pkDfSummary$seq)){
#     df <- pkDf[pkDf$seq == seq,]
#     #df$
#   }
# }

get_replicate_design_data <- function(pkDf,nSubjects, aucICV, cmaxICV){
  # create subset with the same nSubjects random individuals by drug, period
  # for replicate (all periods) and cross-over (period 1 only)


  # pkDfSummary <- summarizePerSeqDrug(pkDf)
  checkIdAppearsOncePerReplicate(pkDf)
  # checkIdsSeqDisjoint(pkDf)

  rep_df <- NULL

  #drugNames <- c("R","T1","T2","T3")
  drugNames <- unique(pkDf$drug)
  nDrugs <- length(drugNames)
  nPeriodsPerReplicate <- nDrugs
  periods <- seq_len(nPeriodsPerReplicate)
  names(periods) <- drugNames

  # #Williams
  # sequencesDf <- expand.grid(replicate(nPeriodsPerReplicate, periods, simplify = FALSE))
  # validSequences <- sapply(seq_len(nrow(sequencesDf)),function(n){ all(periods %in% unlist(sequencesDf[n,]) ) })
  # sequencesDf <- sequencesDf[validSequences,]
  # names(sequencesDf) <- drugNames
  # seqOorder <- order( sapply(seq_len(nrow(sequencesDf)) , function(n){ paste(unlist(sequencesDf[n,]),collapse = "") } ))
  # sequencesDf <- sequencesDf[seqOorder , ]

  drugSequenceDf <- data.frame()
  drugSequence <- seq_len(nDrugs)
  for (drug in drugNames){
    drugSequenceDf <- rbind(drugSequenceDf,data.frame(t(drugSequence)))
    drugSequence <- c(tail(drugSequence, 1), head(drugSequence, -1))
  }
  names(drugSequenceDf) <- drugNames

  uniqueIds <- unique(pkDf$id)
  idsPerSequence <- list()
  for (seq in seq_len(nrow(drugSequenceDf))){
    idsPerSequence[[seq]] <- sample(x = uniqueIds, size = ceiling(nSubjects/nDrugs) )
    uniqueIds <- uniqueIds[ !(uniqueIds %in% idsPerSequence[[seq]]) ]
  }


  for (rep in unique(pkDf$rep)){
    for (seq in seq_len(nrow(drugSequenceDf))){
      for (drug in names(drugSequenceDf)){
        per <- drugSequenceDf[[drug]][seq]
        df <- pkDf[ pkDf$id %in% idsPerSequence[[seq]] & pkDf$drug == drug & pkDf$rep == rep ,]
        df$per <- per + ((rep-1)*nPeriodsPerReplicate)
        df$seq <- seq
        rep_df <- rbind.data.frame(rep_df,df)
      }
    }
  }

  rep_df <-  addICV(pkParameterData = rep_df,
                    aucICV = aucICV,
                    cmaxICV = cmaxICV)

  return(rep_df)
  #check that ids in only appear once per period
  #pkdf will have per col.
  #pkdf will have all subjs taking either R or T in each period.  ie a subject cannot appear more than once in the same period
  #pkdf will
  #user specifies how many periods.  virtual pop generator will add IOV to parameters, ie in each period the simulation will be run with a different set of IOV parameters.
  #crossover has two periods where subjects switch between R and T
  #replicate has more than one
  #Could have multiple treatments
}

get_crossover_design_data <- function(pkDf, aucICV, cmaxICV){
  cross_df <-  addICV(pkParameterData = pkDf,
                      aucICV = aucICV,
                      cmaxICV = cmaxICV)
  return(cross_df)
}

#' @title this_samp
#' @description
#' create subset with n random individuals by drug for parallel design
#' @param pkParameterData The pharmacokinetic parameter data.
#' @param nSubjects The number of subjects to sample.
#' @keywords internal
#' @importFrom nlme lme
#' @importFrom lmerTest lmer
#' @importFrom lme4 lmerControl
this_samp <- function(pkParameterData, nSubjects, aucICV, cmaxICV) {

  par_df <- get_parallel_design_data(pkDf = pkParameterData ,nSubjects =  nSubjects)
  rep_df <- get_replicate_design_data(pkDf = pkParameterData ,nSubjects =  nSubjects, aucICV, cmaxICV)
  cross_df <- rep_df[rep_df$rep %in% 1, ] #get_crossover_design_data(pkDf = rep_df[rep_df$rep %in% 1, ], aucICV, cmaxICV)

  lmPar <- list(rep(NA, 2)) # for parallel
  lmCross <- list(rep(NA, 2)) # for crossover
  lmerRep <- list(rep(NA, 2)) # for replicate

  lmPar[[1]] <- suppressWarnings(lm(log(auc) ~ drug, data = par_df))
  lmPar[[2]] <- suppressWarnings(lm(log(cmax) ~ drug, data = par_df))

  ctrl <- ctrl <- nlme::lmeControl(opt = "optim", msMaxIter = 200)


  lmCross[[1]] <- suppressWarnings(nlme::lme(fixed = log(auc) ~ drug + per + seq, random = ~ 1 | id, data = cross_df,control = ctrl))
  lmCross[[2]] <- suppressWarnings(nlme::lme(fixed = log(cmax) ~ drug + per + seq, random = ~ 1 | id, data = cross_df,control = ctrl))

  replicates <- unique(pkParameterData$rep)
  n_reps <- length(replicates)
  if (n_reps > 1) {
    lmerRep[[1]] <- suppressWarnings(nlme::lme(fixed = log(auc) ~ drug + per + seq, random = ~ 1 | id, data = rep_df,control = ctrl))
    lmerRep[[2]] <- suppressWarnings(nlme::lme(fixed = log(cmax) ~ drug + per + seq, random = ~ 1 | id, data = rep_df,control = ctrl))
    return(list(par = lmPar, cross = lmCross, rep = lmerRep))
  }
  return(list(par = lmPar, cross = lmCross))
}




#' Generates bioequivalence statistics
#'
#' Generates BE statistics for parallel, cross-over and full replication designs.
#' From simulated data with 1000 subjects per drug, the function first calculates
#' area under the concentration-time curve (AUC) and maximum concentration (Cmax)
#' for each simulated profiles. Periods are defined as >1 for each replicated
#' administration (e.g., when intra-occasion varibility is simulated.)
#' Linear regression (`stats::lm()`) is used on `log(auc)` and`log(cmax)`
#' for parallel and cross-over designs, and mixed effect modeling (`lmerTest::lmer()`)
#' is used for the replicated design.
#'
#' @title Calculation of BE statistics
#' @param data A dataframe or tibble  with these columns:
#' subject id, time, concentrations, period, drug.
#' @param n_trials Integer coding the number of random trials per sample size.
#' @param subj_min Integer defining the minimum number of subjects per trial
#' @param subj_max Integer defining the maximum number of subjects per trial
#' @param subj_step Integer defining the step size between `subj_min` and
#' `subj_max`, so that the program will test sample sizes
#' $N = \{subj_min, subj_min + subj_step, subj_min + 2*subj_step,...,subj_max\}$.
#' @param seed Random number seed.
#' @param ci_lcut Lower limit of the confidence interval.
#' @param ci_ucut Upper limit of the confidence interval.
#' @return A list of lists: `[[1:n_samp]][[1:n_trials]][["par", "cross", "rep"]][[1,2]]`,
#' where `n_samp` is the number of sample sizes tested, e.g.,
#' (`subj_max` - `subj-min`) / `subj_step`.  The first item in the
#' design is auc, and the second is cmax.  For example the comparison of auc using
#' the third trial of the second sample size for the replicated design would be
#' `x[[2]][[3]]$rep[[1]]`.
#' @author Michael Neely
#' @author Abdullah Hamadeh
#' @export
calc_be <- function(data,
                    n_trials = 50,
                    subj_min = 5, subj_max = 50, subj_step = 5,
                    seed = -17, ci_lcut = .8, ci_ucut = 1.25,
                    aucICV = 0,#0.147,
                    cmaxICV = 0 #0.217
){

  pkParameterData <- getAucCmaxDf(data)
  # extract data

  periods <- unique(pkParameterData$per)
  # n_pers <- length(periods) # each period is a replicate for that drug
  nSubjects <- seq(from = subj_min, to = subj_max, by = subj_step) # step through sample sizes
  #n_samp <- length(nSubjects)
  set.seed(seed)

  res_list <- vector("list", length(nSubjects))

  #loop through each trial size
  for (ns in seq_along(nSubjects)) {
    print(nSubjects[ns])
    res_list[[ns]] <- get_lmer(nSubjects[[ns]], pkParameterData, n_trials, aucICV = aucICV, cmaxICV = cmaxICV)
  }

  return(res_list)
}
