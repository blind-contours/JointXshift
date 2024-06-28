#' Create Augmented Data for Stochastic Interventions
#'
#' This function creates augmented data for both training and validation datasets by duplicating observations
#' and applying a specified shift (delta) to the exposure variable. The function is useful for implementing
#' stochastic interventions in causal inference studies.
#'
#' @param at Data frame containing the training data.
#' @param av Data frame containing the validation data.
#' @param delta Numeric value specifying the shift to be applied to the exposure variable.
#' @param var Character string specifying the name of the exposure variable to be shifted.
#' @param covars Character vector specifying the names of the covariate columns in the data frames.
#'
#' @return A list containing two elements: \code{at_dup} and \code{av_dup}, which are the augmented training
#' and validation data frames, respectively. Each data frame has an additional column \code{intervention}
#' indicating whether the observation is under the original or shifted exposure.
#'
#' @examples
#' # Example usage:
#' training_data <- data.frame(A = rnorm(100), W1 = rbinom(100, 1, 0.5), W2 = rbinom(100, 1, 0.5))
#' validation_data <- data.frame(A = rnorm(50), W1 = rbinom(50, 1, 0.5), W2 = rbinom(50, 1, 0.5))
#' delta <- -0.5
#' exposure_variable <- "A"
#' covariates <- c("W1", "W2")
#' result <- create_augmented_data(training_data, validation_data, delta, exposure_variable, covariates)
#' augmented_training_data <- result$at_dup
#' augmented_validation_data <- result$av_dup
#'
#' @export
create_augmented_data <- function(data, deltas, exposures, covars) {
  n <- nrow(data)

  data_dup <- data
  data_dup$intervention <- rep(1, times = n)
  data$intervention <- rep(0, times = n)

  for ( i in 1:length(exposures)) {
    v <- exposures[[i]]
    d <- deltas[[i]]
    data_dup[[v]] <- data_dup[[v]] + d
  }

  data_dup <- rbind(data, data_dup)

  return(data_dup)
}

###############################################################################
#' Estimate Density Ratio via Classification
#'
#' This function estimates the density ratio of shifted versus unshifted exposures using a classification approach.
#' It utilizes the `sl3` package to train a Super Learner on the intervention indicator and then predicts the density ratios
#' for both the training and validation datasets.
#'
#' @param at Data frame containing the training data.
#' @param av Data frame containing the validation data.
#' @param delta Numeric value specifying the shift to be applied to the exposure variable.
#' @param var Character string specifying the name of the exposure variable to be shifted.
#' @param covars Character vector specifying the names of the covariate columns in the data frames.
#' @param classifier An `sl3` learner object used for classification.
#'
#' @return A list containing two elements: \code{Hn_at} and \code{Hn_av}, which are data tables of density ratios for the
#' training and validation datasets, respectively. Each data table has two columns: \code{noshift} and \code{shift},
#' representing the density ratios for unshifted and shifted exposures.
#'
#' @examples
#' # Example usage:
#' training_data <- data.frame(A = rnorm(100), W1 = rbinom(100, 1, 0.5), W2 = rbinom(100, 1, 0.5))
#' validation_data <- data.frame(A = rnorm(50), W1 = rbinom(50, 1, 0.5), W2 = rbinom(50, 1, 0.5))
#' delta <- -0.5
#' exposure_variable <- "A"
#' covariates <- c("W1", "W2")
#' classifier <- sl3::Lrnr_glm$new()
#' result <- estimate_density_ratio(training_data, validation_data, delta, exposure_variable, covariates, classifier)
#' training_density_ratios <- result$Hn_at
#' validation_density_ratios <- result$Hn_av
#'
#' @export
estimate_density_ratio <- function(data, deltas, exposures, covars, classifier) {
  augmented_data <- create_augmented_data(data, deltas, exposures, covars)

  sl_task <- sl3::sl3_Task$new(
    data = augmented_data,
    outcome = "intervention",
    covariates = covars,
    outcome_type = "binary"
  )

  sl <- sl3::Lrnr_sl$new(
    learners = classifier,
    metalearner = sl3::Lrnr_nnls$new()
  )


  class_model <- suppressWarnings(suppressMessages(sl$train(sl_task)))

  # at predictions -----------
  class_model_preds <- class_model$predict(sl_task)

  # Compute the density ratios for shifted exposures
  data_u_t_shift <- class_model_preds[augmented_data$intervention == 1]
  data_density_ratio_shift <- data_u_t_shift / (1 - data_u_t_shift)

  # Compute the density ratios for unshifted exposures
  data_density_ratio_unshift <- 1

  # Combine the density ratios into a data.table
  density_ratios <- data.table::as.data.table(
    cbind(data_density_ratio_unshift, data_density_ratio_shift)
  )

  data.table::setnames(density_ratios, c("noshift", "shift"))


  return(density_ratios)
}
###############################################################################
#' @title Calculate the Joint Parameter
#' @description Using the output results for the eif use the delta method
#' to simply subtract the joint estimate from the no shift estimate. Do
#' the same for the eif to get the variance for this new parameter.
#' @param results_element Joint results for formatting
#' @export

calc_joint_results <- function(results_element) {
  psi <- results_element$psi - results_element$noshift_psi
  psi_var <- var(results_element$eif - results_element$noshift_eif) /
    length(results_element$eif)
  se_ests <- sqrt(psi_var)
  CI <- calc_CIs(psi, se_ests)
  p_vals <- calc_pvals(psi, se_ests)

  # Create a data frame with relevant column names
  joint_results <- data.frame(
    psi = psi,
    psi_var = psi_var,
    se_ests = se_ests,
    CI_lower = CI[1],
    CI_upper = CI[2],
    p_vals = p_vals
  )

  return(joint_results)
}


###############################################################################
#' @title Calculates the Interaction Parameter
#' @description Using the output results for the eif use the delta method
#' to simply subtract the two additive estimates from the joint shift. Do
#' the same for the eif to get the variance for this new parameter.
#' @param results_table Table of interaction results
#' @param joint_shift_fold_results Joint shift results
#' @export

calc_intxn_results <- function(results_table, joint_shift_fold_results) {
  intxn_psi <- results_table[[3, 1]] -
    results_table[[2, 1]] -
    results_table[[1, 1]]

  psi_variance <- var((joint_shift_fold_results[[3]]$eif -
    joint_shift_fold_results[[3]]$noshift_eif) -
    (joint_shift_fold_results[[2]]$eif -
      joint_shift_fold_results[[2]]$noshift_eif) -
    (joint_shift_fold_results[[1]]$eif -
      joint_shift_fold_results[[1]]$noshift_eif)) /
    length(joint_shift_fold_results[[1]]$eif)

  psi_se <- sqrt(psi_variance)

  CI <- calc_CIs(intxn_psi, psi_se)
  p_vals <- calc_pvals(intxn_psi, psi_se)

  intxn_results <- c(intxn_psi, psi_variance, psi_se, CI[1], CI[2], p_vals)

  return(intxn_results)
}

###############################################################################
#' @title Calculate Confidence Interals
#' @description Gives information on which variable sets are used in the
#' basis function in the data-adaptive estimation section
#' @param psi The Psi parameter calculated
#' @param psi_se Standard deviation for the Psi parameter
#' @export
calc_CIs <- function(psi, psi_se) {
  psi_CIs <- c(
    round(psi + stats::qnorm(0.05 / 2, lower.tail = T) * psi_se, 4),
    round(psi + stats::qnorm(0.05 / 2, lower.tail = F) * psi_se, 4)
  )
  return(psi_CIs)
}


###############################################################################
#' @title Calculate the frequency basis functions are used in the folds
#' @description Gives information on which variable sets are used in the
#' basis function in the data-adaptive estimation section
#' @param fold_basis_results List of lists of variable sets found in the basis
#' function data-adaptive parameter section
#' @param n_folds Number of folds used
#'
#' @export
calc_basis_freq <- function(fold_basis_results, n_folds) {
  fold_basis_table <- table(unlist(fold_basis_results))
  fold_basis_table <- round(fold_basis_table / n_folds, 2)

  return(fold_basis_table)
}



###############################################################################
#' @title Extract variables from the basis function search
#' @description This simple function creates a list of target variables
#' that were used in the data adpative procedure.
#' @param common_variables Variable sets found in the data adaptive procedure
#' @param i iteration of search
#' @param a_names Names of the exposures
#' @param w_names Names of the covariates
#' @param z_names Mediator names
#'
#' @export
extract_vars_from_basis <- function(common_variables, i, a_names,
                                    w_names, z_names) {
  target <- common_variables[i]
  match_list <- list()

  if (length(unlist(strsplit(target, "-"))) == 2) {
    intxn_vars <- unlist(strsplit(target, "-"))
  }
  for (j in 1:length(c(a_names, z_names, w_names))) {
    var <- c(a_names, z_names, w_names)[j]
    matches <- var == target
    if (matches) {
      match_list[j] <- var
    }
  }

  matches <- unlist(match_list[!match_list %in% ""])
  return(list("matches" = matches, "target" = target))
}

###############################################################################
#' @title Calculate the p-values based on variance of the influence curve
#' @description Calculates the p-value using the standard formula
#' @param psi The estimated target parameter
#' @param variance The variance of that parameter
#'
#' @export
calc_pvals <- function(psi, variance) {
  p_value <- 2 * stats::pnorm(abs(psi / sqrt(variance)), lower.tail = F)
  return(p_value)
}

is.JointXshift <- function(x) {
  class(x) == "JointXshift"
}

is.JointXshift <- function(x) {
  class(x) == "JointXshift_msm"
}
