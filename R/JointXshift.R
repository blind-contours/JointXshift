#' @title JointXshift
#'
#' @description Under a fixed shift to exposures identify using g-computation the joint shift of
#' pairwise exposures in a mixed exposure compared to the additive individual shifts. Positive values indicate
#' synergy and negative antagonism, get the top synergy and antagonism results and use CV-TMLE to efficiently
#' estimate the interaction target parameter.
#'
#' @param w A \code{matrix}, \code{data.frame}, or similar containing a set of
#' baseline covariates. These variables are measured before exposures.
#' @param a \code{matrix}, \code{data.frame}, or similar containing individual or
#' multiple exposures.
#' @param z \code{matrix}, \code{data.frame}, or similar containing individual or
#' multiple mediators (optional).
#' @param y \code{numeric} vector of observed outcomes.
#' @param deltas A \code{numeric} value indicating the shift in exposures to
#' define the target parameter, with respect to the scale of the exposures (A). If adaptive_delta
#' is true, these values will be reduced.
#' @param var_sets A list specifying variable sets for deterministic JointXshift usage.
#' Example: var_sets <- c("A_1", "A_1-Z_2") where the analyst provides variable sets
#' for exposures, exposure-mediator, or exposure-covariate relationships.
#' @param estimator The type of estimator to fit: \code{"tmle"} for targeted
#' maximum likelihood estimation, or \code{"onestep"} for a one-step estimator.
#' @param fluctuation Method used in the targeting step for TML estimation: "standard" or "weighted".
#' This determines where to place the auxiliary covariate in the logistic tilting regression.
#' @param pi_learner Learners for fitting Super Learner ensembles to densities via \pkg{sl3}.
#' @param mu_learner Learners for fitting Super Learner ensembles to the outcome model via \pkg{sl3}.
#' @param g_learner Learners for fitting Super Learner ensembles to the g-mechanism
#' g(A|W) (a probability estimator, not a density estimator) for mediation via \pkg{sl3}.
#' @param e_learner Learners for fitting Super Learner ensembles to the e-mechanism
#' g(A|Z,W) (a probability estimator, not a density estimator) for mediation via \pkg{sl3}.
#' @param zeta_learner Learners for fitting Super Learner ensembles to the outcome model via \pkg{sl3}..
#' @param n_folds Number of folds to use in cross-validation, default is 2.
#' @param outcome_type Data type of the outcome, default is "continuous".
#' @param parallel Whether to parallelize across cores (default: TRUE).
#' @param parallel_type Type of parallelization to use if parallel is TRUE:
#' "multi_session" (default), "multicore", or "sequential".
#' @param num_cores Number of CPU cores to use in parallelization (default: 2).
#' @param seed \code{numeric} seed value to be passed to all functions.
#' @param hn_trunc_thresh Truncation level for the clever covariate (default: 10).
#' @param adaptive_delta If TRUE, reduces the user-specified delta until
#' the Hn calculated for a shift does not have any observation greater
#' than hn_trunc_thresh (default: FALSE).
#'
#' @return An S3 object of class \code{JointXshift} containing the results of the
#' procedure to compute a TML or one-step estimate of the counterfactual mean
#' under a modified treatment policy that shifts a continuous-valued exposure
#' by a scalar amount \code{delta}. These exposures are data-adaptively
#' identified using the CV-TMLE procedure.
#' @export
#' @importFrom MASS mvrnorm
#' @importFrom foreach %dopar%
#' @importFrom magrittr %>%
#' @importFrom stats as.formula glm p.adjust plogis predict qlogis qnorm qunif rnorm runif
#' @importFrom rlang :=
#' @importFrom stringr str_count
#' @import furrr
#' @importFrom purrr map
#' @importFrom data.table rbindlist

JointXshift <- function(w,
                      a,
                      y,
                      deltas,
                      mu_learner = NULL,
                      g_learner = NULL,
                      n_folds = 2,
                      outcome_type = "continuous",
                      parallel = TRUE,
                      parallel_type = "multi_session",
                      num_cores = 2,
                      seed = seed,
                      hn_trunc_thresh = 50) {

  # coerce W to matrix and, if no names in W, assign them generically
  if (!is.data.frame(w)) w <- as.data.frame(w)
  w_names <- colnames(w)
  if (is.null(w_names)) {
    w_names <- paste0("w", seq_len(ncol(w)))
    colnames(w) <- w_names
  }


  # coerce W to matrix and, if no names in W, assign them generically
  a <- data.frame(a)
  a_names <- colnames(a)

  if (is.null(a_names)) {
    a_names <- paste0("a", seq_len(ncol(a)))
    colnames(a) <- a_names
  }


  if (is.null(mu_learner)) {
    sls <- create_sls()
    mu_learner <- sls$mu_learner
  }

  if (is.null(g_learner)) {
    sls <- create_sls()
    g_learner <- sls$g_learner
  }



  if (parallel == TRUE) {
    if (parallel_type == "multi_session") {
      future::plan(future::multisession,
                   workers = num_cores,
                   gc = TRUE
      )
    } else {
      future::plan(future::multicore,
                   workers = num_cores,
                   gc = TRUE
      )
    }
  } else {
    future::plan(future::sequential,
                 gc = TRUE
    )
  }

  data_internal <- data.table::data.table(w, a, y)
  `%notin%` <- Negate(`%in%`)

  joint_Hn <- estimate_density_ratio(data = data_internal, deltas =  deltas, exposures = a_names, covars = c(a_names, w_names), classifier = g_learner)

  joint_Hn$shift <- pmin(joint_Hn$shift, hn_trunc_thresh)

  covars <- c(a_names, w_names)

  joint_qn_estims <- joint_stoch_shift_est_Q(
    exposures = a_names,
    deltas = deltas,
    mu_learner = mu_learner,
    covars = covars,
    data = data_internal,
    outcome_type = outcome_type
  )

  joint_tmle_fit <- tmle_exposhift(
    data_internal = data_internal,
    delta = mean(unlist(deltas)),
    Qn_scaled = joint_qn_estims,
    Hn = joint_Hn,
    y = data_internal$y
  )


  joint_out <- calc_joint_results(joint_tmle_fit)

  multi_Hn <- estimate_density_ratio_multiple(data = data_internal, deltas =  deltas, exposures = a_names, covars = c(a_names, w_names), classifier = mu_learner, interaction = FALSE)

  indiv_stoch_shift_est_Q

  indiv_qn_estims <- indiv_stoch_shift_est_Q(
    exposures = a_names,
    deltas = deltas,
    mu_learner = mu_learner,
    covars = covars,
    data = data_internal,
    outcome_type = outcome_type
  )

  tmle_debias_est_list <- list()
  tmle_ice_list <- list()

  for (i in 1:length(a_names)) {

    indiv_qn_estim_i <- indiv_qn_estims[[i]]
    indiv_hn_estim_i <- multi_Hn[[i]]

    hn_compare <- as.data.frame(cbind(multi_Hn[[1]], indiv_hn_estim_i))
    colnames(hn_compare) <- c("noshift", "shift")

    tmle_fit <- tmle_exposhift(
      data_internal = data_internal,
      delta = deltas[[i]],
      Qn_scaled = indiv_qn_estim_i,
      Hn = hn_compare,
      y = data_internal$y
    )

    tmle_debias_est_list[i] <- tmle_fit$psi - tmle_fit$noshift_psi
    tmle_ice_list[[i]] <- tmle_fit$eif - tmle_fit$noshift_eif
  }

  additive_psi <- sum(unlist(tmle_debias_est_list))
  additive_eif  <- Reduce(`+`, tmle_ice_list)

  joint_vs_additive_psi <- (joint_tmle_fit$psi - joint_tmle_fit$noshift_psi) - additive_psi
  joint_vs_additive_eif <- (joint_tmle_fit$eif - joint_tmle_fit$noshift_eif) - additive_eif

  psi_var <- var(joint_vs_additive_eif) /
    length(joint_vs_additive_eif)

  se_ests <- sqrt(psi_var)
  CI <- calc_CIs(joint_vs_additive_psi, se_ests)
  p_vals <- calc_pvals(joint_vs_additive_psi, se_ests)

  # Create a data frame with relevant column names
  joint_vs_additive_results <- data.frame(
    psi = joint_vs_additive_psi,
    psi_var = psi_var,
    se_ests = se_ests,
    CI_lower = CI[1],
    CI_upper = CI[2],
    p_vals = p_vals
  )

  psi_additive_var <- var(additive_eif) /
    length(additive_eif)
  additive_se_ests <- sqrt(psi_additive_var)
  additive_CI <- calc_CIs(additive_psi, additive_se_ests)
  p_vals_additive <- calc_pvals(additive_psi, additive_se_ests)

  # Create a data frame with relevant column names
  additive_results <- data.frame(
    psi = additive_psi,
    psi_var = psi_additive_var,
    se_ests = additive_se_ests,
    CI_lower = additive_CI[1],
    CI_upper = additive_CI[2],
    p_vals = p_vals_additive
  )

  out <- list("Joint Effects" = joint_out, "Additive Effects" = additive_results, "Joint vs Additive Effects" = joint_vs_additive_results)


  return(out)
}
