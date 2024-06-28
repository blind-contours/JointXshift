#' Estimate the Outcome Mechanism
#'
#' @details Compute the outcome regression for the observed data, including
#'  with the shift imposed by the intervention. This returns the outcome
#'  regression for the observed data (at A) and under the counterfactual shift
#'  shift (at A + delta).
#'
#' @param exposures A \code{character} vector of exposures to be shifted.
#' @param covars A \code{character} vector covariates to adjust for.
#' @param deltas A \code{numeric} indicating the magnitude of the shift to be
#'  computed for the exposure \code{A}. This is passed to the internal
#'  \code{\link{shift_additive}} and is currently limited to additive shifts.
#' @param mu_learner Object containing a set of instantiated learners from the
#'  \pkg{sl3}, to be used in fitting an ensemble model.
#' @param av A \code{dataframe} of validation data specific to the fold
#' @param at A \code{dataframe} of training data specific to the fold
#' @param outcome_type Variable type of the outcome
#' @importFrom stats glm as.formula predict
#' @importFrom data.table as.data.table setnames copy set
#' @importFrom stringr str_detect
#' @importFrom assertthat assert_that
#' @export
#' @return A \code{data.table} with two columns, containing estimates of the
#'  outcome mechanism at the natural value of the exposure Q(A, W) and an
#'  upshift of the exposure Q(A + delta, W).

joint_stoch_shift_est_Q <- function(exposures,
                                    deltas,
                                    mu_learner,
                                    covars,
                                    data,
                                    outcome_type) {
  future::plan(future::sequential, gc = TRUE)

  # scale the outcome for logit transform
  if (outcome_type != "binary") {
    y_star_data <- scale_to_unit(vals = data$y)

    data$y <- y_star_data
  }

  data_shifted <- data

  for ( i in 1:length(exposures)) {
    v <- exposures[[i]]
    d <- deltas[[i]]
    data_shifted[[v]] <- data_shifted[[v]] + d
  }
    # need a data set with the exposure stochastically shifted DOWNWARDS A-delta
  sl <- Lrnr_sl$new(
    learners = mu_learner,
    metalearner = sl3::Lrnr_nnls$new()
  )

    task_noshift <- suppressMessages(sl3::sl3_Task$new(
      data = data,
      covariates = covars,
      outcome = "y",
      outcome_type = "quasibinomial"
    ))

    task_upshift <- suppressMessages(sl3::sl3_Task$new(
      data = data_shifted,
      covariates = covars,
      outcome = "y",
      outcome_type = "quasibinomial"
    ))


    sl_fit <- suppressMessages(sl$train(task_noshift))

    # fit new Super Learner to the natural (no shift) data and predict
    pred_no_shift <- bound_precision(sl_fit$predict(task_noshift))
    pred_shift <- bound_precision(sl_fit$predict(task_upshift))

    # create output data frame and return result
    out <- data.table::as.data.table(cbind(
      pred_no_shift,
      pred_shift
    ))

    data.table::setnames(out, c("noshift", "upshift"))


  return(out)
}
