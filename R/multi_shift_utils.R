#' Create Augmented Data for Multiple Types of Shifts
#'
#' This function creates augmented datasets for multiple types of shifts, including marginal shifts and all two-way interactions. Each type of shift is given a unique and ordinal intervention indicator.
#'
#' @param data A data.frame containing the original dataset.
#' @param deltas A list of shift values for each exposure.
#' @param exposures A list of exposure variables to be shifted.
#' @param covars A list of covariates to be used in the modeling.
#' @return A data.frame containing the augmented dataset with multiple types of shifts and a unique intervention indicator for each type.
#' @examples
#' \dontrun{
#' data <- data.frame(A1 = rnorm(100), A2 = rnorm(100), W = rnorm(100), Y = rnorm(100))
#' deltas <- list(0.5, -0.3)
#' exposures <- list("A1", "A2")
#' covars <- list("W")
#' augmented_data <- create_augmented_data_multiple(data, deltas, exposures, covars)
#' }
#' @export
create_augmented_data_multiple <- function(data, deltas, exposures, covars, interaction = FALSE) {
  n <- nrow(data)
  shift_types <- c("noshift")

  for (i in 1:length(exposures)) {
    shift_types <- c(shift_types, paste0("shift_", exposures[[i]]))
  }

  for (i in 1:(length(exposures) - 1)) {
    for (j in (i + 1):length(exposures)) {
      shift_types <- c(shift_types, paste0("shift_", exposures[[i]], "_", exposures[[j]]))
    }
  }

  augmented_data <- data
  augmented_data$intervention <- rep(0, times = n)
  augmented_data$shift_type <- rep("noshift", times = n)

  shift_data_list <- list()
  shift_data_list[[1]] <- augmented_data

  # Create augmented datasets for each type of shift
  shift_index <- 2
  intervention_id <- 1
  for (i in 1:length(exposures)) {
    v <- exposures[[i]]
    d <- deltas[[i]]
    shifted_data <- data
    shifted_data[[v]] <- shifted_data[[v]] + d
    shifted_data$intervention <- rep(intervention_id, times = n)
    shifted_data$shift_type <- rep(paste0("shift_", v), times = n)
    shift_data_list[[shift_index]] <- shifted_data
    shift_index <- shift_index + 1
    intervention_id <- intervention_id + 1
  }

  if(interaction == TRUE){

  for (i in 1:(length(exposures) - 1)) {
    for (j in (i + 1):length(exposures)) {
      v1 <- exposures[[i]]
      v2 <- exposures[[j]]
      d1 <- deltas[[i]]
      d2 <- deltas[[j]]
      shifted_data <- data
      shifted_data[[v1]] <- shifted_data[[v1]] + d1
      shifted_data[[v2]] <- shifted_data[[v2]] + d2
      shifted_data$intervention <- rep(intervention_id, times = n)
      shifted_data$shift_type <- rep(paste0("shift_", v1, "_", v2), times = n)
      shift_data_list[[shift_index]] <- shifted_data
      shift_index <- shift_index + 1
      intervention_id <- intervention_id + 1
      }
    }
  }

  augmented_data <- do.call(rbind, shift_data_list)
  return(augmented_data)
}

#' Estimate Density Ratios for Multiple Types of Shifts
#'
#' This function estimates density ratios for multiple types of shifts, including marginal shifts and all two-way interactions, using a categorical SuperLearner.
#'
#' @param data A data.frame containing the original dataset.
#' @param deltas A list of shift values for each exposure.
#' @param exposures A list of exposure variables to be shifted.
#' @param covars A list of covariates to be used in the modeling.
#' @param classifier A classifier object to be used in the SuperLearner.
#' @return A list of density ratios for each type of shift.
#' @examples
#' \dontrun{
#' data <- data.frame(A1 = rnorm(100), A2 = rnorm(100), W = rnorm(100), Y = rnorm(100))
#' deltas <- list(0.5, -0.3)
#' exposures <- list("A1", "A2")
#' covars <- list("W")
#' classifier <- sl3::Lrnr_glm$new()
#' density_ratios <- estimate_density_ratio(data, deltas, exposures, covars, classifier)
#' }
#' @export

estimate_density_ratio_multiple <- function(data, deltas, exposures, covars, classifier, interaction = FALSE) {
  augmented_data <- create_augmented_data_multiple(data, deltas, exposures, covars, interaction = FALSE)

  sl_task <- sl3::sl3_Task$new(
    data = augmented_data,
    outcome = "intervention",
    covariates = covars,
    outcome_type = "categorical"
  )

  discrete_sl_metalrn <- sl3::Lrnr_cv_selector$new(sl3::loss_loglik_multinomial)
  discrete_sl_multinomial <- Lrnr_sl$new(learners = classifier, metalearner = discrete_sl_metalrn)


  class_model <- suppressWarnings(suppressMessages(discrete_sl_multinomial$train(sl_task)))

  # Get predictions for each type of shift
  class_model_preds <- class_model$predict(sl_task)

  extracted_vectors <- lapply(class_model_preds, function(x) x[[1]])

  # Convert to data frame
  df_preds <- do.call(rbind, extracted_vectors)
  colnames(df_preds) <- unique(augmented_data$shift_type)

  # Extract density ratios
  density_ratios <- as.data.frame(df_preds)
  density_ratios$shift_type <- augmented_data$shift_type

  # Compute density ratios for each shift type
  shift_types <- names(density_ratios)[names(density_ratios) != "shift_type"]

  density_ratios_list <- list()
  for (shift_type in shift_types) {
    density_ratios_subset <- density_ratios[density_ratios$shift_type == shift_type,]
    density_ratios_list[[paste0("ratio_", shift_type)]] <- density_ratios_subset[[shift_type]] / density_ratios_subset$noshift
  }

  return(density_ratios_list)
}

