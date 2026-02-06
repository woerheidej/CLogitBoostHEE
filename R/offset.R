#' Fit an offset boosting model for stratified/matched data
#'
#' Fits a component-wise gradient boosting model to estimate an offset
#' for the main effect model. Uses cross-validation to determine
#' the optimal number of boosting iterations (`mstop`).
#'
#' @param data A data.frame containing the outcome, covariates, and strata variable.
#' @param formula A `gamboost` formula defining predictors for the offset model.
#' @param mstop Integer maximum number of boosting iterations.
#' @param nu Numeric step size for boosting (learning rate, default 1).
#' @param strata Character string naming the strata variable for matched design.
#' @param K Integer of folds for cross-validation (default 5).
#' @param steady_state_percentage Integer Threshold for minimal improvement in CV risk
#'   to declare the model in steady state (default 0.01, i.e., 0.01% change).
#' @param n_cores Integer number of how many cores are available for CV.
#' @param plot Logical. If the CV of the offset should be printed.
#'
#' @return A fitted `gamboost` object with optimal number of iterations.
#'
#' @examples
#' \dontrun{
#' offset_mod <- gen_offset_model(data, formula = "resp ~ X + Z1 + Z2", mstop = 500, nu = 0.1, strata = "strata")
#' }
#' @export
gen_offset_model <- function(data, formula, mstop, nu, strata, n_cores = 1, K = 5, early_stopping = TRUE, steady_state_percentage = 0.01, plot = TRUE) {

  cores <- min(n_cores, K)
  # Fit initial boosting model
  offset_model <- gamboost(
    formula,
    data = data,
    family = CLogit(),
    control = boost_control(mstop = mstop, nu = nu)
  )

  if(early_stopping){
    # Generate CV folds for strata
    sim.folds <- make_cv_folds(data, strata, K = K)
    # Cross-validation
    early_stopping <- cvrisk(offset_model, folds = sim.folds, mc.cores = cores)
    opt <- mstop(early_stopping)
  # Check if model has reached steady state
  # (Checks the last 5 iterations and calculates a mean rolling change)
  range_idx <- max(1, mstop - 4):mstop
  rolling_change <- mean(
    sapply(range_idx, function(i) {
      (mean(early_stopping[, i]) / mean(early_stopping[, i+1]) - 1)
    })
  ) * 100
  is_steady <- rolling_change < steady_state_percentage

  # Message if optimal mstop is close to or equal maximum
  if (opt %in% (mstop-5):mstop) {
    if(is_steady){
      message("Optimal CV stopping iteration is near maximum, but the model is nearly stable. Should be fine!")
    } else{
      message("Optimal CV stopping iteration is near or equal to the maximum. Consider increasing mstop or adjusting nu.")
       }
  }
  if(plot){
    plot(early_stopping)
  }

  # Return model refitted to optimal mstop
  return(offset_model[opt])  }
  else{ return(offset_model)}
}


#' Create stratified cross-validation folds for matched data
#'
#' Generates fold assignments for cross-validation in matched or stratified
#' designs, keeping all observations within a stratum together.
#'
#' @param data A data.frame containing the strata variable.
#' @param strata Character string of the strata/matching variable.
#' @param K Number of folds for cross-validation (default 10).
#'
#' @return A matrix of fold indicators (0 = held-out, 1 = training set).
#'
#' @examples
#' \dontrun{
#' folds <- make_cv_folds(data, strata = "strata", K = 5)
#' }
#' @export
#'
#'
make_cv_folds <- function(data, strata, K = 10) {

  strata.unique <- sample(unique(data[[strata]]))
  n.strata <- length(strata.unique)

  # Compute number of strata per fold
  n.fold <- rep(floor(n.strata / K), K)
  remainder <- n.strata - sum(n.fold)
  if(remainder > 0){
    n.fold[seq_len(remainder)] <- n.fold[seq_len(remainder)] + 1
  }

  # Initialize fold matrix
  folds <- matrix(1, nrow = nrow(data), ncol = K)

  start <- 0
  for (i in seq_len(K)) {
    strata.i <- strata.unique[(start + 1):(start + n.fold[i])]
    folds[data[[strata]] %in% strata.i, i] <- 0
    start <- start + n.fold[i]
  }

  folds
}
