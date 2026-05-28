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
#' @param early_stopping Logical. Toggle early stopping on or off.
#' @param steady_state_percentage Integer Threshold for minimal improvement in CV risk
#'   to declare the model in steady state (default 0.01, i.e., 0.01% change).
#' @param n_cores Integer number of how many cores are available for CV.
#' @param do_plot Logical. If the CV of the offset should be printed.
#'
#' @return A fitted `gamboost` object with optimal number of iterations.
#'
#' @examples
#' \dontrun{
#' offset_mod <- gen_offset_model(data, formula = "resp ~ X + Z1 + Z2",
#' mstop = 500, nu = 0.1, strata = "strata")
#' }
#' @import parallel
#' @export
gen_offset_model <- function(data, formula, mstop, nu, strata,
                             n_cores = 1, K = 5,
                             early_stopping = TRUE,
                             steady_state_percentage = 0.01,
                             do_plot = TRUE) {

  RhpcBLASctl::blas_set_num_threads(n_cores)

  offset_model <- gamboost(
    formula,
    data = data,
    family = CLogit(),
    control = boost_control(mstop = mstop, nu = nu)
  )

  coefs_initial <- coef(offset_model)

  if (!early_stopping) return(offset_model)

  RhpcBLASctl::blas_set_num_threads(1)
  sim.folds <- make_cv_folds(data, strata, K = K)

  if (.Platform$OS.type == "windows" && n_cores > 1) {
    cores <- min(n_cores, K)
    cl <- parallel::makeCluster(cores)

    parallel::clusterEvalQ(cl, library(mboost))

    myApply <- function(X, FUN, ...) parallel::parLapply(cl, X, FUN, ...)

    cv_stopping <- cvrisk(offset_model, folds = sim.folds, papply = myApply)

    parallel::stopCluster(cl)
  } else {
    cores <- min(n_cores, K)
    cv_stopping <- cvrisk(offset_model, folds = sim.folds, mc.cores = cores)
  }

  opt <- mstop(cv_stopping)

  # check non-finite
  vals <- as.numeric(cv_stopping)
  has_nonfinite <- any(!is.finite(vals))

  if (has_nonfinite) {
    message("Non-finite CV risk values detected; diagnostics saved to cv_infinity.RData")
    save(cv_stopping, opt, coefs_initial, file = "cv_infinity.RData")
  }

  # Robust mean CV risk per iteration (column), ignoring non-finite
  mean_risk <- apply(cv_stopping, 2, function(x) mean(x[is.finite(x)], na.rm = TRUE))

  # Steady-state check on last 5 iterations
  last_i <- length(mean_risk)
  idx <- max(1, last_i - 4): (last_i) # last 5 iterations
  rolling_change <- mean(sapply(idx, function(i) {
    (mean_risk[i-1] / mean_risk[i] - 1)
  }), na.rm = TRUE) * 100

  is_steady <- rolling_change < steady_state_percentage

  if (opt >= (mstop - 5)) {
    if (is_steady) {
      message("Optimal mstop near max, but CV risk looks steady -> probably ok.")
    } else {
      message("Optimal mstop near max and not steady -> consider increasing mstop or reducing nu.")
    }
  }

  if (isTRUE(do_plot)) {
    if (!has_nonfinite) {
      plot(cv_stopping, main = "CV Early Stopping")
    } else {
      # optional: plot mean finite risks instead of failing
      plot(mean_risk, type = "l", main = "CV Early Stopping (finite mean risk)",
           xlab = "mstop", ylab = "mean CV risk (finite only)")
    }
  }

  return(offset_model[opt])
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
  if (remainder > 0) {
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
