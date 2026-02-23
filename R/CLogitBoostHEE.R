#' Detect Heterogeneous Exposure Effects via Boosting and Stability Selection
#'
#' Fits a component-wise gradient boosting model with optional
#' interactions for exposure effect heterogeneity and performs
#' stability selection on the fitted model.
#'
#' This function automates preprocessing (scaling continuous covariates
#' and converting categorical covariates to factors), generates the
#' boosting formula, computes offsets if needed, fits a `gamboost`
#' model, and applies stability selection using `stabsel`.
#'
#' @param data A `data.frame` containing the outcome, exposure, covariates,
#'   strata, and matching variables.
#'
#' @param exposure Character string giving the exposure variable name.
#'   Can be set to `NULL` for detection of heterogeneous effects in general.
#'
#' @param strata Character string giving the strata variable used for a
#'   matched case–control design.
#'
#' @param outcome Character string giving the outcome variable name in
#'   the original data.
#'
#' @param matching Character vector of variable names used for matching.
#'   These variables are excluded from the predictor set.
#'
#' @param q Numeric giving the number of variables to select in each
#'   boosting iteration.
#'
#' @param PFER Numeric giving the per-family error rate for stability
#'   selection.
#'
#' @param cutoff Numeric giving the selection probability cutoff for
#'   stability selection. Exactly two of `q`, `PFER`, and `cutoff`
#'   must be specified.
#'
#' @param nu Numeric giving the boosting step size (learning rate,
#'   default `1`).
#'
#' @param mstop Integer giving the maximum number of boosting iterations
#'   (default `2000`).
#'
#' @param B Integer giving the number of complementary pairs used for
#'   stability selection (default `50`).
#'
#' @param sampling_type Character string specifying the stability selection
#'   scheme: `"SS"` for Shah & Samworth complementary pairs stability
#'   selection or `"MB"` for the original Meinshausen & Bühlmann approach
#'   (default `"SS"`).
#'
#' @param assumption Character string specifying the error-bound assumption
#'   used for complementary pairs stability selection. One of `"none"`,
#'   `"unimodal"`, or `"r-concave"` (default `r-concave`).
#'
#' @param n_cores Integer giving the number of CPU cores used for parallel
#'   stability selection (default `1`). Be careful when choosing more than one core
#'   in scenarios with multiple thousands of observations or hundreds of parameters
#'   due to RAM bottleneck.
#'
#' @param df_bols Integer giving the degrees of freedom for linear base
#'   learners.
#'
#' @param df_bbs Integer giving the degrees of freedom for smooth base
#'   learners.
#'
#' @param intercept Logical indicating whether base learners should include
#'   an intercept (default `FALSE`).
#'
#' @param center Logical indicating whether continuous variables should be
#'   centered during preprocessing (default `TRUE`).
#'
#' @param flexible Logical indicating whether both linear and smooth base
#'   learners are included for continuous covariates (default `TRUE`).
#'
#' @param reduction_scaler Numeric used to scale the reduced `mstop` in
#'   stability selection according to
#'   \code{q * 10 * (1 / nu) * reduction_scaler}.
#'
#' @param early_stopping Logical indicating whether the offset model should
#'   be refined via early stopping. Can be disabled for high-dimensional
#'   data to reduce runtime (default `TRUE`).
#'
#' @return A `stabsel` object containing the selected variables and their
#'   selection probabilities.
#'
#' @examples
#' \dontrun{
#' data(sim)
#' sim_results <- CLogitBoostHEE(
#'   sim$data,
#'   exposure = "X",
#'   strata = "strata",
#'   outcome = "y",
#'   matching = "s",
#'   q = 5,
#'   PFER = 0.1,
#'   cutoff = NULL
#' )
#' }
#'
#' @import stats
#' @import mboost
#' @import stabs
#' @import RhpcBLASctl
#'
#' @export

CLogitBoostHEE <- function(data,
                              exposure = NULL,
                              strata = "strata",
                              outcome = "y",
                              matching = NULL,
                              q = NULL,
                              PFER = NULL,
                              cutoff = NULL,
                              nu = 1,
                              mstop = 2000,
                              B = 50,
                              sampling_type = "SS",
                              assumption = "r-concave",
                              n_cores = 1,
                              df_bols = 1,
                              df_bbs = 1,
                              intercept = FALSE,
                              center = TRUE,
                              flexible = TRUE,
                              reduction_scaler = 1,
                              early_stopping = TRUE,
                              only_boosting = FALSE,
                              boosting_interactions = NULL) {

  if(only_boosting == FALSE){ # Stability selection parameter check:
  provided <- c(
    q = !is.null(q),
    PFER = !is.null(PFER),
    cutoff = !is.null(cutoff)
  )
  n_provided <- sum(provided)
  if (n_provided != 2) {
    stop(
      sprintf(
        "Exactly two of 'q', 'PFER', and 'cutoff' must be specified, but you provided %d:\n  q = %s, PFER = %s, cutoff = %s",
        n_provided,
        ifelse(is.null(q), "NULL", q),
        ifelse(is.null(PFER), "NULL", PFER),
        ifelse(is.null(cutoff), "NULL", cutoff)
      ),
      call. = FALSE
    )
  }}

  # Detect variable types
  vars_info <- detect_continuous(data, exclude = c(strata, outcome))
  cont_vars <- vars_info$cont_vars
  cat_vars  <- vars_info$cat_vars

  # Preprocess data
  data[cat_vars] <- lapply(data[cat_vars], factor)
  data[cont_vars] <- scale(data[cont_vars])
  # Create response variable:
  data$resp <- cbind(data[[outcome]], data[[strata]])

  # Generate offset model and predictions
  offset_formula <- generate_formula(
    data = data,
    exposure = exposure,
    response = "resp",
    strata = strata,
    df_bols = df_bols,
    df_bbs = df_bbs,
    outcome = outcome,
    matching = matching,
    flexible = flexible,
    include_interactions = FALSE,
    boosting_interactions = boosting_interactions
  )
  offset.cv <- gen_offset_model(
    data = data,
    formula = offset_formula$form,
    mstop = mstop,
    nu = nu,
    strata = strata,
    n_cores = n_cores,
    early_stopping = early_stopping
  )
  if(only_boosting){
    return(offset.cv)
  }
  offset_pred <- predict(offset.cv, type = "link")

  # Create stratified folds
  strata_vec <- data[[strata]]
  folds <- matrix(0, ncol = B, nrow = nrow(data))
  for (j in 1:B) {
    smp <- sample(unique(strata_vec), floor(length(unique(strata_vec)) / 2))
    folds[strata_vec %in% smp, j] <- 1
  }

  singularity <- detect_singular_cols(
    data = data,
    exposure = exposure,
    outcome = outcome,
    response = "resp",
    strata = strata,
    folds = folds
  )
  singular_cols <- names(singularity[singularity > 0])

  if (length(singular_cols) > 0) {
    # Generate informative message
    singular_msg <- paste0(singular_cols,
                           " (",
                           singularity[singular_cols],
                           " fold",
                           ifelse(singularity[singular_cols] > 1, "s", ""),
                           ")")
    warning(
      sprintf(
        "Some columns are constant (singular) in at least one fold: %s. They are removed from HEE analysis.",
        paste(singular_msg, collapse = ", ")
      ),
      call. = FALSE
    )

  }



  # Main formula
  main_formula <- generate_formula(
    data = data,
    exposure = exposure,
    response = "resp",
    strata = strata,
    df_bols = df_bols,
    df_bbs = df_bbs,
    outcome = outcome,
    matching = matching,
    flexible = flexible,
    singular = singular_cols,
    include_interactions = TRUE
  )


  # The start of stability selection part:
  mstop_reduced <- q * 10 * (1 / nu) * reduction_scaler # reduce to gain efficiency in computation

  # Fit initial boosting model
  initial_model <- gamboost(
    main_formula$form,
    data = data,
    family = CLogit(),
    control = boost_control(mstop = mstop_reduced, nu = nu),
    offset = offset_pred
  )

  stabsel_args <- list(
    initial_model,
    folds = folds,
    sampling.type = sampling_type,
    assumption = assumption,
    B = B,
    mc.cores = n_cores
  )

  if (!is.null(q))
    stabsel_args$q      <- q
  if (!is.null(PFER))
    stabsel_args$PFER   <- PFER
  if (!is.null(cutoff))
    stabsel_args$cutoff <- cutoff

  RhpcBLASctl::blas_set_num_threads(1)

  if (.Platform$OS.type == "windows" && n_cores > 1) {
    cores <- n_cores

    # Run with the chunk.size argument
    cl <- parallel::makeCluster(cores)
    on.exit(parallel::stopCluster(cl), add = TRUE)

    parallel::clusterEvalQ(cl, {
      library(mboost)
    })

    myApply <- function(X, FUN, ...) {
      parallel::parLapply(cl, X, FUN, ...)
    }

    # Switch from mc.cores to papply
    stabsel_args$mc.cores <- NULL
    stabsel_args$papply   <- myApply
  }

  stabsel_model <- tryCatch({
    withCallingHandlers(
      do.call(stabsel, stabsel_args),
      warning = function(w) {
        if (grepl("Lapack routine dgesv|Original error message",
                  conditionMessage(w))) {
          invokeRestart("muffleWarning")
        }
      }
    )
  }, error = function(e) {
    warning("Stability selection failed: ", conditionMessage(e), call. = FALSE)
    NULL
  })

  stabsel_model$call <- NULL

  return(stabsel_model)
}
