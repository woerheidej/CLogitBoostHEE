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
#' @param exposure Character string giving the exposure variable name. Can be
#'   set to `NULL` for detection of heterogeneous effects in general.
#' @param response Character string for the response variable used in
#'   boosting (default `"resp"`). Typically a two-column matrix with
#'   outcome and strata.
#' @param strata Character string giving the strata variable used for
#'   matched case-control design.
#' @param outcome Character string giving the outcome variable name
#'   in the original data.
#' @param matching Character vector of variable names used for matching.
#'   These are excluded from predictors.
#' @param q Number of variables to select in each boosting iteration.
#' @param PFER Per-family error rate for stability selection.
#' @param cutoff Selection probability cutoff for stability selection.
#'   Exactly two of `q`, `PFER`, `cutoff` must be specified.
#' @param nu Step size for boosting (learning rate, default 1).
#' @param mstop Maximum number of boosting iterations (default 1000).
#' @param B Number of complementary pairs for stability selection (default 50).
#' @param sampling_type Either `"SS"` for Shah & Samworth complementary pairs stability selection
#'  or `"MB"` for the original stability selection approach bei Meinshausen & Bühlmann
#'  (default `"SS"`).
#' @param assumption Error bound assumption for complementary pairs stability selection,
#'   one of `"none"`, `"unimodal"`, or `"rconcave"`.
#' @param n_cores Number of cores for parallel stability selection (default 1).
#' @param df_bols Degrees of freedom for linear base learners.
#' @param df_bbs Degrees of freedom for smooth base learners.
#' @param intercept Logical, should base learners include an intercept (default FALSE).
#' @param center Logical, should continuous variables be centered (default TRUE).
#' @param flexible Logical, if TRUE both linear and
#'   smooth base learners are included for continuous covariates.
#' @param reduction_scaler Number used for scaling of the reduced stabsel mstop using
#' the formula \code{q \* 5 \* (1/nu) \* reduction}.
#' @param early_stopping Logical, should the offset model be refined through
#' early stopping? Can be deactivated for high-dimensional data to save runtime. (default TRUE)
#'
#' @param remove Logical, if TRUE, removes variables that are singular in any fold.
#'
#' @return A `stabsel` object containing selected variables and selection
#'   probabilities.
#'
#' @examples
#' \dontrun{
#' data_sim <- create_data() # simulation data
#' stab_model <- CLogitBoostingHEE(data_sim$data,
#'                          exposure = "X",
#'                          response = "resp",
#'                          strata = "strata",
#'                          outcome = "y",
#'                          q = 5, PFER = 0.1, cutoff = NULL)
#' }
#'
#' @import stats
#' @import mboost
#' @import stabs

#' @export
CLogitBoostingHEE <- function(
    data,
    exposure = NULL,
    response = "resp",
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
    assumption = "none",
    n_cores = 1,
    df_bols = 1,
    df_bbs = 1,
    intercept = FALSE,
    center = TRUE,
    flexible = TRUE,
    reduction_scaler = 1,
    early_stopping = TRUE,
    remove = FALSE
) {

  # Stability selection parameter check:
  provided <- c(q = !is.null(q), PFER = !is.null(PFER), cutoff = !is.null(cutoff))
  n_provided <- sum(provided)
  if (n_provided != 2) {
    stop(sprintf(
      "Exactly two of 'q', 'PFER', and 'cutoff' must be specified, but you provided %d:\n  q = %s, PFER = %s, cutoff = %s",
      n_provided,
      ifelse(is.null(q), "NULL", q),
      ifelse(is.null(PFER), "NULL", PFER),
      ifelse(is.null(cutoff), "NULL", cutoff)
    ), call. = FALSE)
  }

  # Helper: detect continuous variables
  detect_continuous <- function(df, exclude = c()) {
    vars <- setdiff(names(df), exclude)
    is_cont <- sapply(df[, vars, drop = FALSE], function(x)
      is.numeric(x) && length(unique(x)) > 2
    )
    list(
      cont_vars = vars[is_cont],
      cat_vars  = vars[!is_cont]
    )
  }

  # Detect variable types
  vars_info <- detect_continuous(data, exclude = c(response, strata, outcome))
  cont_vars <- vars_info$cont_vars
  cat_vars  <- vars_info$cat_vars

  # Preprocess data
  data_proc <- data
  data_proc[cat_vars] <- lapply(data_proc[cat_vars], factor)
  data_proc[cont_vars] <- scale(data_proc[cont_vars])
  # Create response variable:
  data_proc$resp <- cbind(data_proc[[outcome]], data_proc[[strata]])

  # Generate offset model and predictions
  offset_formula <- generate_formula(
    data = data_proc,
    exposure = exposure,
    response = response,
    strata = strata,
    df_bols = df_bols,
    df_bbs = df_bbs,
    outcome = outcome,
    matching = matching,
    flexible = flexible,
    include_interactions = FALSE
  )
  offset.cv <- gen_offset_model(data = data_proc, formula = offset_formula$form, mstop = mstop, nu = nu, strata = strata, n_cores = n_cores, early_stopping = early_stopping)
  offset_pred <- predict(offset.cv, type = "link")

  # Create stratified folds
  strata_vec <- data_proc[[strata]]
  folds <- matrix(0, ncol = B, nrow = nrow(data_proc))
  for (j in 1:B) {
    smp <- sample(unique(strata_vec), floor(length(unique(strata_vec)) / 2))
    folds[strata_vec %in% smp, j] <- 1
  }

  # Check for constant variables in complementary folds
  # Check for constant columns and interactions fold-wise
  const_by_pair <- sapply(setdiff(names(data_proc), c(response, strata)), function(v) {
    # Count number of folds where the variable (or its interaction) is constant
    sum(sapply(seq_len(ncol(folds)), function(j) {
      idx1 <- folds[, j] == 1
      idx2 <- folds[, j] == 0

      x1 <- data_proc[[v]][idx1]
      x2 <- data_proc[[v]][idx2]
      x3 <- data_proc[[exposure]][idx1]
      x4 <- data_proc[[exposure]][idx2]

      v1_const <- length(unique(x1[!is.na(x1)])) <= 1
      v2_const <- length(unique(x2[!is.na(x2)])) <= 1

      # interaction check (safe for numeric/factor)
      v3_const <- length(unique(interaction(x1[!is.na(x1)], x3[!is.na(x3)], drop = TRUE))) <= 1
      v4_const <- length(unique(interaction(x2[!is.na(x2)], x4[!is.na(x4)], drop = TRUE))) <= 1

      sum((v1_const | v3_const),(v2_const|v4_const))
    }))
  })

  singular_cols <- names(const_by_pair[const_by_pair > 0])

  if (length(singular_cols) > 0) {
    # Generate informative message
    singular_msg <- paste0(
      singular_cols,
      " (", const_by_pair[singular_cols], " fold", ifelse(const_by_pair[singular_cols] > 1, "s", ""), ")"
    )

    if(remove){
      warning(sprintf(
        "Some columns are constant (singular) in at least one fold: %s. They are removed from analysis.",
        paste(singular_msg, collapse = ", ")
      ), call. = FALSE)
      data_proc <- data_proc %>%
        select(-singular_cols)
    }
    else{
    warning(sprintf(
      "Some columns are constant (singular) in at least one fold: %s. If that poses a problem, consider removing or adjusting these variables before rerunning.",
      paste(singular_msg, collapse = ", ")
    ), call. = FALSE)}
  }




  # Main formula
  main_formula <- generate_formula(
    data = data_proc,
    exposure = exposure,
    response = response,
    strata = strata,
    df_bols = df_bols,
    df_bbs = df_bbs,
    outcome = outcome,
    matching = matching,
    flexible = flexible,
    include_interactions = TRUE
  )


  # The start of stability selection part:
  mstop_reduced <- q * 5 * (1/nu) * reduction_scaler # reduce to gain efficiency in computation

  # Fit initial boosting model
  initial_model <- gamboost(
    main_formula$form,
    data = data_proc,
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
  if (!is.null(q))      stabsel_args$q      <- q
  if (!is.null(PFER))   stabsel_args$PFER   <- PFER
  if (!is.null(cutoff)) stabsel_args$cutoff <- cutoff

  stabsel_model <- tryCatch({
    withCallingHandlers(
      do.call(stabsel, stabsel_args),
      warning = function(w) {
        # supress repetitive warnings
        if (grepl("Lapack routine dgesv|Original error message", conditionMessage(w))) {
          invokeRestart("muffleWarning") # suppress warnings
        }
      }
    )
  }, error = function(e) {
    warning("Stability selection failed: ", conditionMessage(e), call. = FALSE)
    NULL
  })
  return(stabsel_model)
}
