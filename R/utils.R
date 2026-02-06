#' Helper: detect continuous variables
#' @param df A data.frame containing all data.
#'
#' @return Two lists, continuous variables, and categorical ones.
#'
#' @export
#'
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

#' Detect constant or singular columns (including interactions) per fold
#'
#' Checks each variable (and its interaction with exposure) for being constant
#' within complementary folds. Works for numeric, factor, and mixed variables.
#'
#' @param data_proc Data frame containing the preprocessed dataset.
#' @param exposure Character name of exposure variable.
#' @param response Character name of response variable.
#' @param strata Character name of strata variable.
#' @param folds Binary matrix of complementary folds (rows = obs, cols = folds).
#' @param sampling_type Character string specifying the stability selection
#'   scheme: `"SS"` for Shah & Samworth complementary pairs stability
#'   selection or `"MB"` for the original Meinshausen & Bühlmann approach
#'   (default `"SS"`).
#'
#' @return Character vector of column names that are constant in at least one fold.
#'
detect_singular_cols <- function(data_proc, exposure, response,outcome, strata, folds, sampling_type = "SS") {

  vars_to_check <- setdiff(names(data_proc), c(response, strata,outcome))

  const_by_pair <- sapply(vars_to_check, function(v) {
    # Count number of folds where variable or its interaction is constant
    sum(sapply(seq_len(ncol(folds)), function(j) {
      idx1 <- folds[, j] == 1
      x1 <- data_proc[[v]][idx1]
      x3 <- data_proc[[exposure]][idx1]

      # Check whether the variable is constant in the fold
      v1_const <- length(unique(x1[!is.na(x1)])) <= 1
      # Check whether the interaction with the exposure is constant in the fold
      v3_const <- interaction_const(x1, x3)
      # Default to v2_const & v4_const FALSE to comply with the logic
      v2_const <- FALSE
      v4_const <- FALSE

      # Do the same for the complementary pair if appropriate
      if(sampling_type == "SS"){
      idx2 <- folds[, j] == 0
      x2 <- data_proc[[v]][idx2]
      x4 <- data_proc[[exposure]][idx2]
      v2_const <- length(unique(x2[!is.na(x2)])) <= 1
      v4_const <- interaction_const(x2, x4)}

      # 0 = not constant, 1 = constant in one fold, 2 = constant in both folds
      sum((v1_const | v3_const), (v2_const | v4_const))
    }))
  })

  return(const_by_pair)
}

#' Fit an offset boosting model for stratified/matched data
#' Checks if the two variables lead to a constant interaction.
#'
#' @param x Character name of first column.
#' @param z Character name of second column.
#'
#' @return TRUE/FALSE, if the two variables lead to constant interaction.
#'
#' @export
#'
interaction_const <- function(x, z) {
  ok <- !is.na(x) & !is.na(z)
  if (!any(ok)) return(TRUE)

  x <- x[ok]
  z <- z[ok]

  # factor × factor
  if (is.factor(x) && is.factor(z)) {
    ref_x <- levels(x)[1]
    ref_z <- levels(z)[1]
    return(!any(x != ref_x & z != ref_z))
  }

  # numeric × factor
  if (is.numeric(x) && is.factor(z)) {
    ref_z <- levels(z)[1]
    return(!any(x != 0 & z != ref_z))
  }

  if (is.factor(x) && is.numeric(z)) {
    ref_x <- levels(x)[1]
    return(!any(x != ref_x & z != 0))
  }

  # numeric × numeric
  if (is.numeric(x) && is.numeric(z)) {
    prod <- x * z
    return(length(unique(prod)) <= 1)
  }

  stop("Unsupported type combination")
}

#' Plot Stability Selection Results with Cleaned Variable Names
#'
#' Plots a `stabsel` object while cleaning variable names by removing
#' technical annotations such as `df = 1`, `intercept = FALSE/TRUE`,
#' and `center = TRUE/FALSE`. The plot retains the original `stabsel`
#' gradient bars and layout, with a customizable title.
#'
#' @param sim_results A `stabsel` object returned from a stability
#'   selection procedure (e.g., the output of `CLogitBoostingHEE`).
#'
#' @param main_title Character string specifying the plot title. Default
#'   is `"Results of Stability Selection"`.
#'
#' @return Invisibly returns `NULL`. Generates a plot of the stability
#'   selection results.
#'
#' @examples
#' \dontrun{
#' # Suppose stab_model is a stabsel object from your analysis
#' plot_results(stab_model)
#' plot_results(stab_model, main_title = "HEE Selection Probabilities")
#' }
#'
#' @import stabs
#' @export
plot_results <- function(sim_results, main_title = "Results of Stability Selection") {

  # Check that sim_results is a stabsel object
  if (!inherits(sim_results, "stabsel")) {
    stop("sim_results must be a 'stabsel' object.")
  }

  # Clean variable names:
  new_names <- names(sim_results$max)

  # Remove df = <number>
  new_names <- gsub("df = [0-9]+", "", new_names)

  # Remove center = TRUE/FALSE
  new_names <- gsub("center = (TRUE|FALSE)", "", new_names)

  # Remove intercept = TRUE/FALSE
  new_names <- gsub("intercept = (TRUE|FALSE)", "", new_names)

  # Remove commas and extra =
  new_names <- gsub(",", "", new_names)
  new_names <- gsub(" =", "", new_names)

  # Remove trailing double spaces before closing parentheses
  new_names <- gsub("  \\)", ")", new_names)

  # Trim leading/trailing whitespace
  new_names <- trimws(new_names)

  # Assign cleaned names back
  names(sim_results$max) <- new_names

  # Plot with custom title and adjusted y-axis label
  plot(sim_results, main = main_title, xlab = "Selection Frequency")
  invisible(NULL)
}



