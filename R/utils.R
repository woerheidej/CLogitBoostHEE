#' Fit an offset boosting model for stratified/matched data
#' Checks if the two variables lead to a constant interaction.
#'
#' @param x First column
#' @param z Second column.
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

#' Helper: detect continuous variables
#' @param df data frame
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
#' @param data_proc Data frame containing the preprocessed dataset
#' @param exposure Name of exposure variable (character)
#' @param response Name of response variable (character)
#' @param strata Name of strata variable (character)
#' @param folds Binary matrix of complementary folds (rows = obs, cols = folds)
#'
#' @return Character vector of column names that are constant in at least one fold
#'
detect_singular_cols <- function(data_proc, exposure, response,outcome, strata, folds) {

  vars_to_check <- setdiff(names(data_proc), c(response, strata,outcome))

  const_by_pair <- sapply(vars_to_check, function(v) {
    # Count number of folds where variable or its interaction is constant
    sum(sapply(seq_len(ncol(folds)), function(j) {
      idx1 <- folds[, j] == 1
      idx2 <- folds[, j] == 0

      x1 <- data_proc[[v]][idx1]
      x2 <- data_proc[[v]][idx2]
      x3 <- data_proc[[exposure]][idx1]
      x4 <- data_proc[[exposure]][idx2]

      v1_const <- length(unique(x1[!is.na(x1)])) <= 1
      v2_const <- length(unique(x2[!is.na(x2)])) <= 1

      v3_const <- interaction_const(x1, x3)
      v4_const <- interaction_const(x2, x4)

      sum((v1_const | v3_const), (v2_const | v4_const))
    }))
  })

  return(const_by_pair)
}
