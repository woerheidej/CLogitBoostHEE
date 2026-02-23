#' Construct a \code{gamboost} formula with main effects and interactions
#'
#' Programmatically builds a model formula for component-wise gradient
#' boosting using linear (\code{bols}) and optional smooth (\code{bbs})
#' base learners. Optionally includes exposure–covariate interaction
#' terms to enable detection of heterogeneous exposure effects.
#'
#' @param data A data.frame containing the response, exposure, strata,
#'   matching variables, and covariates.
#'
#' @param exposure Character string giving the name of the exposure
#'   variable. If supplied and \code{include_interactions = TRUE},
#'   interactions between the exposure and all other covariates are added. If
#'   \code{NULL}, interactions between all covariates are added.
#'
#' @param response Character string giving the name of the response
#'   tuple used on the left-hand side of the formula (default: \code{"resp"}).
#'
#' @param strata Character string giving the name of the strata variable.
#'   This variable is excluded from the set of predictors.
#'
#' @param outcome Character string giving the name of the outcome variable.
#' This variable is excluded from the predictor set.
#'
#' @param matching Character vector of variable names used for matching.
#'   These variables are excluded from the main effects.
#'
#' @param include_interactions Logical. If \code{TRUE}, include interaction
#'   terms between the exposure and all other covariates.
#'
#' @param df_bols Integer specifying the degrees of freedom for linear
#'   base learners (\code{bols}).
#'
#' @param df_bbs Integer specifying the degrees of freedom for smooth
#'   base learners (\code{bbs}).
#'
#' @param intercept Logical. Whether linear base learners should include
#'   an intercept term (usually \code{FALSE} for boosting).
#'
#' @param center Logical. Whether continuous covariates should be centered
#'   in smooth base learners.
#'
#' @param flexible Logical. If \code{TRUE}, both linear and smooth base
#'   learners are included for continuous covariates; otherwise, only
#'   linear base learners are used.
#'
#' @param singular Character vector, don't include these as interactions, as they lead to singular cols.
#'
#' @return A list with components:
#'   \describe{
#'     \item{form}{A model formula suitable for \code{mboost::gamboost}.}
#'     \item{all_terms}{A character vector of all base learner terms used
#'       to construct the formula.}
#'   }
#'
#' @import utils
#' @export
generate_formula <- function(data,
                             exposure = NULL,
                             response = "resp",
                             strata = "strata",
                             outcome = "y",
                             matching = "s",
                             include_interactions = TRUE,
                             df_bols = 1,
                             df_bbs = 1,
                             intercept = FALSE,
                             center = TRUE,
                             flexible = TRUE,
                             singular = NULL,
                             boosting_interactions = NULL){
  vars <- setdiff(names(data), c(response, strata, outcome, matching))
  is_cont <- sapply(data[, vars, drop = FALSE], function(x)
    is.numeric(x) && length(unique(x)) > 2)

  cont_vars <- vars[is_cont]
  cat_vars  <- vars[!is_cont]

  make_term <- function(type, var) {
    switch(
      type,
      bols = paste0("bols(", var, ", intercept = ", intercept, ", df = ", df_bols, ")"),
      bbs  = paste0("bbs(", var, ", df = ", df_bbs, ", center = ", center, ")")
    )
  }

  base_terms <- c(vapply(cat_vars, make_term, "", type = "bols"), if (flexible) {
    c(
      vapply(cont_vars, make_term, "", type = "bols"),
      vapply(cont_vars, make_term, "", type = "bbs")
    )
  } else {
    vapply(cont_vars, make_term, "", type = "bols")
  })

  interaction_terms <- character(0)
  if (include_interactions && is.null(boosting_interactions)) {
    if (!is.null(exposure) && exposure %in% vars) {
      # Interactions only with exposure
      inter_vars <- setdiff(c(vars, matching), c(exposure, singular))

      interaction_terms <- vapply(inter_vars, function(v)
        paste0(
          "bols(",
          exposure,
          ", by = ",
          v,
          ", intercept = ",
          intercept,
          ", df = ",
          df_bols,
          ")"
        ), character(1))

    } else {
      # exposure is NULL => include all pairwise interactions among vars + matching
      all_vars <- c(vars, matching)

      # generate all unique combinations
      combs <- t(combn(all_vars, 2))

      interaction_terms <- apply(combs, 1, function(pair) {
        paste0("bols(",
               pair[1],
               ", by = ",
               pair[2],
               ", intercept = ",
               intercept,
               ", df = ",
               df_bols,
               ")")
      })
    }
  }
  else if (!is.null(boosting_interactions)){
    interaction_terms <- vapply(boosting_interactions, function(v)
      paste0(
        "bols(",
        exposure,
        ", by = ",
        v,
        ", intercept = ",
        intercept,
        ", df = ",
        df_bols,
        ")"
      ), character(1))
  }

  terms <- c(base_terms, interaction_terms)

  list(form = as.formula(paste(
    response, "~", paste(terms, collapse = " + ")
  )), all_terms = terms)
}
