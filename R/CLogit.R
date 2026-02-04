#' Implements the CLogit Family for matched
#' case-control data, to be used with \code{mboost::gamboost}.
#'
#' @param risktype Character string. Either \code{"classic"} (mean risk)
#'   or \code{"median"} (median risk).
#'
#' @return An object of class \code{Family} usable in \code{gamboost}.
#'
#' @export
CLogit <- function(risktype = c("classic", "median"))
{
  risktype <- match.arg(risktype, choices = c("classic", "median"))
  ngradient <- function(y, f, w = 1) {

    event <- y[, 1]
    strata <- y[, 2]
    n <- length(event)
    if (length(f) == 1)
      f <- rep(f, n)
    if (length(w) == 1)
      w <- rep(w, n)

    strata <- factor(strata, levels = unique(strata))

    fj <- split(f, strata)
    wj <- split(w, strata)
    ej <- split(event, strata)

    ngrad <- c(unlist(mapply(
      FUN = gradfun_help, ej, fj, wj, USE.NAMES = FALSE
    )))
    ngrad
  }
  risk <- function(y, f, w = 1) {
    event <- y[, 1]
    strata <- y[, 2]
    n <- length(event)
    if (length(f) == 1)
      f <- rep(f, n)
    if (length(w) == 1)
      w <- rep(w, n)
    indx <- rep(1:n, w)

    event <- event[indx]
    strata <- strata[indx]
    f <- f[indx]

    risk <- rep(0, length(event))
    for (i in 1:length(event)) {
      risk[i] <- exp(f[i]) / sum(exp(f[strata == strata[i]]))
    }
    risk <- log(risk[event == 1])

    if (risktype == "median") {
      -median(risk)
    } else{
      -mean(risk)
    }
  }

  Family(
    ngradient =  ngradient,
    risk = risk,
    gradfun_help <- function(y, f, w) {
      sum(y) * w * (y - exp(f) / sum(exp(f)))
    },
    check_y = function(y) {
      stopifnot("Stratum information within response variable has to be sorted!" = sort(y[, 2]) ==
                  y[, 2])
      return(y)
    },
    name = "Negative Conditional Log-Likelihood",
    offset = function(y, w) {
      0
    }
  )
}
