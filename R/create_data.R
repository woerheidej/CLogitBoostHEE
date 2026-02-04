#' Simulate a matched case-control dataset
#'
#' Generates a matched case-control dataset with binary exposure,
#' continuous and binary confounders, district strata, and optional
#' heterogeneous exposure effects.
#'
#' @param set Character. Name of the base simulation scenario for outcome generation.
#' @param set.het Character. Name of the scenario used for heterogeneous exposure effects.
#' @param snr Character. Signal-to-noise ratio for main outcome effects. Options: "normal", "strong", "weak".
#' @param snr.x Character. Signal-to-noise ratio for exposure model. Options: "normal", "strong", "weak".
#' @param snr.het Character. Signal-to-noise ratio for heterogeneous exposure effects. Options: "normal", "strong", "weak".
#' @param p Integer. Number of covariates (confounders) to simulate.
#' @param N Integer. Total number of individuals in the population.
#' @param d Integer. Number of districts (clusters) in the population.
#' @param exp.effect Numeric. Baseline log-odds effect of exposure on outcome.
#' @param het.exposure.strength Numeric. Strength of effect modification by exposure heterogeneity.
#' @param n_strata Integer. Number of matched strata (cases) to sample.
#' @param n_het_effects Integer. Number of covariates that contribute to heterogeneous effects.
#' @param k Integer. Number of individuals per stratum (1 case + k-1 controls).
#' @param rr Numeric. Correlation coefficient for compound symmetry covariance structure.
#' @param rho Numeric. Correlation parameter for Toeplitz covariance structure.
#' @param corr Character. Covariate correlation structure. Options: "iid", "cs", "Toeplitz".
#' @param alpha.fac Numeric. Scaling factor for random district and sex effects.
#' @param int Numeric. Intercept for outcome model.
#' @param int.x Numeric. Intercept for exposure model.
#' @param clr_mode Logical. If TRUE, output coefficient names in CLR-style; otherwise use boosting-style names.
#'
#' @return A list with two elements:
#' \describe{
#'   \item{data}{A data.frame containing the matched case-control dataset. Columns include:
#'               exposure (X), strata, covariates (Z1,...,Zp), outcome (y), and sex (s).}
#'   \item{true_coefs}{A named numeric vector of the true coefficients for main and heterogeneous effects.}
#' }
#'
#' @examples
#' \dontrun{
#' sim <- create_data(
#'   set = "SetA",
#'   set.het = "SetHet",
#'   p = 10,
#'   n_strata = 500,
#'   k = 4
#' )
#' }
#'
#' @export
create_data <- function(set = "SetA",
                        set.x = "SetA",
                        set.het = "SetHet",
                        snr = "normal",
                        snr.x = "normal",
                        snr.het = "normal",
                        p = 10,
                        N = 5e5,
                        d = 1e3,
                        exp.effect = 0.6,
                        het.exposure.strength = 0.6,
                        n_strata = 500,
                        n_het_effects = 2,
                        k = 4,
                        rr = 0,
                        rho = 0,
                        corr = c("Toeplitz","iid"),
                        alpha.fac = 1,
                        int = -2,
                        int.x = -1,
                        clr_mode = FALSE) {
  # sample Sex as a binary var
  s <- rbinom(N, 1, 0.5)
  s.effect <- runif(1, -1, 1)/2
  s.effect.x <- runif(1, -1, 1)/2
  alpha.Sex <- c(-s.effect*alpha.fac, s.effect*alpha.fac) # add effect for one Sex
  alpha.Sex.x <- c(-s.effect.x*alpha.fac, s.effect.x*alpha.fac) # add effect for one Sex

  # sample districts
  districts <- sample(1:d, N, replace = TRUE)
  alpha <- runif(d, -2, 2) * alpha.fac
  alpha.x <- runif(d, -2, 2) * alpha.fac

  #  sample confounders
  corr <- match.arg(corr)
  # independent and identically distributed
  if (corr == "iid") {
    Sigma <- diag(p)
  } else if (corr == "cs") {
    Sigma <- matrix(rr, nrow = p, ncol = p); diag(Sigma) <- 1
  } else if (corr == "Toeplitz") {
    if (abs(rho) >= 1) stop("rho must be in (-1, 1) for a valid Toeplitz covariance")
    Sigma <- toeplitz(rho^(0:(p - 1)))
  }
  Z <- mvtnorm::rmvnorm(N, sigma = Sigma)
  colnames(Z) <- paste0("Z", seq(1:ncol(Z)))

  ## create some binary variables
  if(p>2) {
  Z[, 1] <- as.numeric(Z[, 1] > median(Z[,1]))
  Z[, 2] <- as.numeric(Z[, 2] > median(Z[,2]))
  }
  ## define correct SNR factor for x
  snr.fac.x <- 1
  if (snr.x == "strong") {
    snr.fac.x <- 1.5
  }
  if (snr.x == "weak") {
    snr.fac.x <- 0.5
  }

  ## create how confounders affect exposure
  fun.x.obj <- do.call(paste0("fun", set.x), list(Z = Z, exp.effect = exp.effect, snr = snr.fac.x))
  fun.x <- if (is.list(fun.x.obj)) fun.x.obj$effect else fun.x.obj
  true.coefs <- if (is.list(fun.x.obj)) fun.x.obj$coefs else rep(NA, p)

  eta.x <- int.x + alpha.x[districts] + alpha.Sex.x[s + 1] + fun.x
  prob.x <- exp(eta.x) / (1 + exp(eta.x))

  ## draw exposure variable
  X <- rbinom(n = N, size = 1, prob = prob.x)

  ## define correct SNR factor
  snr.fac <- 1
  if (snr == "strong") {
    snr.fac <- 1.5
  }
  if (snr == "weak") {
    snr.fac <- 0.5
  }

  ## create function how confounders affect outcome
  fun.y.obj <- do.call(paste0("fun", set), list(Z = Z, snr = snr.fac, exp.effect = exp.effect))
  fun.y <- if (is.list(fun.y.obj))
    fun.y.obj$effect
  else
    fun.y.obj
  true.coefs <- if (is.list(fun.y.obj))
    fun.y.obj$coefs
  else
    rep(NA, p)
  names(true.coefs) <- colnames(Z)

  ## define correct SNR factor
  snr.fac.het <- 1
  if (snr.het == "strong") {
    snr.fac.het <- 1.5
  }
  if (snr.het == "weak") {
    snr.fac.het <- 0.5
  }

  ## create function how confounders + Sex interact with exposure to affect outcome
  fun.het.obj <- do.call(
    paste0("fun", set.het),
    list(
      X = X,
      Z = Z,
      s = s,
      het.exposure.strength = het.exposure.strength,
      snr = snr.fac.het,
      n_het_effects = n_het_effects,
      signs = fun.y.obj$signs
    )
  )
  fun.het <- if (is.list(fun.het.obj))
    fun.het.obj$effect
  else
    fun.het.obj
  true.het.coefs <- if (is.list(fun.het.obj))
    fun.het.obj$coefs
  else
    rep(NA, p)
  if(clr_mode){
    names(true.het.coefs) <- paste0("X * ", colnames(cbind(Z, s)))
  }else{
    names(true.het.coefs) <- paste0("bols(X, by = ", colnames(cbind(Z, s)), ", intercept = FALSE, df = 1)")
  }

  ## add effect modification for outcome
  sign <- sample(c(-1, 1), 1) #sample if -1 or 1 for effect of X
  eta <- int + alpha[districts] + alpha.Sex[s + 1] + sign * exp.effect * X + fun.het * X  + fun.y
  ## old: eta <- int + Xd%*%alpha + Sex %*% alpha.Sex + exp.effect*X  + fun.y

  probs <- exp(eta) / (1 + exp(eta))


  ## draw outcome variable
  Y <- rbinom(n = N, size = 1, prob = probs)
  ## separate cases and controls in population
  Cases <- (1:N)[Y == 1]
  Controls <- (1:N)[Y == 0]


  ## find cases and controls for study
  cases <- sample(Cases, n_strata)
  matches <- matrix(NA, ncol = k - 1, nrow = n_strata)

  for (i in seq_len(n_strata)) {
    con.i <- Controls[districts[Controls] == districts[cases[i]] &
                        s[Controls]   == s[cases[i]]]
    match.i <- sample(con.i, k - 1)
    matches[i, ] <- match.i
    Controls <- setdiff(Controls, match.i)
  }

  obs.list <- c(t(cbind(cases, matches)))


  ## create study data set
  y <- Y[obs.list]
  X <- X[obs.list]
  z <- Z[obs.list, ]
  d <- districts[obs.list]
  s <- s[obs.list]

  strata <- rep(1:n_strata, each = k)
  sim.data <- as.data.frame(cbind(X, strata, z))
  names(sim.data)[3:(p + 2)] <- paste0("Z", 1:p)
  sim.data$y <- y
  sim.data$s <- s

  values <- c(sign * exp.effect, true.coefs, true.het.coefs)
  name_vec <- c("X", names(true.coefs), names(true.het.coefs))

  true_coefs <- setNames(values, name_vec)

  return(list(data = sim.data, true_coefs = true_coefs))
}

#' Simulate effects for scenario SetC
#'
#' Generates sparse (30%) exposure affecting effects for a set of covariates.
#'
#' @param Z Numeric matrix of covariates (columns = variables, rows = observations).
#' @param exp.effect Numeric. Base effect size for covariates.
#' @param snr Numeric. Signal-to-noise ratio multiplier.
#'
#' @return A list with components:
#' \describe{
#'   \item{effect}{Numeric vector of simulated effects for each observation.}
#'   \item{coefs}{Numeric vector of true coefficients assigned to each covariate.}
#'   \item{signs}{Sign of the effect for each coefficient (+/-).}
#' }
#'
#' @export
funSetC <- function(Z, exp.effect, snr) {
  p <- ncol(Z)
  par <- c(
    sample(c(-exp.effect, exp.effect), 2, replace = FALSE),
    sample(rep(c(exp.effect, -exp.effect) /(1/0.5),p-2), p - 2, replace = FALSE))
  par[sample(p, 0.5 * p)] <- 0
  list(
    effect = as.vector(Z %*% par) * snr,
    coefs  = par * snr
  )
}


#' Simulate main effects for scenario SetA
#'
#' Generates confounding main effects for a set of covariates,
#' including binary and continuous variables.
#'
#' @param Z Numeric matrix of covariates (columns = variables, rows = observations).
#' @param exp.effect Numeric. Base effect size for covariates.
#' @param snr Numeric. Signal-to-noise ratio multiplier.
#'
#' @return A list with components:
#' \describe{
#'   \item{effect}{Numeric vector of simulated effects for each observation.}
#'   \item{coefs}{Numeric vector of true coefficients assigned to each covariate.}
#'   \item{signs}{Sign of the effect for each coefficient (+/-).}
#' }
#'
#' @export
funSetA <- function(Z, exp.effect, snr) {
  par <- rep(0, ncol(Z))
  par[1:2] <- sample(c(-exp.effect, exp.effect), replace = FALSE)
  par[3:ncol(Z)] <- sample(rep(c(exp.effect / (1/0.5), -exp.effect / (1/0.5)),
                               times = (ncol(Z) - 2) / 2), replace = FALSE)
  signs <- ifelse(par < 0, -1, +1)
  signs <- c(signs, 1)
  effect <- (Z %*% par)*snr
  list(
    effect = effect,
    coefs = par * snr,
    signs = signs
  )
}

#' Simulate heterogeneous exposure effects for scenario SetHet
#'
#' Generates heterogeneous exposure effects for covariates and a binary variable (s).
#'
#' @param Z Numeric matrix of covariates.
#' @param s Numeric vector of binary variable (e.g., Sex) of length nrow(Z).
#' @param het.exposure.strength Numeric. Strength of effect modification.
#' @param snr Numeric. Signal-to-noise ratio multiplier.
#' @param n_het_effects Integer. Number of covariates with heterogeneous effects.
#' @param signs Numeric vector of effect signs (+1/-1) corresponding to columns of Z.
#'
#' @return A list with components:
#' \describe{
#'   \item{effect}{Numeric vector of simulated heterogeneous effects.}
#'   \item{coefs}{Numeric vector of true heterogeneous coefficients.}
#'   \item{n_het_effects}{Number of heterogeneous effects simulated.}
#' }
#'
#' @export
funSetHet <- function(X, Z, s, het.exposure.strength, snr, n_het_effects, signs) {
  target_sd <- (sd(X*Z[,1]) + sd(X*Z[,2])) / 2 # mean between the two,
  # exposure ~0.3, Z[,1:2] ~ 0.5 = 0.15 -> SD 0.85/0.15 = ~ 0.35

  Z <- cbind(Z, s)
  par <- rep(0, ncol(Z))

  is_binary <- apply(Z, 2, function(x) length(unique(x)) < 3)
  effect_cols <- sample(1:ncol(Z), n_het_effects)

  for (i in seq_along(effect_cols)) {
    col_idx <- effect_cols[i]
    sign <- signs[col_idx]
    if (is_binary[col_idx]) {
      par[col_idx] <- sign * het.exposure.strength
    } else if (is.na(sign)) {
      sign_sub <- sample(c(-1, 1), 1)
      par[col_idx] <- sign_sub * het.exposure.strength / (sd(X*Z[,col_idx])/target_sd)
    } else {
      par[col_idx] <- sign * het.exposure.strength / (sd(X*Z[,col_idx])/target_sd)
    }
  }

  effect <- Z %*% par
  ret <- effect * snr
  list(
    effect = ret,
    coefs = par * snr,
    n_het_effects = n_het_effects
  )
}

#' Simulate main effects for scenario SetB with nonlinear components
#'
#' Generates a mixture of binary and nonlinear covariate effects.
#'
#' @param Z Numeric matrix of covariates.
#' @param snr Numeric. Signal-to-noise ratio multiplier.
#' @param exp.effect Numeric. Base effect size for covariates.
#'
#' @return A list with components:
#' \describe{
#'   \item{effect}{Numeric vector of simulated effects for each observation.}
#'   \item{coefs}{Numeric vector of true coefficients for linear terms; nonlinear effects are NA.}
#'   \item{signs}{Sign of linear coefficients (+1/-1); nonlinear terms are NA.}
#' }
#'
#' @export
funSetB <- function(Z, snr = 1, exp.effect = 0.6, scale_nonlinear = TRUE) {

  binary_coefs <- sample(c(-exp.effect, exp.effect), size = 2, replace = FALSE) * snr
  binary_effect <- Z[, 1] * binary_coefs[1] + Z[, 2] * binary_coefs[2]

  # Target mean effect (0.5 = mean(binary_cols))
  target_sd <- 0.5 * exp.effect  # 0.3 for 50/50 var

  n_cont <- ncol(Z) - 2

  base_funcs <- list(
    f1 = function(x) (0.6 * sin(2 * x) + 0.6 * cos(x)),
    f2 = function(x) (0.5 * x^2),
    f3 = function(x) (-0.5 * x^2),
    f4 = function(x) (0.25 * x + sin(1.5 * x) * 0.7 * cos(x)),
    f5 = function(x) ((1 / (1 + exp(-2 * x)))-0.5),
    f6 = function(x) ((1 / (1 + exp(2 * x)))-0.5),
    f7 = function(x) (cos(2 * x) * 0.8),
    f8 = function(x) ifelse(x < 0, 0.8 * x, 0.4 * x)
  )

  funcs_n <- length(base_funcs)
  selected_cols <- if (n_cont >= funcs_n) 3:ncol(Z) else sample(3:ncol(Z), funcs_n, replace = TRUE)
  funcs_final <- base_funcs
  for (i in seq_along(base_funcs)) {
    funcs_final[[i]] <- local({
      f_i <- base_funcs[[i]]
      x_i <- Z[, selected_cols[i]]
      current_sd <- sd(f_i(x_i))
      scale_factor <- target_sd / current_sd
      function(x) scale_factor * f_i(x)
    })
  }

  nonlinear_effect <- 0
  for (i in seq_along(funcs_final)) {
    nonlinear_effect <- nonlinear_effect + funcs_final[[i]](Z[, selected_cols[i]])
  }

  total_effect <- (binary_effect + nonlinear_effect)

  coefs <- c(binary_coefs, rep(NA, n_cont))
  signs <- c(ifelse(coefs < 0, -1, 1),1)
  list(
    effect = total_effect,
    coefs = coefs,
    signs = signs
  )
}




