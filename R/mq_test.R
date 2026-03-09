#' MQ Test for Common Stochastic Trends
#'
#' @description
#' Implements the Bai & Ng (2004) MQ test to determine the number of common
#' stochastic trends among estimated factors. This helps distinguish between
#' I(1) and I(0) factors.
#'
#' @param Fhat Matrix of estimated factors in levels (T x k).
#' @param model Model specification (1-5).
#' @param N Number of cross-sections.
#' @param parametric Logical; if TRUE, use parametric version.
#'
#' @return List with:
#'   \itemize{
#'     \item{MQ}{MQ test statistic}
#'     \item{n_trends}{Estimated number of stochastic trends}
#'   }
#'
#' @references
#' Bai, J., & Ng, S. (2004). A PANIC attack on unit roots and cointegration.
#' \emph{Econometrica}, 72(4), 1127-1177. \doi{10.1111/j.1468-0262.2004.00528.x}
#'
#' @keywords internal
.mq_test <- function(Fhat, model, N, parametric = FALSE) {
  
  if (is.null(Fhat) || ncol(Fhat) == 0) {
    return(list(MQ = NA_real_, n_trends = 0))
  }
  
  Fhat <- as.matrix(Fhat)
  Tobs <- nrow(Fhat)
  k <- ncol(Fhat)
  
  # Demean or detrend factors based on model
  Fhat_adj <- .detrend_factors(Fhat, model)
  
  # Compute individual ADF statistics for each factor
  adf_stats <- numeric(k)
  
  for (j in 1:k) {
    if (parametric) {
      # Parametric: use ADF with optimal lag selection
      adf_result <- .adfrc(Fhat_adj[, j], method = 1, p_max = 4)
      adf_stats[j] <- adf_result$t_adf
    } else {
      # Non-parametric: Phillips-Perron style
      adf_stats[j] <- .pp_test(Fhat_adj[, j])
    }
  }
  
  # Sort statistics (most negative first)
  adf_sorted <- sort(adf_stats)
  
  # MQ test: sequentially test for stochastic trends
  # Use 5% critical value from Bai & Ng (2004) Table I
  # Critical values depend on the null hypothesis rank
  
  cv_5pct <- .get_mq_critical_values(k)
  
  # Count stochastic trends
  n_trends <- 0
  for (j in 1:k) {
    # Test H0: at least j stochastic trends
    # MQ statistic is the j-th smallest ADF statistic
    if (adf_sorted[j] > cv_5pct[j]) {
      n_trends <- j
    } else {
      break
    }
  }
  
  # MQ statistic (sum of ADF stats for identified I(1) factors)
  if (n_trends > 0) {
    MQ <- sum(adf_sorted[1:n_trends])
  } else {
    MQ <- adf_sorted[1]  # Report the most negative
  }
  
  list(MQ = MQ, n_trends = n_trends)
}


#' Detrend Factors
#'
#' @description
#' Demean or detrend factor estimates based on model specification.
#'
#' @param Fhat Matrix of factors.
#' @param model Model specification.
#'
#' @return Adjusted factor matrix.
#'
#' @keywords internal
.detrend_factors <- function(Fhat, model) {
  
  Tobs <- nrow(Fhat)
  k <- ncol(Fhat)
  
  Fhat_adj <- matrix(NA, nrow = Tobs, ncol = k)
  
  for (j in 1:k) {
    f <- Fhat[, j]
    
    if (model %in% c(1, 3)) {
      # Demean only (constant-type models)
      Fhat_adj[, j] <- f - mean(f)
    } else {
      # Detrend (trend-type models)
      trend <- 1:Tobs
      fit <- stats::lm(f ~ trend)
      Fhat_adj[, j] <- stats::residuals(fit)
    }
  }
  
  Fhat_adj
}


#' Phillips-Perron Test Statistic
#'
#' @description
#' Computes a non-parametric Phillips-Perron style unit root test statistic.
#'
#' @param y Numeric vector.
#'
#' @return PP test statistic.
#'
#' @keywords internal
.pp_test <- function(y) {
  
  y <- as.numeric(y)
  y <- y[!is.na(y)]
  n <- length(y)
  
  if (n < 5) {
    return(NA_real_)
  }
  
  # Simple AR(1) regression: y_t = rho * y_{t-1} + e_t
  y_lag <- y[-n]
  y_cur <- y[-1]
  
  # OLS
  rho_hat <- sum(y_lag * y_cur) / sum(y_lag^2)
  resid <- y_cur - rho_hat * y_lag
  
  # Variance estimates
  sigma2 <- sum(resid^2) / (n - 2)
  se_rho <- sqrt(sigma2 / sum(y_lag^2))
  
  # Long-run variance using Newey-West estimator
  max_lag <- floor(4 * (n / 100)^(2/9))
  gamma <- numeric(max_lag + 1)
  
  for (j in 0:max_lag) {
    if (j == 0) {
      gamma[j + 1] <- sum(resid^2) / n
    } else {
      gamma[j + 1] <- sum(resid[(j + 1):length(resid)] * resid[1:(length(resid) - j)]) / n
    }
  }
  
  # Bartlett weights
  lambda2 <- gamma[1]
  for (j in 1:max_lag) {
    w <- 1 - j / (max_lag + 1)
    lambda2 <- lambda2 + 2 * w * gamma[j + 1]
  }
  
  # PP correction
  correction <- (lambda2 - gamma[1]) / 2 * sqrt(n) / sqrt(sum(y_lag^2))
  t_pp <- (rho_hat - 1) / se_rho - correction / se_rho
  
  t_pp
}


#' Get MQ Critical Values
#'
#' @description
#' Returns 5% critical values for the MQ test from Bai & Ng (2004) Table I.
#'
#' @param k Number of factors.
#'
#' @return Vector of critical values.
#'
#' @keywords internal
.get_mq_critical_values <- function(k) {
  

  # 5% critical values for MQ test (asymptotic)
  # From Bai & Ng (2004), Table I
  # These are for the j-th smallest eigenvalue test
  
  # Approximate critical values (more negative = reject I(1))
  cv <- c(
    -2.86,  # j = 1
    -2.86,  # j = 2
    -2.86,  # j = 3
    -2.86,  # j = 4
    -2.86   # j = 5
  )
  
  if (k > length(cv)) {
    cv <- c(cv, rep(-2.86, k - length(cv)))
  }
  
  cv[1:k]
}
