#' ADF Test with Automatic Lag Selection
#'
#' @description
#' Performs an Augmented Dickey-Fuller (ADF) test on a time series with
#' automatic lag selection via BIC or fixed lag order.
#'
#' This implements the ADFRC procedure from the original GAUSS code.
#'
#' @param y Numeric vector of the time series (in levels).
#' @param method Lag selection method: 1 = automatic (BIC), 0 = fixed.
#' @param p_max Maximum lag order to consider.
#'
#' @return List with:
#'   \itemize{
#'     \item{t_adf}{ADF t-statistic}
#'     \item{rho_adf}{Estimated rho coefficient}
#'     \item{p_sel}{Selected lag order}
#'   }
#'
#' @keywords internal
.adfrc <- function(y, method = 1, p_max = 4) {
  
  y <- as.numeric(y)
  y <- y[!is.na(y)]
  n <- length(y)
  
  if (n < p_max + 3) {
    return(list(t_adf = NA_real_, rho_adf = NA_real_, p_sel = 0))
  }
  
  # Compute first difference
  dy <- diff(y)
  y_lag <- y[-n]  # y_{t-1}
  
  # Automatic lag selection via BIC
  if (method == 1) {
    bic_vals <- numeric(p_max + 1)
    
    for (p in 0:p_max) {
      adf_fit <- .adf_regression(dy, y_lag, p)
      if (is.null(adf_fit)) {
        bic_vals[p + 1] <- Inf
      } else {
        k <- p + 1  # Number of parameters (rho + p lags)
        n_eff <- length(adf_fit$residuals)
        ssr <- sum(adf_fit$residuals^2)
        bic_vals[p + 1] <- n_eff * log(ssr / n_eff) + k * log(n_eff)
      }
    }
    
    p_sel <- which.min(bic_vals) - 1
  } else {
    p_sel <- p_max
  }
  
  # Final ADF regression with selected lag
  adf_fit <- .adf_regression(dy, y_lag, p_sel)
  
  if (is.null(adf_fit)) {
    return(list(t_adf = NA_real_, rho_adf = NA_real_, p_sel = p_sel))
  }
  
  # Extract t-statistic for rho (coefficient on y_{t-1})
  # In ADF regression: dy_t = rho * y_{t-1} + sum(gamma_j * dy_{t-j}) + e_t
  # H0: rho = 0 (unit root)
  
  rho_adf <- adf_fit$coefficients[1]
  se_rho <- adf_fit$se[1]
  t_adf <- rho_adf / se_rho
  
  list(t_adf = t_adf, rho_adf = rho_adf, p_sel = p_sel)
}


#' ADF Regression Helper
#'
#' @description
#' Runs the ADF regression with specified lag order.
#'
#' @param dy First differences of y.
#' @param y_lag Lagged levels y_{t-1}.
#' @param p Lag order for augmentation.
#'
#' @return List with coefficients, standard errors, and residuals.
#'
#' @keywords internal
.adf_regression <- function(dy, y_lag, p) {
  
  n <- length(dy)
  
  if (n < p + 2) {
    return(NULL)
  }
  
  # Build lagged differences for augmentation
  if (p > 0) {
    # Create lagged differences matrix
    dy_lags <- matrix(NA, nrow = n, ncol = p)
    for (j in 1:p) {
      dy_lags[(j + 1):n, j] <- dy[1:(n - j)]
    }
    
    # Trim to valid observations
    valid_start <- p + 1
    dy_valid <- dy[valid_start:n]
    y_lag_valid <- y_lag[valid_start:n]
    dy_lags_valid <- dy_lags[valid_start:n, , drop = FALSE]
    
    X <- cbind(y_lag_valid, dy_lags_valid)
  } else {
    dy_valid <- dy
    y_lag_valid <- y_lag
    X <- matrix(y_lag_valid, ncol = 1)
  }
  
  y_vec <- dy_valid
  n_eff <- length(y_vec)
  k <- ncol(X)
  
  if (n_eff <= k) {
    return(NULL)
  }
  
  # OLS estimation
  XtX <- crossprod(X)
  Xty <- crossprod(X, y_vec)
  
  # Check for singularity
  XtX_inv <- tryCatch(solve(XtX), error = function(e) NULL)
  if (is.null(XtX_inv)) {
    return(NULL)
  }
  
  beta <- as.vector(XtX_inv %*% Xty)
  residuals <- y_vec - X %*% beta
  
  # Standard errors
  sigma2 <- sum(residuals^2) / (n_eff - k)
  var_beta <- sigma2 * XtX_inv
  se <- sqrt(diag(var_beta))
  
  list(coefficients = beta, se = se, residuals = as.vector(residuals))
}
