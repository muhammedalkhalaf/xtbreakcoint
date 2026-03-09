#' Iterative Factor-Break Estimation
#'
#' @description
#' Jointly estimates common factors and structural break dates using an
#' iterative procedure based on Banerjee & Carrion-i-Silvestre (2015).
#'
#' @param data Data frame with panel data.
#' @param depvar Name of dependent variable.
#' @param indepvars Names of independent variables.
#' @param id Panel identifier column name.
#' @param time Time identifier column name.
#' @param model Model specification (1-5).
#' @param max_factors Maximum number of factors.
#' @param trim Trimming parameter.
#' @param max_iter Maximum iterations.
#' @param tolerance Convergence tolerance.
#'
#' @return List with estimated factors, residuals, breaks, and diagnostics.
#'
#' @keywords internal
.factcoint_iter <- function(data, depvar, indepvars, id, time,
                            model, max_factors, trim,
                            max_iter, tolerance) {
  
  panels <- unique(data[[id]])
  N <- length(panels)
  time_vals <- sort(unique(data[[time]]))
  Tobs <- length(time_vals)
  Tm1 <- Tobs - 1
  
  # ---- Step 1: Initial cointegrating regression (no breaks, no factors) ----
  # Run pooled OLS on first differences to avoid spurious regression
  
  # Build first-differenced data
  Dy <- matrix(NA, nrow = Tm1, ncol = N)
  Dx <- array(NA, dim = c(Tm1, N, length(indepvars)))
  
  for (i in seq_along(panels)) {
    idx <- data[[id]] == panels[i]
    y_i <- data[[depvar]][idx]
    Dy[, i] <- diff(y_i)
    
    for (k in seq_along(indepvars)) {
      x_ik <- data[[indepvars[k]]][idx]
      Dx[, i, k] <- diff(x_ik)
    }
  }
  
  # Pooled first-differenced regression
  Dy_vec <- as.vector(Dy)
  Dx_mat <- matrix(NA, nrow = Tm1 * N, ncol = length(indepvars))
  for (k in seq_along(indepvars)) {
    Dx_mat[, k] <- as.vector(Dx[, , k])
  }
  
  # Add deterministic components based on model
  det_mat <- .build_deterministics(Tm1, N, model, break_pos = NULL)
  X_full <- cbind(Dx_mat, det_mat)
  
  # OLS estimation
  fit <- stats::lm.fit(x = X_full, y = Dy_vec)
  beta <- fit$coefficients[seq_along(indepvars)]
  
  # Compute first-differenced residuals De (T-1 x N)
  De <- matrix(fit$residuals, nrow = Tm1, ncol = N)
  
  # ---- Initialize breaks ----
  breaks <- rep(0, N)
  
  # ---- Iterative estimation ----
  old_ssr <- sum(De^2)
  converged <- FALSE
  iter <- 0
  n_factors <- 0
  Fhat <- NULL
  Lambda <- NULL
  
  while (!converged && iter < max_iter) {
    iter <- iter + 1
    
    # ---- Step 2: Estimate common factors from residuals ----
    if (max_factors > 0) {
      factor_result <- .estimate_factors(De, max_factors)
      n_factors <- factor_result$n_factors
      Fhat <- factor_result$Fhat
      Lambda <- factor_result$Lambda
      
      # Defactor residuals: De_idio = De - Fhat * Lambda'
      if (n_factors > 0) {
        De_common <- Fhat %*% t(Lambda)
        De_idio <- De - De_common
      } else {
        De_idio <- De
      }
    } else {
      De_idio <- De
      n_factors <- 0
    }
    
    # ---- Step 3: Estimate break dates for each unit ----
    if (model >= 3) {
      for (i in seq_along(panels)) {
        breaks[i] <- .estimate_break(De_idio[, i], Tobs, model, trim)
      }
      
      # For model 5, use common break (average, rounded)
      if (model == 5) {
        common_break <- round(mean(breaks[breaks > 0]))
        if (is.na(common_break)) common_break <- 0
        breaks <- rep(common_break, N)
      }
    }
    
    # ---- Step 4: Re-estimate cointegrating regression with breaks ----
    # Re-run the regression with structural break dummies
    det_mat_new <- .build_deterministics(Tm1, N, model, breaks)
    X_full_new <- cbind(Dx_mat, det_mat_new)
    
    fit_new <- stats::lm.fit(x = X_full_new, y = Dy_vec)
    De_new <- matrix(fit_new$residuals, nrow = Tm1, ncol = N)
    
    # Check convergence
    new_ssr <- sum(De_new^2)
    rel_change <- abs(new_ssr - old_ssr) / (abs(old_ssr) + 1e-10)
    
    if (rel_change < tolerance) {
      converged <- TRUE
    }
    
    De <- De_new
    old_ssr <- new_ssr
  }
  
  # Final factor estimation
  if (max_factors > 0 && n_factors > 0) {
    factor_result <- .estimate_factors(De, max_factors)
    n_factors <- factor_result$n_factors
    Fhat <- factor_result$Fhat
    Lambda <- factor_result$Lambda
    
    De_common <- Fhat %*% t(Lambda)
    De_idio <- De - De_common
  } else {
    De_idio <- De
  }
  
  list(
    Dres = De_idio,
    Fhat = Fhat,
    Lambda = Lambda,
    breaks = breaks,
    n_factors = n_factors,
    iterations = iter,
    ssr = old_ssr
  )
}


#' Estimate Common Factors via Principal Components
#'
#' @description
#' Estimates common factors from a panel of residuals using principal
#' components analysis (PCA). The number of factors is selected using the
#' Bai & Ng (2002) IC criterion.
#'
#' @param De Matrix of first-differenced residuals (T-1 x N).
#' @param max_factors Maximum number of factors to consider.
#'
#' @return List with Fhat (factors), Lambda (loadings), and n_factors.
#'
#' @references
#' Bai, J., & Ng, S. (2002). Determining the number of factors in approximate
#' factor models. \emph{Econometrica}, 70(1), 191-221.
#' \doi{10.1111/1468-0262.00273}
#'
#' @keywords internal
.estimate_factors <- function(De, max_factors) {
  
  Tm1 <- nrow(De)
  N <- ncol(De)
  
  if (max_factors == 0) {
    return(list(Fhat = NULL, Lambda = NULL, n_factors = 0))
  }
  
  # Cap at min(T-1, N) - 1
  k_max <- min(max_factors, min(Tm1, N) - 1)
  if (k_max <= 0) {
    return(list(Fhat = NULL, Lambda = NULL, n_factors = 0))
  }
  
  # PCA on the covariance matrix
  # GAUSS uses eigenvalue decomposition of De'De / (T*N)
  cov_mat <- crossprod(De) / (Tm1 * N)
  
  eig <- eigen(cov_mat, symmetric = TRUE)
  eigenvalues <- eig$values
  eigenvectors <- eig$vectors
  
  # Select number of factors using IC_p2 criterion (Bai & Ng, 2002)
  # IC(k) = log(V(k)) + k * g(N, T)
  # g(N,T) = (N + T)/(N*T) * log(min(N, T))
  
  g_NT <- (N + Tm1) / (N * Tm1) * log(min(N, Tm1))
  
  ic_vals <- numeric(k_max + 1)
  
  # k = 0: no factors
  V0 <- sum(De^2) / (Tm1 * N)
  ic_vals[1] <- log(V0)
  
  for (k in seq_len(k_max)) {
    # Factors: sqrt(T) * first k eigenvectors
    Lambda_k <- eigenvectors[, 1:k, drop = FALSE] * sqrt(N)
    Fhat_k <- De %*% eigenvectors[, 1:k, drop = FALSE] / sqrt(N)
    
    # Residual variance
    De_fitted <- Fhat_k %*% t(Lambda_k)
    resid_k <- De - De_fitted
    Vk <- sum(resid_k^2) / (Tm1 * N)
    
    ic_vals[k + 1] <- log(Vk) + k * g_NT
  }
  
  # Select k that minimizes IC
  k_opt <- which.min(ic_vals) - 1
  
  if (k_opt == 0) {
    return(list(Fhat = NULL, Lambda = NULL, n_factors = 0))
  }
  
  # Extract optimal factors and loadings
  Lambda <- eigenvectors[, 1:k_opt, drop = FALSE] * sqrt(N)
  Fhat <- De %*% eigenvectors[, 1:k_opt, drop = FALSE] / sqrt(N)
  
  list(Fhat = Fhat, Lambda = Lambda, n_factors = k_opt)
}


#' Estimate Break Date for a Single Unit
#'
#' @description
#' Estimates the structural break date by minimizing the sum of squared
#' residuals over all possible break points.
#'
#' @param resid Vector of residuals.
#' @param Tobs Total number of time periods.
#' @param model Model specification.
#' @param trim Trimming parameter.
#'
#' @return Estimated break position (0 if no break).
#'
#' @keywords internal
.estimate_break <- function(resid, Tobs, model, trim) {
  
  Tm1 <- length(resid)
  
  # Trimming: exclude first and last trim*T observations
  trim_obs <- max(1, floor(trim * Tm1))
  start_pos <- trim_obs + 1
  end_pos <- Tm1 - trim_obs
  
  if (start_pos >= end_pos) {
    return(0)
  }
  
  # Grid search over break positions
  ssr_min <- Inf
  break_opt <- 0
  
  for (tb in start_pos:end_pos) {
    # Build break dummy
    D_break <- c(rep(0, tb), rep(1, Tm1 - tb))
    
    # Regression of residuals on break dummy
    X <- cbind(1, D_break)
    fit <- stats::lm.fit(x = X, y = resid)
    ssr <- sum(fit$residuals^2)
    
    if (ssr < ssr_min) {
      ssr_min <- ssr
      break_opt <- tb
    }
  }
  
  break_opt
}


#' Build Deterministic Component Matrix
#'
#' @description
#' Constructs the matrix of deterministic components (constant, trend,
#' break dummies) for the panel regression.
#'
#' @param Tm1 Number of first-differenced observations.
#' @param N Number of cross-sections.
#' @param model Model specification (1-5).
#' @param break_pos Vector of break positions (or NULL).
#'
#' @return Matrix of deterministic components.
#'
#' @keywords internal
.build_deterministics <- function(Tm1, N, model, break_pos) {
  
  # Model 1: constant only (no deterministics in first diffs)
  # Model 2: constant + trend (trend in first diffs = constant)
  # Model 3: constant + level shift (level shift in diffs = impulse dummy)
  # Model 4: constant + trend + level shift
  # Model 5: constant + trend + level + slope shift
  
  n_obs <- Tm1 * N
  
  det_list <- list()
  
  # Trend component (appears as constant in first differences)
  if (model >= 2) {
    det_list$trend <- rep(1, n_obs)
  }
  
  # Break components
  if (model >= 3 && !is.null(break_pos)) {
    # Level shift in first differences = impulse at break
    D_level <- numeric(n_obs)
    
    for (i in seq_len(N)) {
      tb <- break_pos[i]
      if (tb > 0 && tb <= Tm1) {
        row_idx <- (i - 1) * Tm1 + tb
        if (row_idx <= n_obs) {
          D_level[row_idx] <- 1
        }
      }
    }
    det_list$level_shift <- D_level
    
    # Trend shift (model 4 and 5)
    if (model >= 4) {
      D_trend <- numeric(n_obs)
      for (i in seq_len(N)) {
        tb <- break_pos[i]
        if (tb > 0 && tb < Tm1) {
          for (t in (tb + 1):Tm1) {
            row_idx <- (i - 1) * Tm1 + t
            if (row_idx <= n_obs) {
              D_trend[row_idx] <- 1
            }
          }
        }
      }
      det_list$trend_shift <- D_trend
    }
    
    # Slope shift (model 5) - affects the relationship, not just deterministics
    if (model == 5) {
      D_slope <- numeric(n_obs)
      for (i in seq_len(N)) {
        tb <- break_pos[i]
        if (tb > 0 && tb < Tm1) {
          for (t in (tb + 1):Tm1) {
            row_idx <- (i - 1) * Tm1 + t
            if (row_idx <= n_obs) {
              D_slope[row_idx] <- t - tb
            }
          }
        }
      }
      det_list$slope_shift <- D_slope
    }
  }
  
  if (length(det_list) == 0) {
    return(matrix(1, nrow = n_obs, ncol = 1))  # At least a constant
  }
  
  do.call(cbind, det_list)
}
