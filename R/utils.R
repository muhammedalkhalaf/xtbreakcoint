#' Parse Model Specification
#'
#' @description
#' Parses model specification from string or integer to internal format.
#'
#' @param model Model specification (string or integer).
#'
#' @return List with model number and label.
#'
#' @keywords internal
.parse_model <- function(model) {
  
  if (is.numeric(model)) {
    model_num <- as.integer(model)
    if (model_num < 1 || model_num > 5) {
      stop("Model number must be between 1 and 5")
    }
  } else {
    model_str <- tolower(as.character(model))
    model_num <- switch(model_str,
      "constant" = 1L,
      "1" = 1L,
      "trend" = 2L,
      "2" = 2L,
      "levelshift" = 3L,
      "3" = 3L,
      "trendshift" = 4L,
      "4" = 4L,
      "regimeshift" = 5L,
      "5" = 5L,
      stop("Invalid model specification. Choose: constant, trend, ",
           "levelshift, trendshift, or regimeshift")
    )
  }
  
  model_label <- switch(model_num,
    "Constant (Model 1)",
    "Constant + Trend (Model 2)",
    "Constant + Level Shift (Model 3)",
    "Constant + Trend + Level Shift (Model 4)",
    "Constant + Trend + Level + Slope Shift (Model 5)"
  )
  
  list(num = model_num, label = model_label)
}


#' Get Empirical Moments
#'
#' @description
#' Returns empirical moments (mean and variance) of the ADF t-statistic
#' under the null hypothesis of no cointegration. These are from Monte Carlo
#' simulations in the original GAUSS code (T=100, 100,000 replications).
#'
#' @param model Model specification (1-5).
#'
#' @return List with mean, variance, and lambda_dependent flag.
#'
#' @keywords internal
.get_moments <- function(model) {
  
 # From GAUSS brkfactors_heterog.gss (lines 79-106)
  # Monte Carlo moments for T=100, 100,000 replications
  
  if (model %in% c(1, 3)) {
    # Constant-type models (no trend in DGP)
    mean_t <- -0.41632799
    var_t <- 0.98339487
    lambda_dependent <- FALSE
  } else if (model %in% c(2, 4)) {
    # Trend-type models
    mean_t <- -1.5377067
    var_t <- 0.35005403
    lambda_dependent <- FALSE
  } else if (model == 5) {
    # Regime shift: moments depend on break fraction (lambda)
    # Default values (will be updated based on estimated lambda)
    mean_t <- -1.998013
    var_t <- 0.35833849
    lambda_dependent <- TRUE
  } else {
    stop("Invalid model number")
  }
  
  list(mean = mean_t, var = var_t, lambda_dependent = lambda_dependent)
}


#' Get Lambda-Dependent Moments for Model 5
#'
#' @description
#' Returns empirical moments for Model 5 (regime shift) that depend on the
#' break fraction lambda. Values from GAUSS code for lambda = 0.1 to 0.9.
#'
#' @param lambda Break fraction (between 0 and 1).
#'
#' @return List with mean and variance.
#'
#' @keywords internal
.get_lambda_moments <- function(lambda) {
  
  # GAUSS mean_t and var_t vectors (9 values, lambda 0.1 to 0.9)
  # Interpolate or select closest bin
  
  lambda_grid <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
  
  mean_grid <- c(
    -1.6803178,  # 0.1
    -1.8163351,  # 0.2
    -1.9198423,  # 0.3
    -1.9805257,  # 0.4
    -1.998013,   # 0.5
    -1.9752734,  # 0.6
    -1.9125286,  # 0.7
    -1.816865,   # 0.8
    -1.6755147   # 0.9
  )
  
  var_grid <- c(
    0.40488013,  # 0.1
    0.41454518,  # 0.2
    0.40165997,  # 0.3
    0.36829752,  # 0.4
    0.35833849,  # 0.5
    0.36808259,  # 0.6
    0.39040626,  # 0.7
    0.4229098,   # 0.8
    0.39749512   # 0.9
  )
  
  # Clamp lambda to [0.1, 0.9]
  lambda <- max(0.1, min(0.9, lambda))
  
  # Find closest bin
  idx <- which.min(abs(lambda_grid - lambda))
  
  list(mean = mean_grid[idx], var = var_grid[idx])
}
