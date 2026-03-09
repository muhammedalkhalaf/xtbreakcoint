#' Panel Cointegration Test with Structural Breaks
#'
#' @description
#' Implements the panel cointegration test of Banerjee and Carrion-i-Silvestre
#' (2015), allowing for structural breaks and cross-section dependence through
#' common factors.
#'
#' @details
#' This function tests for panel cointegration in the presence of structural
#' breaks and cross-sectional dependence. The methodology follows these steps:
#'
#' 1. **Iterative Factor-Break Estimation**: Jointly estimates common factors
#'    and individual break dates using an iterative procedure.
#' 2. **Individual ADF Tests**: Applies ADF tests to the defactored
#'    (idiosyncratic) residuals for each cross-section unit.
#' 3. **Panel Test Statistic**: Combines individual ADF statistics into a
#'    standardized panel statistic using Monte Carlo moments.
#' 4. **MQ Test**: If factors are detected, tests whether they are stationary
#'    or represent stochastic trends (Bai & Ng, 2004).
#'
#' @section Model Specifications:
#' \describe{
#'   \item{1 - constant}{Constant only}
#'   \item{2 - trend}{Constant plus linear trend}
#'   \item{3 - levelshift}{Constant plus level shift at break}
#'   \item{4 - trendshift}{Constant, trend, plus level shift (default)}
#'   \item{5 - regimeshift}{Constant, trend, level shift, plus slope shift}
#' }
#'
#' @param formula A formula of the form `y ~ x1 + x2 + ...` specifying the
#'   cointegrating relationship.
#' @param data A data frame containing panel data with columns for the panel
#'   identifier, time identifier, and all variables in the formula.
#' @param id Character string naming the panel (cross-section) identifier.
#' @param time Character string naming the time identifier.
#' @param model Model specification for deterministic components. One of:
#'   `"constant"` (1), `"trend"` (2), `"levelshift"` (3), `"trendshift"` (4,
#'   default), or `"regimeshift"` (5). Can also specify as integer 1-5.
#' @param max_factors Maximum number of common factors to estimate (default: 5).
#'   Set to 0 to skip factor estimation.
#' @param max_lag Maximum lag order for ADF tests (default: 4).
#' @param lag_method Method for selecting ADF lag order: `"auto"` for automatic
#'   selection via BIC (default) or `"fixed"` to use `max_lag`.
#' @param trim Trimming parameter for break estimation, proportion of sample
#'   excluded from endpoints (default: 0.15). Must be in (0, 0.5).
#' @param max_iter Maximum iterations for factor-break estimation (default: 20).
#' @param tolerance Convergence tolerance for iterative estimation (default:
#'   0.001).
#'
#' @return An object of class `"xtbreakcoint"` containing:
#' \describe{
#'   \item{Z_t}{Panel test statistic (standard normal under H0)}
#'   \item{p_value}{One-sided p-value for Z_t}
#'   \item{tbar}{Average of individual ADF t-statistics}
#'   \item{mean_t}{Expected value of t-statistic under H0}
#'   \item{var_t}{Variance of t-statistic under H0}
#'   \item{N}{Number of cross-section units}
#'   \item{T}{Number of time periods}
#'   \item{n_factors}{Number of estimated common factors}
#'   \item{n_trends}{Number of stochastic trends (from MQ test)}
#'   \item{MQ_np}{Non-parametric MQ test statistic}
#'   \item{MQ_p}{Parametric MQ test statistic}
#'   \item{iterations}{Number of iterations for factor-break convergence}
#'   \item{reject_pct}{Percentage of individual units rejecting H0 at 5\%}
#'   \item{adf_stats}{Vector of individual ADF t-statistics}
#'   \item{lag_orders}{Vector of selected lag orders for each unit}
#'   \item{breaks}{Vector of estimated break dates (periods from start)}
#'   \item{factors}{Matrix of estimated common factors (T-1 x n_factors)}
#'   \item{model}{Model specification used}
#'   \item{call}{The matched call}
#' }
#'
#' @references
#' Banerjee, A., & Carrion-i-Silvestre, J. L. (2015). Cointegration in panel
#' data with structural breaks and cross-section dependence. \emph{Journal of
#' Applied Econometrics}, 30(1), 1-22. \doi{10.1002/jae.2348}
#'
#' Bai, J., & Ng, S. (2004). A PANIC attack on unit roots and cointegration.
#' \emph{Econometrica}, 72(4), 1127-1177. \doi{10.1111/j.1468-0262.2004.00528.x}
#'
#' @examples
#' # Generate example panel data
#' set.seed(42)
#' N <- 10  # panels
#' T <- 50  # time periods
#' 
#' # Create cointegrated data with a structural break
#' panel_data <- data.frame(
#'   id = rep(1:N, each = T),
#'   time = rep(1:T, N),
#'   y = NA,
#'   x = NA
#' )
#' 
#' for (i in 1:N) {
#'   idx <- panel_data$id == i
#'   x <- cumsum(rnorm(T))
#'   u <- rnorm(T, sd = 0.5)
#'   # Cointegrating relationship with break at t=25
#'   beta <- ifelse(1:T <= 25, 1, 1.5)
#'   y <- 1 + beta * x + u
#'   panel_data$x[idx] <- x
#'   panel_data$y[idx] <- y
#' }
#' 
#' # Test for cointegration
#' result <- xtbreakcoint(y ~ x, data = panel_data, id = "id", time = "time")
#' print(result)
#'
#' @export
xtbreakcoint <- function(formula, data, id, time,
                         model = "trendshift",
                         max_factors = 5,
                         max_lag = 4,
                         lag_method = c("auto", "fixed"),
                         trim = 0.15,
                         max_iter = 20,
                         tolerance = 0.001) {
  
  # Capture call
  cl <- match.call()
  
  # Match arguments
  lag_method <- match.arg(lag_method)
  

  # Validate model specification
  model_info <- .parse_model(model)
  model_num <- model_info$num
  model_label <- model_info$label
  
  # Validate inputs
  if (!is.data.frame(data)) {
    stop("'data' must be a data frame")
  }
  
  if (!id %in% names(data)) {
    stop("Panel identifier '", id, "' not found in data")
  }
  
  if (!time %in% names(data)) {
    stop("Time identifier '", time, "' not found in data")
  }
  
  if (trim <= 0 || trim >= 0.5) {
    stop("'trim' must be between 0 and 0.5")
  }
  
  if (max_lag < 0) {
    stop("'max_lag' must be non-negative")
  }
  
  if (max_factors < 0) {
    stop("'max_factors' must be non-negative")
  }
  
  # Parse formula and extract variables
  mf <- model.frame(formula, data = data, na.action = na.pass)
  depvar <- names(mf)[1]
  indepvars <- names(mf)[-1]
  
  if (length(indepvars) == 0) {
    stop("At least one independent variable is required")
  }
  
  # Get panel structure
  data <- data[order(data[[id]], data[[time]]), ]
  panels <- unique(data[[id]])
  N <- length(panels)
  
  time_vals <- unique(data[[time]])
  time_vals <- sort(time_vals)
  Tmin <- min(time_vals)
  Tmax <- max(time_vals)
  Tobs <- Tmax - Tmin + 1
  
  # Validate balanced panel
  for (i in seq_along(panels)) {
    panel_times <- data[[time]][data[[id]] == panels[i]]
    if (length(panel_times) != Tobs) {
      stop("Unbalanced panel detected for unit ", panels[i],
           ". Expected ", Tobs, " observations, found ", length(panel_times))
    }
  }
  
  # ---- Get empirical moments (from GAUSS Monte Carlo) ----
  moments <- .get_moments(model_num)
  mean_t <- moments$mean
  var_t <- moments$var
  moments_depend_on_lambda <- moments$lambda_dependent
  
  # ---- Step 1: Iterative Factor-Break Estimation ----
  factor_result <- .factcoint_iter(
    data = data,
    depvar = depvar,
    indepvars = indepvars,
    id = id,
    time = time,
    model = model_num,
    max_factors = max_factors,
    trim = trim,
    max_iter = max_iter,
    tolerance = tolerance
  )
  
  n_factors <- factor_result$n_factors
  n_iters <- factor_result$iterations
  final_ssr <- factor_result$ssr
  Dres <- factor_result$Dres  # First-differenced idiosyncratic residuals (T-1 x N)
  Fhat <- factor_result$Fhat  # Estimated factors (T-1 x n_factors)
  breaks <- factor_result$breaks  # Break dates (N x 1)
  
  # ---- Model 5: Lambda-dependent moments ----
  if (model_num == 5 && moments_depend_on_lambda) {
    # Use the common break point
    lambda <- breaks[1] / Tobs
    new_moments <- .get_lambda_moments(lambda)
    mean_t <- new_moments$mean
    var_t <- new_moments$var
  }
  
  # ---- Step 2: Individual ADF Tests ----
  adf_stats <- numeric(N)
  lag_orders <- integer(N)
  
  method_num <- if (lag_method == "auto") 1 else 0
  
  for (i in seq_len(N)) {
    # Cumulate first-differenced residuals to get levels
    # GAUSS: e = cumsumc(De)
    resid_levels <- cumsum(Dres[, i])
    
    # ADF test on levels
    adf_result <- .adfrc(resid_levels, method = method_num, p_max = max_lag)
    
    adf_stats[i] <- adf_result$t_adf
    lag_orders[i] <- adf_result$p_sel
  }
  
  # ---- Step 3: Panel Test Statistic ----
  # GAUSS formula:
  # test_t = (N^(-1/2)*sumc(m_adf) - mean_t*sqrt(N)) / sqrt(var_t)
  # = sqrt(N)*(tbar - mean_t) / sqrt(var_t)
  
  valid_adf <- adf_stats[!is.na(adf_stats)]
  n_valid <- length(valid_adf)
  
  if (n_valid == 0) {
    stop("No valid ADF statistics computed")
  }
  
  tbar <- mean(valid_adf)
  Z_t <- sqrt(n_valid) * (tbar - mean_t) / sqrt(var_t)
  p_value <- stats::pnorm(Z_t)  # One-sided (left tail)
  
  # Rejection rate at 5% (GAUSS: meanc(m_adf .lt -1.95))
  n_reject <- sum(valid_adf < -1.95)
  reject_pct <- 100 * n_reject / n_valid
  
  # ---- Step 4: MQ Test for Stochastic Trends ----
  n_trends <- 0
  MQ_np <- NA_real_
  n_trends_p <- 0
  MQ_p <- NA_real_
  Fhat_cumul <- NULL
  
  if (n_factors > 0) {
    # Cumulate Fhat from first diffs to levels
    Fhat_cumul <- apply(Fhat, 2, cumsum)
    
    # Non-parametric MQ test
    mq_result_np <- .mq_test(Fhat_cumul, model_num, N, parametric = FALSE)
    MQ_np <- mq_result_np$MQ
    n_trends <- mq_result_np$n_trends
    
    # Parametric MQ test
    mq_result_p <- .mq_test(Fhat_cumul, model_num, N, parametric = TRUE)
    MQ_p <- mq_result_p$MQ
    n_trends_p <- mq_result_p$n_trends
  }
  
  # ---- Build result object ----
  result <- list(
    Z_t = Z_t,
    p_value = p_value,
    tbar = tbar,
    mean_t = mean_t,
    var_t = var_t,
    N = N,
    T = Tobs,
    n_factors = n_factors,
    n_trends = n_trends,
    n_trends_p = n_trends_p,
    MQ_np = MQ_np,
    MQ_p = MQ_p,
    iterations = n_iters,
    reject_pct = reject_pct,
    adf_stats = adf_stats,
    lag_orders = lag_orders,
    breaks = breaks,
    factors = Fhat_cumul,
    model = model_label,
    model_num = model_num,
    depvar = depvar,
    indepvars = indepvars,
    lag_method = lag_method,
    panels = panels,
    Tmin = Tmin,
    call = cl
  )
  
  class(result) <- "xtbreakcoint"
  return(result)
}


#' Print Method for xtbreakcoint Objects
#'
#' @param x An object of class `"xtbreakcoint"`.
#' @param ... Additional arguments (ignored).
#'
#' @return Invisibly returns `x`.
#'
#' @export
print.xtbreakcoint <- function(x, ...) {
  
  cat("\n")
  cat(strrep("=", 70), "\n")
  cat("Panel Cointegration Test with Structural Breaks\n")
  cat(strrep("=", 70), "\n")
  cat("Reference: Banerjee & Carrion-i-Silvestre (2015, JAE)\n")
  cat(strrep("-", 70), "\n")
  cat("Dep. var:    ", x$depvar, "\n")
  cat("Indep. vars: ", paste(x$indepvars, collapse = ", "), "\n")
  cat("Model:       ", x$model, "\n")
  cat("N (panels):  ", x$N, "\n")
  cat("T (periods): ", x$T, "\n")
  cat("Factors:     ", x$n_factors, "\n")
  cat("Lag method:  ", x$lag_method, "\n")
  cat(strrep("-", 70), "\n\n")
  
  cat("H0: No cointegration (unit root in residuals)\n")
  cat("H1: Cointegration exists\n\n")
  
  cat(sprintf("  Average t-ADF (tbar):      %9.4f\n", x$tbar))
  cat(sprintf("  E[t] under H0:             %9.4f\n", x$mean_t))
  cat(sprintf("  Var[t] under H0:           %9.4f\n", x$var_t))
  cat(strrep("-", 50), "\n")
  cat(sprintf("  Panel Z_t statistic:       %9.4f\n", x$Z_t))
  cat(sprintf("  p-value (one-sided):       %9.4f\n", x$p_value))
  cat(strrep("-", 50), "\n\n")
  
  # Decision
  if (x$p_value < 0.01) {
    cat("Decision: Reject H0 at 1% -- strong evidence of cointegration\n")
  } else if (x$p_value < 0.05) {
    cat("Decision: Reject H0 at 5% -- evidence of cointegration\n")
  } else if (x$p_value < 0.10) {
    cat("Decision: Reject H0 at 10% -- weak evidence of cointegration\n")
  } else {
    cat("Decision: Fail to reject H0 -- no evidence of cointegration\n")
  }
  
  cat(sprintf("\nIndividual rejection rate (5%%): %.1f%% (%d/%d units)\n",
              x$reject_pct, round(x$reject_pct * x$N / 100), x$N))
  
  if (x$n_factors > 0) {
    cat(sprintf("Common factors detected: %d\n", x$n_factors))
    cat(sprintf("Stochastic trends (MQ test): %d (non-parametric), %d (parametric)\n",
                x$n_trends, x$n_trends_p))
  }
  
  cat(strrep("=", 70), "\n\n")
  
  invisible(x)
}


#' Summary Method for xtbreakcoint Objects
#'
#' @param object An object of class `"xtbreakcoint"`.
#' @param ... Additional arguments (ignored).
#'
#' @return Invisibly returns the object.
#'
#' @export
summary.xtbreakcoint <- function(object, ...) {
  
  x <- object
  
  # Print main results
  print(x)
  
  # Individual unit results
  cat("Individual ADF Test Results:\n")
  cat(strrep("-", 60), "\n")
  cat(sprintf("%12s  %10s  %8s  %10s\n", "Panel", "t-ADF", "Lag", "Break"))
  cat(strrep("-", 60), "\n")
  
  for (i in seq_along(x$panels)) {
    brk <- x$breaks[i]
    brk_str <- if (brk > 0) as.character(brk + x$Tmin) else "---"
    cat(sprintf("%12s  %10.4f  %8d  %10s\n",
                as.character(x$panels[i]),
                x$adf_stats[i],
                x$lag_orders[i],
                brk_str))
  }
  cat(strrep("-", 60), "\n\n")
  
  # MQ test results if factors detected
  if (x$n_factors > 0) {
    cat("MQ Test for Common Stochastic Trends (Bai & Ng, 2004):\n")
    cat(strrep("-", 60), "\n")
    cat(sprintf("  Non-parametric MQ: %10.4f  (%d/%d stochastic trends)\n",
                x$MQ_np, x$n_trends, x$n_factors))
    cat(sprintf("  Parametric MQ:     %10.4f  (%d/%d stochastic trends)\n",
                x$MQ_p, x$n_trends_p, x$n_factors))
    cat(strrep("-", 60), "\n\n")
  }
  
  invisible(x)
}
