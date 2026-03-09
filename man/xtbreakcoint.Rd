\name{xtbreakcoint}
\alias{xtbreakcoint}
\title{Panel Cointegration Test with Structural Breaks}
\description{
Implements the panel cointegration test of Banerjee and Carrion-i-Silvestre
(2015), allowing for structural breaks and cross-section dependence through
common factors.
}
\usage{
xtbreakcoint(
  formula,
  data,
  id,

  time,
  model = "trendshift",
  max_factors = 5,
  max_lag = 4,

  lag_method = c("auto", "fixed"),
  trim = 0.15,
  max_iter = 20,
  tolerance = 0.001
)
}
\arguments{
\item{formula}{A formula of the form \code{y ~ x1 + x2 + ...} specifying the
cointegrating relationship.}

\item{data}{A data frame containing panel data with columns for the panel
identifier, time identifier, and all variables in the formula.}

\item{id}{Character string naming the panel (cross-section) identifier.}

\item{time}{Character string naming the time identifier.}

\item{model}{Model specification for deterministic components. One of:
\code{"constant"} (1), \code{"trend"} (2), \code{"levelshift"} (3),
\code{"trendshift"} (4, default), or \code{"regimeshift"} (5). Can also
specify as integer 1-5.}

\item{max_factors}{Maximum number of common factors to estimate (default: 5).
Set to 0 to skip factor estimation.}

\item{max_lag}{Maximum lag order for ADF tests (default: 4).}

\item{lag_method}{Method for selecting ADF lag order: \code{"auto"} for
automatic selection via BIC (default) or \code{"fixed"} to use
\code{max_lag}.}

\item{trim}{Trimming parameter for break estimation, proportion of sample
excluded from endpoints (default: 0.15). Must be in (0, 0.5).}

\item{max_iter}{Maximum iterations for factor-break estimation (default: 20).}

\item{tolerance}{Convergence tolerance for iterative estimation (default:
0.001).}
}
\value{
An object of class \code{"xtbreakcoint"} containing:
\describe{
  \item{Z_t}{Panel test statistic (standard normal under H0)}
  \item{p_value}{One-sided p-value for Z_t}
  \item{tbar}{Average of individual ADF t-statistics}
  \item{mean_t}{Expected value of t-statistic under H0}
  \item{var_t}{Variance of t-statistic under H0}
  \item{N}{Number of cross-section units}
  \item{T}{Number of time periods}
  \item{n_factors}{Number of estimated common factors}
  \item{n_trends}{Number of stochastic trends (from MQ test)}
  \item{MQ_np}{Non-parametric MQ test statistic}
  \item{MQ_p}{Parametric MQ test statistic}
  \item{iterations}{Number of iterations for factor-break convergence}
  \item{reject_pct}{Percentage of individual units rejecting H0 at 5\%}
  \item{adf_stats}{Vector of individual ADF t-statistics}
  \item{lag_orders}{Vector of selected lag orders for each unit}
  \item{breaks}{Vector of estimated break dates (periods from start)}
  \item{factors}{Matrix of estimated common factors (T-1 x n_factors)}
  \item{model}{Model specification used}
  \item{call}{The matched call}
}
}
\details{
This function tests for panel cointegration in the presence of structural
breaks and cross-sectional dependence. The methodology follows these steps:

\enumerate{
  \item \strong{Iterative Factor-Break Estimation}: Jointly estimates common
    factors and individual break dates using an iterative procedure.
  \item \strong{Individual ADF Tests}: Applies ADF tests to the defactored
    (idiosyncratic) residuals for each cross-section unit.
  \item \strong{Panel Test Statistic}: Combines individual ADF statistics
    into a standardized panel statistic using Monte Carlo moments.
  \item \strong{MQ Test}: If factors are detected, tests whether they are
    stationary or represent stochastic trends (Bai & Ng, 2004).
}
}
\section{Model Specifications}{
\describe{
  \item{1 - constant}{Constant only}
  \item{2 - trend}{Constant plus linear trend}
  \item{3 - levelshift}{Constant plus level shift at break}
  \item{4 - trendshift}{Constant, trend, plus level shift (default)}
  \item{5 - regimeshift}{Constant, trend, level shift, plus slope shift}
}
}

\references{
Banerjee, A., & Carrion-i-Silvestre, J. L. (2015). Cointegration in panel
data with structural breaks and cross-section dependence. \emph{Journal of
Applied Econometrics}, 30(1), 1-22. \doi{10.1002/jae.2348}

Bai, J., & Ng, S. (2004). A PANIC attack on unit roots and cointegration.
\emph{Econometrica}, 72(4), 1127-1177. \doi{10.1111/j.1468-0262.2004.00528.x}
}
\examples{
# Generate example panel data
set.seed(42)
N <- 10  # panels
T <- 50  # time periods

# Create cointegrated data with a structural break
panel_data <- data.frame(
  id = rep(1:N, each = T),
  time = rep(1:T, N),
  y = NA,
  x = NA
)

for (i in 1:N) {
  idx <- panel_data$id == i
  x <- cumsum(rnorm(T))
  u <- rnorm(T, sd = 0.5)
  # Cointegrating relationship with break at t=25
  beta <- ifelse(1:T <= 25, 1, 1.5)
  y <- 1 + beta * x + u
  panel_data$x[idx] <- x
  panel_data$y[idx] <- y
}

# Test for cointegration
result <- xtbreakcoint(y ~ x, data = panel_data, id = "id", time = "time")
print(result)
}
\author{Muhammad Alkhalaf}

Based on original GAUSS code by A. Banerjee and J.L. Carrion-i-Silvestre.
}
