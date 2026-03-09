#' @keywords internal
"_PACKAGE"

#' xtbreakcoint: Panel Cointegration Tests with Structural Breaks
#'
#' @description
#' The xtbreakcoint package implements panel cointegration tests that allow for
#' structural breaks and cross-section dependence, following the methodology of
#' Banerjee and Carrion-i-Silvestre (2015). The package provides a comprehensive
#' framework for testing cointegration in panel data when standard assumptions
#' may be violated.
#'
#' @section Main Function:
#' \describe{
#'   \item{\code{\link{xtbreakcoint}}}{Panel cointegration test with structural
#'     breaks and common factors}
#' }
#'
#' @section Key Features:
#' \itemize{
#'   \item Accounts for cross-section dependence through common factors
#'   \item Allows for structural breaks in the cointegrating relationship
#'   \item Supports five model specifications with varying deterministic
#'     components
#'   \item Provides individual and panel test statistics
#'   \item Includes the Bai & Ng (2004) MQ test for stochastic trends
#' }
#'
#' @section Model Specifications:
#' \describe{
#'   \item{Model 1 (constant)}{Constant only in the deterministic component}
#'   \item{Model 2 (trend)}{Constant plus linear trend}
#'   \item{Model 3 (levelshift)}{Constant plus level shift at the break}
#'   \item{Model 4 (trendshift)}{Constant, trend, plus level shift (default)}
#'   \item{Model 5 (regimeshift)}{Constant, trend, level shift, plus slope
#'     shift}
#' }
#'
#' @section Methodology:
#' The test proceeds in four steps:
#' \enumerate{
#'   \item Iterative estimation of common factors and structural break dates
#'   \item ADF tests on defactored (idiosyncratic) residuals for each unit
#'   \item Construction of a standardized panel test statistic
#'   \item MQ test to determine whether common factors are I(1) or I(0)
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
#' Bai, J., & Ng, S. (2002). Determining the number of factors in approximate
#' factor models. \emph{Econometrica}, 70(1), 191-221.
#' \doi{10.1111/1468-0262.00273}
#'
#' @docType package
#' @name xtbreakcoint-package
NULL
