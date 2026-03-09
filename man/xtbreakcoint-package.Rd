\name{xtbreakcoint-package}
\alias{xtbreakcoint-package}
\docType{package}
\title{Panel Cointegration Tests with Structural Breaks}
\description{
Implements panel cointegration tests allowing for structural breaks and
cross-section dependence following the methodology of Banerjee and
Carrion-i-Silvestre (2015). The package provides iterative factor-break
estimation, individual ADF tests on defactored residuals, standardized panel
test statistics, and the Bai and Ng (2004) MQ test for identifying common
stochastic trends.
}
\details{
The main function is \code{\link{xtbreakcoint}}, which performs the panel
cointegration test.

\strong{Key Features:}
\itemize{
  \item Accounts for cross-section dependence through common factors
  \item Allows for structural breaks in the cointegrating relationship
  \item Supports five model specifications with varying deterministic
    components
  \item Provides individual and panel test statistics
  \item Includes the Bai & Ng (2004) MQ test for stochastic trends
}

\strong{Model Specifications:}
\describe{
  \item{Model 1 (constant)}{Constant only}
  \item{Model 2 (trend)}{Constant plus linear trend}
  \item{Model 3 (levelshift)}{Constant plus level shift at break}
  \item{Model 4 (trendshift)}{Constant, trend, plus level shift (default)}
  \item{Model 5 (regimeshift)}{Constant, trend, level shift, plus slope shift}
}
}
\references{
Banerjee, A., & Carrion-i-Silvestre, J. L. (2015). Cointegration in panel
data with structural breaks and cross-section dependence. \emph{Journal of
Applied Econometrics}, 30(1), 1-22. \doi{10.1002/jae.2348}

Bai, J., & Ng, S. (2004). A PANIC attack on unit roots and cointegration.
\emph{Econometrica}, 72(4), 1127-1177. \doi{10.1111/j.1468-0262.2004.00528.x}

Bai, J., & Ng, S. (2002). Determining the number of factors in approximate
factor models. \emph{Econometrica}, 70(1), 191-221.
\doi{10.1111/1468-0262.00273}
}
\author{
Dr Merwan Roudane \email{merwanroudane920@gmail.com}

Based on original GAUSS code by A. Banerjee and J.L. Carrion-i-Silvestre.
}
\keyword{package}
