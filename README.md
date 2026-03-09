# xtbreakcoint

<!-- badges: start -->
[![CRAN status](https://www.r-pkg.org/badges/version/xtbreakcoint)](https://CRAN.R-project.org/package=xtbreakcoint)
<!-- badges: end -->

Panel cointegration tests allowing for structural breaks and cross-section 
dependence, implementing the methodology of Banerjee and Carrion-i-Silvestre 
(2015).

## Installation

You can install the development version from GitHub:

```r
# install.packages("devtools")
```

## Overview

The `xtbreakcoint` package tests for panel cointegration when:

- **Structural breaks** may be present in the cointegrating relationship
- **Cross-section dependence** exists through common factors

The methodology follows four steps:

1. **Iterative Factor-Break Estimation**: Jointly estimates common factors
   and individual break dates
2. **Individual ADF Tests**: Applies ADF tests to defactored residuals
3. **Panel Test Statistic**: Combines individual statistics into a
   standardized panel statistic
4. **MQ Test**: Tests whether common factors are I(1) stochastic trends

## Example

```r
library(xtbreakcoint)

# Generate example panel data
set.seed(42)
N <- 10   # panels
T <- 50   # time periods

panel_data <- data.frame(
  id = rep(1:N, each = T),
  time = rep(1:T, N),
  y = NA,
  x = NA
)

# Create cointegrated data with a structural break at t=25
for (i in 1:N) {
  idx <- panel_data$id == i
  x <- cumsum(rnorm(T))
  u <- rnorm(T, sd = 0.5)
  beta <- ifelse(1:T <= 25, 1, 1.5)  # Break in slope
  y <- 1 + beta * x + u
  panel_data$x[idx] <- x
  panel_data$y[idx] <- y
}

# Test for cointegration
result <- xtbreakcoint(y ~ x, data = panel_data, id = "id", time = "time")
print(result)
```

## Model Specifications

| Model | Name | Deterministics |
|-------|------|----------------|
| 1 | constant | Constant only |
| 2 | trend | Constant + trend |
| 3 | levelshift | Constant + level shift |
| 4 | trendshift | Constant + trend + level shift (default) |
| 5 | regimeshift | Constant + trend + level + slope shift |

## References

- Banerjee, A., & Carrion-i-Silvestre, J. L. (2015). Cointegration in panel
  data with structural breaks and cross-section dependence. *Journal of
  Applied Econometrics*, 30(1), 1-22. 
  [doi:10.1002/jae.2348](https://doi.org/10.1002/jae.2348)

- Bai, J., & Ng, S. (2004). A PANIC attack on unit roots and cointegration.
  *Econometrica*, 72(4), 1127-1177.
  [doi:10.1111/j.1468-0262.2004.00528.x](https://doi.org/10.1111/j.1468-0262.2004.00528.x)

- Bai, J., & Ng, S. (2002). Determining the number of factors in approximate
  factor models. *Econometrica*, 70(1), 191-221.
  [doi:10.1111/1468-0262.00273](https://doi.org/10.1111/1468-0262.00273)

## Author

Based on original GAUSS code by A. Banerjee and J.L. Carrion-i-Silvestre.

## License

GPL (>= 3)
