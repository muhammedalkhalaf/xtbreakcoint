# xtbreakcoint 1.0.6

* Bug fix: the estimation did not follow Banerjee and Carrion-i-Silvestre (2015). It pooled the units in one first-differenced regression, dated the breaks by regressing residuals on a step dummy, chose the number of factors with IC2, selected ADF lags by BIC and used different MQ computations. The engine is now a port of the authors' GAUSS replication code (factcoint_iter, factcoint, ADFRC and MQ_test): unit-by-unit cointegrating regressions in first differences, break dates by minimum SSR over the central 70% of the sample (a common break for model 5), the Bai and Ng (2002) IC1 criterion, the iteration between factors and breaks, ADF regressions without deterministic terms with general-to-specific lag selection at |t| = 1.645, and the MQ tests of Bai and Ng (2004) on the detrended factors.
* The GAUSS code dates the model 4 break by minimising the squared coefficient vector rather than the sum of squared residuals; the package uses the sum of squared residuals for every model, as described in the paper.
* Without factors, the estimated break dates agree with the Stata command xtbreakcoint (SSC) on the same data.
* Break dates in `summary()` are now reported as the last period of the first regime.

# xtbreakcoint 1.0.5

* Corrected two DOIs in NEWS.md: Banerjee and Carrion-i-Silvestre (2015) is 10.1002/jae.2348 and Bai and Ng (2004) is 10.1111/j.1468-0262.2004.00528.x. No changes to code.

# xtbreakcoint 1.0.1

Initial CRAN release.

## Features

* Panel cointegration test with structural breaks following Banerjee & 
  Carrion-i-Silvestre (2015, JAE)
* Cross-section dependence handled via common factor estimation
* Five model specifications for deterministic components
* Automatic factor number selection using Bai & Ng (2002) IC criterion
* Structural break date estimation via SSR minimization
* Individual ADF tests on defactored residuals
* Standardized panel test statistic with Monte Carlo moments
* Bai & Ng (2004) MQ test for identifying common stochastic trends
* Automatic or fixed lag selection for ADF regressions

## References

* Banerjee, A., & Carrion-i-Silvestre, J. L. (2015). Cointegration in panel
  data with structural breaks and cross-section dependence. *Journal of
  Applied Econometrics*, 30(1), 1-22. doi:10.1002/jae.2348

* Bai, J., & Ng, S. (2004). A PANIC attack on unit roots and cointegration.
  *Econometrica*, 72(4), 1127-1177. doi:10.1111/j.1468-0262.2004.00528.x
