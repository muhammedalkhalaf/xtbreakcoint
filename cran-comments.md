## xtbreakcoint 1.0.6

This release corrects the computations below; the 1.0.5 submission (reference metadata only) should be discarded in favour of this one.

* Bug fix: the estimation did not follow Banerjee and Carrion-i-Silvestre (2015). It pooled the units in one first-differenced regression, dated the breaks by regressing residuals on a step dummy, chose the number of factors with IC2, selected ADF lags by BIC and used different MQ computations. The engine is now a port of the authors' GAUSS replication code (factcoint_iter, factcoint, ADFRC and MQ_test): unit-by-unit cointegrating regressions in first differences, break dates by minimum SSR over the central 70% of the sample (a common break for model 5), the Bai and Ng (2002) IC1 criterion, the iteration between factors and breaks, ADF regressions without deterministic terms with general-to-specific lag selection at |t| = 1.645, and the MQ tests of Bai and Ng (2004) on the detrended factors.
* The GAUSS code dates the model 4 break by minimising the squared coefficient vector rather than the sum of squared residuals; the package uses the sum of squared residuals for every model, as described in the paper.
* Without factors, the estimated break dates agree with the Stata command xtbreakcoint (SSC) on the same data.
* Break dates in `summary()` are now reported as the last period of the first regime.

## Test environments

* Ubuntu 24.04, R 4.3.3 and R-devel, R CMD check --as-cran

## R CMD check results

0 errors | 0 warnings | 0 notes
