pif: **Potential Impact Fraction**
================

Installing the package
----------------------

The R package `pif` was developed to calculate the Population Attributable Fraction (PAF) and the Population Impact Fraction (PIF) via the empirical, kernel, and approximate methods. Confidence intervals for the PAF and the PIF using different approaches have been programmed. Along with the confidence intervals several sensitivity analysis can be conducted.

We suggest installing the package from Github to get the latest version by using the following code:

``` r
install.packages("devtools")
devtools::install_github("INSP-RH/pif")
```

Recall that, once installed, to use a package in R you need to call the associated `library`:

``` r
library(pif)
```

In the following sections we show how the package can be used to estimate both `paf`and `pif`

Potential Attributable Fraction (PAF)
-------------------------------------

Let X denote exposure values to something. And define the exponential relative risk function RR(X|ø)=exp(ø X), where the parameter ø is estimated to be 0.11.

``` r
rr       <- function(x, theta){exp(theta*x)}
thetahat <- 0.11                            
X        <- data.frame(Exposure = rbeta(1000, 1, 3))
```

We use the `paf` function to estimate the population atributable fraction:

``` r
paf(X = X, thetahat = thetahat, rr = rr)
```

    ## [1] 0.0285131

The PAF calculated above considers the sampling weights of the exposure values to be equal. In case the observed exposure values have different sampling weights, the PAF is calculated by:

``` r
weights <- c(rep(1/1500, 500), rep(2/1500, 500))
paf(X = X, thetahat = thetahat, rr = rr, weights = weights)
```

    ## [1] 0.02853325

The default method is empirical, but it can be changed to kernel method when the exposure values for one observation are unidimensional and there are no covariates.

``` r
paf(X = X, thetahat = thetahat, rr = rr, method = "kernel")
```

    ## [1] 0.02943677

Confidence intervals for the PAF can be calculated if the variance of the estimator of ø is known. For example if the variance of the estimator of ø is 0.0025:

``` r
thetavar <- 0.0025
paf.confidence(X = X, thetahat = thetahat, thetavar = thetavar, rr = rr)
```

    ##           Lower_CI     Point_Estimate           Upper_CI 
    ##       0.0042985131       0.0285130990       0.0537475882 
    ## Estimated_Variance 
    ##       0.0001663288

The default method is `bootstrap`; however, other confidence interval methods are available:

``` r
paf.confidence(X = X, thetahat = thetahat, thetavar = thetavar, rr = rr,
               confidence_method = "inverse")
```

    ##       Lower_CI Point_Estimate       Upper_CI 
    ##    0.003213858    0.028513099    0.053170225

The previous examples consider a random sample of exposure values is available, however we developed a method for estimating the PAF when only mean and variance of the exposure values are known.

``` r
Xmean <- data.frame(mean(X[,"Exposure"]))
Xvar  <- var(X)
paf(X = Xmean, Xvar = Xvar, thetahat = thetahat, rr = rr, method = "approximate")
```

    ## [1] 0.02851213

We remind the reader that if the whole exposure X is available the approximate method is the last resource.

### Potential Impact Fraction (PIF)

Consider the same problem as before. We are now interested in estimating the effects of a counterfactual scenario, which is not necessarily the no-exposure counterfactual. As an example, consider a transformation of X given by: aX + b. Consider a = 0.5 and b = -1, to estimate the PIF first we code the counterfactual function:

``` r
cft <- function(X){0.5*X - 1}
```

Then we estimate the PIF via the empirical method:

``` r
pif(X = X, thetahat = thetahat, rr = rr, cft = cft)
```

    ## [1] 0.1170813

The PIF can also be calculated by estimating a distribution for exposure values using a kernel:

``` r
pif(X = X, thetahat = thetahat, rr = rr, cft = cft, 
    method = "kernel")
```

    ## [1] 0.1170839

confidence intervals are available through `pif.confidence`:

``` r
pif.confidence(X = X, thetahat = thetahat, rr = rr, cft = cft, 
               method = "kernel", thetavar = thetavar)
```

    ##           Lower_CI     Point_Estimate           Upper_CI 
    ##        0.020677535        0.117083916        0.216560276 
    ## Estimated_Variance 
    ##        0.002495226

Additional functions in the package include sensitivity analysis plots `pif.sensitivity`, `paf.sensitivity`, `pif.plot`, `paf.plot`, `pif.heatmap` and `counterfactual.plot` for display of the counterfactual function. An exploration of this functions as well as additional examples of usage and utilization of advanced options can be found in the package's vignettes:

``` r
browseVignettes("pif")
```
