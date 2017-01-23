pif: **Potential Impact Fraction**
================

Installing the package
----------------------

The R package `pif` was developed to calculate the Population Attributable Fraction (PAF) and the Population Impact Fraction (PIF) via the empirical, kernel, and approximate methods. Confidence intervals for the PAF and the PIF using different approaches have been programmed. Along with the confidence intervals several sensitivity analysis can be conducted.

We suggest installing the package from Github to get the latest version by using the following code:

``` r
install.packages("devtools")
devtools::install_github("INSP-RH/pif", build_vignettes = TRUE)
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
X        <- rnorm(1000, 9, 3)
```

We use the `paf` function to estimate the population atributable fraction:

``` r
paf(X = X, thetahat = thetahat, rr = rr)
```

    ## Warning in check.exposure(.X): Some exposure values are less than zero,
    ## verify this is correct.

    ## [1] 0.6426382

The PAF calculated above considers the sampling weights of the exposure values to be equal. In case the observed exposure values have different sampling weights, the PAF is calculated by:

``` r
weights <- c(rep(1/1500, 500), rep(2/1500, 500))
paf(X = X, thetahat = thetahat, rr = rr, weights = weights)
```

    ## Warning in check.exposure(.X): Some exposure values are less than zero,
    ## verify this is correct.

    ## [1] 0.6433246

The default method is empirical, but it can be changed to kernel method when the exposure values for one observation are unidimensional and there are no covariates.

``` r
paf(X = X, thetahat = thetahat, rr = rr, method = "kernel")
```

    ## Warning in check.exposure(.X): Some exposure values are less than zero,
    ## verify this is correct.

    ## [1] 0.6441464

Confidence intervals for the PAF can be calculated if the variance of the estimator of ø is known. For example if the variance of the estimator of ø is 0.0025:

``` r
thetavar <- 0.0025
paf.confidence(X = X, thetahat = thetahat, thetavar = thetavar, rr = rr)
```

    ##          Lower Point_Estimate          Upper 
    ##      0.2171866      0.6426382      0.8368609

The default method is `bootstrap`; however, other confidence interval methods are available:

``` r
paf.confidence(X = X, thetahat = thetahat, thetavar = thetavar, rr = rr,
               confidence_method = "inverse")
```

    ##          Lower Point_Estimate          Upper 
    ##      0.2390676      0.6426382      0.8321697

The previous examples consider a random sample of exposure values is available, however we developed a method for estimating the PAF when only mean and variance of the exposure values are known.

``` r
Xmean <- mean(X)
Xvar  <- var(X)
paf(X = Xmean, Xvar = Xvar, thetahat = thetahat, rr = rr, method = "approximate")
```

    ## [1] 0.6419899

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

    ## Warning in check.exposure(.X): Some exposure values are less than zero,
    ## verify this is correct.

    ## [1] 0.4718027

The PIF can also be calculated by estimating a distribution for exposure values using a kernel:

``` r
pif(X = X, thetahat = thetahat, rr = rr, cft = cft, 
    method = "kernel")
```

    ## Warning in check.exposure(.X): Some exposure values are less than zero,
    ## verify this is correct.

    ## [1] 0.4730906

confidence intervals are available through `pif.confidence`:

``` r
pif.confidence(X = X, thetahat = thetahat, rr = rr, cft = cft, 
               method = "kernel", thetavar = thetavar)
```

    ## Warning in check.exposure(.X): Some exposure values are less than zero,
    ## verify this is correct.

    ##           Lower_CI     Point_Estimate           Upper_CI 
    ##         0.27723664         0.47309056         0.94536599 
    ## Estimated_Variance 
    ##         0.03366958

Additional functions in the package include sensitivity analysis plots `pif.sensitivity`, `paf.sensitivity`, `pif.plot`, `paf.plot`, `pif.heatmap` and `counterfactual.plot` for display of the counterfactual function. An exploration of this functions as well as additional examples of usage and utilization of advanced options can be found in the package's vignettes:

``` r
browseVignettes("pif")
```
