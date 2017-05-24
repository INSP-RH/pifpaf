pifpaf: **Potential Impact Fraction and Population Attributable Fraction Estimation**
================

Installing the package
----------------------

The R package `pifpaf` was developed to calculate the Population Attributable Fraction (PAF) and the Population Impact Fraction (PIF) via the empirical, kernel, and approximate methods. Confidence intervals for the PAF and the PIF using different approaches have been programmed. Along with the confidence intervals several sensitivity analysis can be conducted.

We suggest installing the package from Github to get the latest version by using the following code:

``` r
if(!require(devtools)){install.packages("devtools")}
devtools::install_github("INSP-RH/pifpaf")
```

Recall that, once installed, to use a package in R you need to call the associated `library`:

``` r
library(pifpaf)
```

In the following sections we show how the package can be used to estimate both `paf`and `pif`

Potential Attributable Fraction (PAF)
-------------------------------------

Let X denote exposure values to something. And define the exponential relative risk function RR(X;ø)=exp(ø X), where the parameter ø is estimated to be 0.11.

``` r
rr       <- function(x, theta){exp(theta*x)}
thetahat <- 0.11                            
X        <- data.frame(Exposure = rbeta(1000, 1, 3))
```

We use the `paf` function to estimate the population atributable fraction:

``` r
paf(X = X, thetahat = thetahat, rr = rr)
```

    ## [1] 0.02776485

The PAF calculated above considers the sampling weights of the exposure values to be equal. In case the observed exposure values have different sampling weights, the PAF is calculated by:

``` r
weights <- c(rep(1/1500, 500), rep(2/1500, 500))
paf(X = X, thetahat = thetahat, rr = rr, weights = weights)
```

    ## [1] 0.02795618

The default method is empirical, but it can be changed to kernel method when the exposure values for one observation are unidimensional and there are no covariates.

``` r
paf(X = X, thetahat = thetahat, rr = rr, method = "kernel")
```

    ## [1] 0.02868922

Confidence intervals for the PAF can be calculated if the variance of the estimator of ø is known. For example if the variance of the estimator of ø is 0.0025:

``` r
thetavar <- 0.0025
paf.confidence(X = X, thetahat = thetahat, thetavar = thetavar, rr = rr)
```

    ##           Lower_CI     Point_Estimate           Upper_CI 
    ##       0.0040516572       0.0277648519       0.0535040581 
    ## Estimated_Variance 
    ##       0.0001574953

The default method is `bootstrap`; however, other confidence interval methods are available:

``` r
paf.confidence(X = X, thetahat = thetahat, thetavar = thetavar, rr = rr,
               confidence_method = "inverse")
```

    ##       Lower_CI Point_Estimate       Upper_CI 
    ##    0.003122144    0.027764852    0.051798395

The previous examples consider a random sample of exposure values is available, however we developed a method for estimating the PAF when only mean and variance of the exposure values are known.

``` r
Xmean <- data.frame(mean(X[,"Exposure"]))
Xvar  <- var(X)
paf(X = Xmean, Xvar = Xvar, thetahat = thetahat, rr = rr, method = "approximate")
```

    ## [1] 0.0277637

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

    ## [1] 0.1167395

The PIF can also be calculated by estimating a distribution for exposure values using a kernel:

``` r
pif(X = X, thetahat = thetahat, rr = rr, cft = cft, 
    method = "kernel")
```

    ## [1] 0.1167428

confidence intervals are available through `pif.confidence`:

``` r
pif.confidence(X = X, thetahat = thetahat, rr = rr, cft = cft, 
               method = "kernel", thetavar = thetavar)
```

    ##           Lower_CI     Point_Estimate           Upper_CI 
    ##        0.024585690        0.116742828        0.212889259 
    ## Estimated_Variance 
    ##        0.002492649

Additional functions in the package include sensitivity analysis plots `pif.sensitivity`, `paf.sensitivity`, `pif.plot`, `paf.plot`, `pif.heatmap` and `counterfactual.plot` for display of the counterfactual function. An exploration of this functions as well as additional examples of usage and utilization of advanced options can be found in the [package's vignettes](vignettes/Introduction_to_pifpaf_package.html):

``` r
browseVignettes("pif")
```
