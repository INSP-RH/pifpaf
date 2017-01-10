pif: **Potential Impact Fraction**
================

Introduction
------------

The *pif* package works for the estimation of the [Potential Impact Fraction](http://www.who.int/publications/cra/chapters/volume2/2129-2140.pdf) as well as the [**Population Attributable Fraction**](http://www.who.int/healthinfo/global_burden_disease/metrics_paf/en/) from cross-sectional survey data.

Installing
----------

*pif* is still at its developmental stage. You need to install it from Github:

``` r
install.packages(devtools)
devtools::install_github("INSP-RH/pif", build_vignettes = TRUE)
```

You can also run our app without installation by copying this text into the R console:

``` r
install.packages(shiny)
shiny::runGitHub("pif", "INSP-RH", subdir = "inst/shiny-examples/PIFapp")
```

if you are only interested in the app please jump directly to [this section](#shiny-app).

Examples
--------

### **Population Attributable Fraction**

Consider that exposure to a certain phenomena is modeled by a random variable X which is normally distributed with mean μ and standard deviation σ. In addition, consider the relative risk function associated to the exposure to be:

![](README_files/figure-markdown_github/unnamed-chunk-2-1.png) where X denotes the exposure and θ the length (in time) exposed.

This function can be coded as:

``` r
rr <- function(x, theta){exp(theta*x)}
```

The **Population Attributable Fraction** is given by [the equation](https://github.com/INSP-RH/pif/blob/master/Theoretical/Worked_formulas.pdf):

![](README_files/figure-markdown_github/unnamed-chunk-4-1.png)

In particular, if we let θ = 1/9, μ = 9 and σ = 3:

![](README_files/figure-markdown_github/unnamed-chunk-5-1.png)

#### Using the package

Let θ be estimated by θ' = 0.11. Consider a random sample of the X's:

``` r
set.seed(242)
X <- rnorm(1000, 9, 3)
```

We use our package to estimate the population atributable fraction:

``` r
library(pif)
paf(X, 0.11, rr)
```

    ## [1] 0.6548373

The default can be changed to using kernels:

``` r
paf(X, 0.11, rr, method = "kernel")
```

    ## [1] 0.6557163

where the kernel-type, adjustment and bandwidth can be changed:

``` r
paf(X, 0.11, rr, method = "kernel", ktype = "cosine", adjust = 0.6)
```

    ## [1] 0.6551543

Confidence intervals can be obtained if the variance of θ' is, say 0.0025:

``` r
paf.confidence(X, 0.11, 0.0025, rr = rr)
```

    ##          Lower Point_Estimate          Upper 
    ##     0.01320037     0.65483733     0.87926904

Several methods have been implemented to give the variance. In particular, because `exp(θ X)` is convex, we can use the `one2one` method by specifying the lower and upper bounds of θ's confidence interval:

``` r
paf.confidence(X, 0.11, thetalow = 0.09, thetaup = 0.13, rr = rr, method = "one2one")
```

    ## Lower.Lower       Point Upper.Upper 
    ##   0.3538432   0.6548373   0.8420717

A method for estimating the PAF when only the estimators of the exposure mean and variance are known is given by:

``` r
pif.approximate(mean(X), var(X), 0.11, rr)
```

    ## [1] 0.6543426

### **Potential Impact Fraction**

Consider the same problem as before. We are now interested in estimating the effects of a counterfactual. As an example, consider a transformation of X given by: `T(X) = aX + b`. In this case, the **Potential Impact Fraction** is given by [the equation](https://github.com/INSP-RH/pif/blob/master/Theoretical/Worked_formulas.pdf):

![](README_files/figure-markdown_github/unnamed-chunk-13-1.png)

Assuming the same values as above and considering `a = 1/2` and `b = -1` we have:

![](README_files/figure-markdown_github/unnamed-chunk-14-1.png)

#### Using the package

We use our package to estimate the **Potential Impact Fraction**. First: we code the counterfactual function:

``` r
cft <- function(X){0.5*X - 1}
```

Then we estimate the PIF:

``` r
pif(X, 0.11, rr, cft)
```

    ## [1] 0.480594

which can also be done with a kernel estimator:

``` r
pif(X, 0.11, rr, cft, method = "kernel", ktype = "gaussian")
```

    ## [1] 0.4815858

### Sensitivity analysis

#### Plots

The command `pif.plot` allows us to analyze how the potential impact fraction varies as the values of θ change:

``` r
pif.plot(X, 0, 0.3, rr, cft)
```

![](README_files/figure-markdown_github/unnamed-chunk-18-1.png)

#### Sensitivity function

The command `pif.sensitivity` allows us to analyze how our estimates for the potential impact fraction would vary had we excluded some part of the exposure sample.

``` r
pif.sensitivity(X, 0.11, rr, cft)
```

![](README_files/figure-markdown_github/unnamed-chunk-19-1.png)

#### Counterfactual Heatmap

We can also generate a heatmap showing how distinct counterfactual assumptions result in different impact fractions

``` r
pif.counterfactual.heatmap(X, 0.11, rr)
```

![](README_files/figure-markdown_github/unnamed-chunk-20-1.png)

Shiny App
---------

If you have installed our package the function `PIFApp` runs the shiny-app in your computer:

``` r
PIFApp()
```

This will open a new window with the App:

<img alt = "ShinyApp" src = "README_files/shinyapp.png">

What is missing?
----------------

We are working for decent confidence intervals for the Impact Fraction.

In addition, confidence intervals for the approximate PAF need be coded.

The User's Interface `ui.R` of the Shiny-App will be updated

Please feel free to contribute to the project.

Contributor Code of Conduct
---------------------------

As contributors and maintainers of this project, we pledge to respect all people who contribute through reporting issues, posting feature requests, updating documentation, submitting pull requests or patches, and other activities.

We are committed to making participation in this project a harassment-free experience for everyone, regardless of level of experience, gender, gender identity and expression, sexual orientation, disability, personal appearance, body size, race, ethnicity, age, or religion.

Examples of unacceptable behavior by participants include the use of sexual language or imagery, derogatory comments or personal attacks, trolling, public or private harassment, insults, or other unprofessional conduct.

Project maintainers have the right and responsibility to remove, edit, or reject comments, commits, code, wiki edits, issues, and other contributions that are not aligned to this Code of Conduct. Project maintainers who do not follow the Code of Conduct may be removed from the project team.

Instances of abusive, harassing, or otherwise unacceptable behavior may be reported by opening an issue or contacting one or more of the project maintainers.

This Code of Conduct is adapted from the Contributor Covenant, version 1.0.0, available from <http://contributor-covenant.org/version/1/0/0/>

Licence
-------

This package is free and open source software, licensed under [GPL-3](https://www.gnu.org/licenses/gpl-3.0.html).

If you use this package please don't forget to cite it.

Authors
-------

-   Rodrigo Zepeda-Tello <rodrigo.zepeda@insp.mx>
-   Dalia Camacho-García-Formentí <daliaf172@gmail.com>
