#' pifpaf: A package for estimating the Potential Impact Fraction and the
#' Population Attributable Fraction from cross-sectional data
#'
#' Uses a generalized method to estimate the 
#' Potential Impact Fraction, (PIF) and the Population Attributable Fraction (PAF) 
#' from cross-sectional data. It creates point-estimates (\code{\link{pif}}, 
#' \code{\link{paf}}), confidence intervals (\code{\link{pif.confidence}}, 
#' \code{\link{paf.confidence}}),and estimates of variance. 
#' In addition it generates plots for conducting sensitivity analysis (\code{\link{pif.sensitivity}}, 
#' \code{\link{pif.heatmap}}, \code{\link{pif.plot}}). 
#' The estimation method corresponds to  Zepeda-Tello,  Camacho-García-Formentí, et al,
#' 'Nonparametric Methods to Estimate the Potential Impact Fraction from Cross-sectional Data',
#' Unpublished manuscript. 2017.
#' This package was developed by the National Institute of Public Health of Mexico
#' under funding by Bloomberg Philanthropies. Please see the package's vignette for more information:
#' \code{browseVignettes("pifpaf")}.
#' 
#'
#' @docType package
#' @name pifpaf
NULL