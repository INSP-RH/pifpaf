% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pif_heatmap.R
\name{pif.heatmap}
\alias{pif.heatmap}
\title{Graphical Sensitivity Analysis of Potential Impact Fraction's Counterfactual}
\usage{
pif.heatmap(X, thetahat, rr, cft = function(X, a, b) {     a * X + b },
  method = "empirical", weights = rep(1/nrow(as.matrix(X)),
  nrow(as.matrix(X))), Xvar = var(X), deriv.method.args = list(),
  deriv.method = "Richardson", adjust = 1, n = 512, ktype = "gaussian",
  bw = "SJ", legendtitle = "PIF", mina = 0, maxa = 1, minb = -1,
  maxb = 0, nmesh = 10,
  title = paste0("Potential Impact Fraction (PIF) with counterfactual",
  "\\nf(X)= aX+b"), xlab = "a", ylab = "b",
  colors = rev(heat.colors(nmesh)), check_exposure = TRUE,
  check_rr = TRUE, check_integrals = TRUE)
}
\arguments{
\item{X}{Random sample (\code{data.frame}) which includes exposure 
and covariates or sample \code{mean} if \code{"approximate"} method is 
selected.}

\item{thetahat}{Asymptotically consistent or Fisher consistent estimator (\code{vector}) of \code{theta} for the Relative 
Risk function.}

\item{rr}{\code{function} for Relative Risk which uses parameter 
\code{theta}. The order of the parameters shound be \code{rr(X, theta)}.

\strong{**Optional**}}

\item{cft}{\code{function} \code{cft(X, a, b)} for counterfactual dependent on
one dimensional parameters \code{a} and \code{b}. Default counterfactual is 
affine: \code{aX + b}.}

\item{method}{Either \code{"empirical"} (default), \code{"kernel"} or 
\code{"approximate"}. For details on estimation methods see 
\code{\link{pif}}.}

\item{weights}{Normalized survey \code{weights} for the sample \code{X}.}

\item{Xvar}{Variance of exposure levels (for \code{"approximate"} 
method).}

\item{deriv.method.args}{\code{method.args} for 
\code{\link[numDeriv]{hessian}} (for \code{"approximate"} method).}

\item{deriv.method}{\code{method} for \code{\link[numDeriv]{hessian}}. 
Don't change this unless you know what you are doing (for 
\code{"approximate"} method).}

\item{adjust}{Adjust bandwith parameter (for \code{"kernel"} 
method) from \code{\link[stats]{density}}.}

\item{n}{Number of equally spaced points at which the density (for 
\code{"kernel"} method) is to be estimated (see 
\code{\link[stats]{density}}).}

\item{ktype}{\code{kernel} type:  \code{"gaussian"}, 
\code{"epanechnikov"}, \code{"rectangular"}, \code{"triangular"}, 
\code{"biweight"}, \code{"cosine"}, \code{"optcosine"} (for \code{"kernel"}
method). Additional information on kernels in \code{\link[stats]{density}}.}

\item{bw}{Smoothing bandwith parameter (for 
\code{"kernel"} method) from \code{\link[stats]{density}}. Default 
\code{"SJ"}.}

\item{legendtitle}{\code{string} title for the legend of plot.}

\item{mina}{Minimum for parameter \code{a} for the counterfactual 
(default \code{0}).}

\item{maxa}{Maximum for parameter \code{a} for the counterfactual 
(default \code{1}).}

\item{minb}{Minimum for parameter \code{b} for the counterfactual 
(default \code{-1}).}

\item{maxb}{Maximum for parameter \code{b} for the counterfactual 
(default \code{0}).}

\item{nmesh}{Number of tiles in plot (default \code{10}).}

\item{title}{\code{string} title for the plot.}

\item{xlab}{\code{string} label for the X-axis of the plot (corresponding
to "a").}

\item{ylab}{\code{string} label for the Y-axis of the plot (corresponding
to "b").}

\item{colors}{\code{vector} of colours for the heatmap.}

\item{check_exposure}{\code{boolean}  Check that exposure \code{X} is 
positive and numeric.}

\item{check_rr}{\code{boolean} Check that Relative Risk function
\code{rr} equals \code{1} when evaluated at \code{0}.}

\item{check_integrals}{\code{boolean}  Check that counterfactual \code{cft} 
and relative risk's \code{rr} expected values are well defined for this 
scenario.}
}
\value{
plotpif      \code{\link[ggplot2]{ggplot}} object plotting a heatmap
  with sensitivity analysis of the counterfactual.
}
\description{
Provides a graphical sensitivity analysis for \code{\link{pif}} by
  varying the parameters of a bivariate counterfactual function \code{cft}.
  By default it evaluates the counterfactual:
  \deqn{
  \textrm{cft}(X) = aX + b.
  }{
  cft(X) = aX + b.
  }
}
\examples{
\dontrun{
#Example 1
#------------------------------------------------------------------
X  <- data.frame(rnorm(25,3))            #Sample
rr <- function(X,theta){exp(X*theta)}    #Relative risk
theta <- 0.01                            #Estimate of theta
pif.heatmap(X, theta = theta, rr = rr)

#Save file using ggplot2
#require(ggplot2)
#ggsave("My Potential Impact Fraction Heatmap Analysis.pdf")

#Change pif estimation method to kernel

pif.heatmap(X, theta = theta, rr = rr, method = "kernel")

#Example 2
#------------------------------------------------------------------
X     <- data.frame(Exposure = rbeta(100, 1, 0.2))
theta <- c(0.12, 1)
rr    <- function(X,theta){X*theta[1] + theta[2]}
cft   <- function(X, a, b){sin(a*X)*b}
pif.heatmap(X, theta = theta, rr = rr, cft = cft, 
     nmesh = 15, colors = rainbow(30), method = "empirical",
     title = "PIF with counterfactual cft(X) = sin(a*X)*b")

#Change estimation method to approximate
Xmean <- data.frame(mean(X[,"Exposure"]))
Xvar  <- var(X)
pif.heatmap(Xmean, Xvar = Xvar, theta = theta, rr = rr, cft = cft, 
     nmesh = 15, colors = rainbow(30), method = "approximate",
     title = "PIF with counterfactual cft(X) = sin(a*X)*b")


#Example 3: Plot univariate counterfactuals
#------------------------------------------------------------------
X       <- data.frame(rgamma(100, 1, 0.2))
theta   <- c(0.12, 1)
rr      <- function(X,theta){X*theta[1] + theta[2]}
cft     <- function(X, a, b){sqrt(a*X)}   #Leave two variables in it
pifplot <- pif.heatmap(X, theta = theta, rr = rr, 
 cft = cft, mina = 0, maxa = 1, minb = 0, maxb = 0, 
 legendtitle = "Potential Impact Fraction",
 title ="Univariate counterfactual", ylab = "", colors = topo.colors(10))
pifplot

#You can also add additional ggplot objects
#require(ggplot2)
#pifplot + annotate("text", x = 0.25, y = 0.4, label = "10yr scenario") + 
#geom_vline(aes(xintercept = 0.5), linetype = "dashed") +
#geom_segment(aes(x = 0.25, y = 0.38, xend = 0.5, yend = 0.3), 
#arrow = arrow(length = unit(0.25, "cm")))

#Example 4: Plot counterfactual with categorical risks
#------------------------------------------------------------------
set.seed(18427)
X        <- data.frame(Exposure = 
              sample(c("Normal","Overweight","Obese"), 100, 
                      replace = TRUE, prob = c(0.4, 0.1, 0.5)))
thetahat <- c(1, 1.7, 2)

#Categorical relative risk function
rr <- function(X, theta){

   #Create return vector with default risk of 1
   r_risk <- rep(1, nrow(X))
   
   #Assign categorical relative risk
   r_risk[which(X[,"Exposure"] == "Normal")]      <- thetahat[1]
   r_risk[which(X[,"Exposure"] == "Overweight")]  <- thetahat[2]
   r_risk[which(X[,"Exposure"] == "Obese")]       <- thetahat[3]
   
   return(r_risk)
}

#Counterfactual of reducing a certain percent of obesity and overweight cases
#to normality
cft <- function(X, per.over, per.obese){

   #Find the overweight and obese individuals
   which_obese <- which(X[,"Exposure"] == "Obese")
   which_over  <- which(X[,"Exposure"] == "Overweight")
   
   #Reduce per.over \% of overweight and per.obese \% of obese
   X[sample(which_obese, length(which_obese)*per.obese),
       "Exposure"] <- "Normal"
   X[sample(which_over,  length(which_over)*per.over),
       "Exposure"] <- "Normal"
   
   return(X)
}

pif.heatmap(X, thetahat = thetahat, rr = rr, cft = cft, mina = 0, minb = 0, maxa = 1, 
            maxb = 1, title = "PIF of excess-weight reduction",
            xlab = "\% Overweight cases", ylab = "\% Obese cases",
            check_exposure = FALSE, check_rr = FALSE)
}
}
\seealso{
See \code{\link{pif}} for Potential Impact Fraction estimation,
  \code{\link{pif.sensitivity}} for sensitivity analysis of the convergence
 process,  \code{\link{pif.plot}} for a plot of Potential Impact
  Fraction as a function of the relative risk's parameter \code{theta}.
}
\author{
Rodrigo Zepeda-Tello \email{rzepeda17@gmail.com}

Dalia Camacho-García-Formentí \email{daliaf172@gmail.com}
}
