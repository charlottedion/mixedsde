% Generated by roxygen2 (4.1.0): do not edit by hand
% Please edit documentation in R/mixedsde.fit.R
\docType{methods}
\name{plot,Bayes.fit,ANY-method}
\alias{plot,Bayes.fit,ANY-method}
\title{Plot method for the Bayesian estimation class object}
\usage{
\S4method{plot}{Bayes.fit,ANY}(x, plot.priorMean = FALSE, reduced = FALSE,
  style = c("chains", "acf", "density", "cred.int"), level = 0.05, true.phi,
  newwindow = FALSE, ...)
}
\arguments{
\item{x}{Bayes.fit class}

\item{plot.priorMean}{logical(1), if TRUE, prior means are added to the plots}

\item{reduced}{logical(1), if TRUE, the chains are reduced with the burn-in and thin rate}

\item{style}{one out of 'chains', 'acf', 'density' or 'cred.int'}

\item{level}{alpha for the credibility intervals, only for style 'cred.int', default = 0.05}

\item{true.phi}{only for style 'cred.int', for the case of known true values, e.g. for simulation}

\item{newwindow}{logical(1), if TRUE, a new window is opened for the plot}

\item{...}{optional plot parameters}
}
\description{
Plot method for the S4 class Bayes.fit
}
\examples{
random <- c(1,2); sigma <- 0.1; param <- c(3, 0.5, 5, 0.2)

sim <- mixedsde.sim(M = 20, T = 1, N = 50, model = 'OU', random = random,
       density.phi = 'normalnormal', param = param, sigma = sigma, X0 = 0, op.plot = 1)

# here: only 100 iterations for example - should be much more!
prior <- list(m = param[c(1,3)], v = param[c(1,3)], alpha.omega = c(11,11),
           beta.omega = param[c(2,4)]^2*10, alpha.sigma = 10, beta.sigma = sigma^2*9)
estim_Bayes <- mixedsde.fit(times = sim$times, X = sim$X, model = 'OU', random = random,
                estim.method = 'paramBayes', prior = prior, nMCMC = 100)
plot(estim_Bayes)
plot(estim_Bayes, style = 'cred.int', true.phi = sim$phi)
plot(estim_Bayes, style = 'acf')
plot(estim_Bayes, style = 'density')
}
\references{
Dion, C., Hermann, S. and Samson, A. (2016). Mixedsde: a R package to fit mixed stochastic differential equations.
}

