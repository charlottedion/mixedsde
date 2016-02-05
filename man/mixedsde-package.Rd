\name{mixedsde-package}
\alias{mixedsde-package}
\alias{mixedsde}
\docType{package}
\title{
Density estimation in mixed stochastic differential models
}
\description{
This package proposes 3 methods for density estimation in the special context of stochastic differential equation with linear random effects in the drift.
}
\details{
\tabular{ll}{
Package: \tab mixedsde\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2016-01-19\cr
License: \tab What license is it under?\cr
}
 An overview of how to use the package, including the most important functions 
}
\author{
Charlotte Dion, Simone Hermann, Adeline Samson

Maintainer: Charlotte Dion <charlotte.dion@imag.fr>
}
\references{
See Bidimensional random effect estimation in mixed stochastic differential model, C. Dion and V. Genon-Catalot, \emph{Stochastic Inference for Stochastic Processes 2015, Springer Netherlands} \bold{1--28}

Maximum likelihood estimation for stochastic differential equations with random effects, M. Delattre, V. Genon-Catalot and A. Samson, \emph{Scandinavian Journal of Statistics 2012}, Vol 40, \bold{322--343}

Bayesian Prediction of Crack Growth Based on a Hierarchical Diffusion Model.  Hermann, S., Ickstadt, K. and C. M{\"u}ller (2015)

}
\keyword{ package }
\seealso{
Optional links to other man pages, e.g. 
 \code{\link[mixedsde.fit]{misdesde.fit}} 
}
\examples{
model = 'OU'; M <- 50;  T <- 100 ; delta <- 0.01; N <- floor(T/delta); sigma <- 0.01 ; random <- c(1,2); density.phi <- 'normalnormal';
param<- c(0.4, 0.5, 1, 0.2);
simu <- mixedsde.sim(M=M,T=T,N=N,model=model,random=random,density.phi=density.phi,param=param,sigma=sigma, invariant = 1)
X <- simu$X ; phi <- simu$phi; times <- simu$times
estim.method<- 'nonparam'
estim <- mixedsde.fit(times=times, X=X, model=model, random=random, estim.method= 'nonparam')

# maybe the examples from mixedsde.fit here ?!?
}
