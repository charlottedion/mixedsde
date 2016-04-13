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

Bayesian Prediction of Crack Growth Based on a Hierarchical Diffusion Model.  Hermann, S., Ickstadt, K. and C. Muller (2016)

}
\keyword{ package }

\examples{
# Frequentist estimation, two random effects

model = 'CIR'; M <- 200;  T <- 10 ; delta <- 0.001; N <- floor(T/delta); sigma <- 0.01
random <- c(1,2); density.phi <- 'gammainvgamma2'
param<- c(1.8, 0.8, 8, 0.05);
simu <- mixedsde.sim(M=M,T=T,N=N,model=model,random=random,density.phi=density.phi,param=param,
              sigma=sigma, invariant = 1)
X <- simu$X ; phi <- simu$phi; times <- simu$times
estim.method<- 'nonparam'
estim <- mixedsde.fit(times=times, X=X, model=model, random=random, estim.method= 'nonparam')
outputsNP <- out(estim)

summary(estim)
print(estim)
\dontrun{
plot(estim)
validation <- valid(estim, numj=floor(runif(1,1,M)))
}
estim.method<-'paramML'
estim_param <- mixedsde.fit(times= times, X= X, model= model, random= random, 
	estim.method = 'paramML')
outputsP <- out(estim_param)

summary(estim_param)

\dontrun{
plot(estim_param)
test1 <- pred(estim, invariant = 1)
test2 <- pred(estim_param, invariant = 1)
}
cutoff <- outputsNP$cutoff
phihat <- outputsNP$estimphi
phihat_trunc <- outputsNP$estimphi_trunc
par(mfrow=c(1,2))
plot.ts(phi[1,], phihat[1,], xlim=c(0, 15), ylim=c(0,15), pch=18); abline(0,1)
points(phi[1,]*(1-cutoff), phihat[1,]*(1-cutoff), xlim=c(0, 20), ylim=c(0,20), pch=18, col='red')
abline(0,1)
plot.ts(phi[2,], phihat[2,], xlim=c(0, 15), ylim=c(0,15),pch=18); abline(0,1)
points(phi[2,]*(1-cutoff), phihat[2,]*(1-cutoff), xlim=c(0, 20), ylim=c(0,20), pch=18, col='red')
abline(0,1)



# Parametric Bayesian estimation one random effect

model <- 'OU'; random <- 1; sigma <- 0.1; fixed <- 5
M <- 50 ; T <- 1; N <- 100
density.phi <- 'normal'; param <- c(3, 0.5)

simu <- mixedsde.sim(M, T = T, N = N, model= model, random = random, fixed = fixed, 
      density.phi= density.phi, param= param, sigma= sigma, X0 = 0)
X <- simu$X; phi <- simu$phi; times <- simu$times
#plot(times, X[1,], ylim = range(X), type = 'l'); for(i in 2:M) lines(times, X[i,])

estim_Bayes_withoutprior <- mixedsde.fit(times, X= X, model = model, random = random, 
            estim.method = 'paramBayes', nMCMC = 100)  # nMCMC should be much larger
plot(estim_Bayes_withoutprior)

prior <- list(m = c(param[1], fixed), v = c(param[1], fixed), alpha.omega = 11, 
                beta.omega = param[2]^2*10, alpha.sigma = 10, beta.sigma = sigma^2*9)
estim_Bayes <- mixedsde.fit(times, X = X, model = model, random = random, 
              estim.method = 'paramBayes', prior = prior, nMCMC = 100)

plot(estim_Bayes)
outputBayes <- out(estim_Bayes)
summary(outputBayes)
(results_Bayes <- summary(estim_Bayes))
plot(estim_Bayes, style = 'cred.int', true.phi = phi)
plot(estim_Bayes_withoutprior, style = 'cred.int', true.phi = phi, reduced = TRUE)

plot2compare(estim_Bayes, estim_Bayes_withoutprior, names = c('with prior', 'without prior'))

print(estim_Bayes)

\dontrun{
pred.result <- pred(estim_Bayes)
summary(out(pred.result))
plot(pred.result)

pred.result.trajectories <- pred(estim_Bayes, trajectories = TRUE)

validbayes <- valid(estim_Bayes, numj = 1)
}
}
