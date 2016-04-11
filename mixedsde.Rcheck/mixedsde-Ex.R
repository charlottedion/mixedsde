pkgname <- "mixedsde"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
options(pager = "console")
library('mixedsde')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
cleanEx()
nameEx("mixedsde-package")
### * mixedsde-package

flush(stderr()); flush(stdout())

### Name: mixedsde-package
### Title: Density estimation in mixed stochastic differential models
### Aliases: mixedsde-package mixedsde
### Keywords: package

### ** Examples

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
plot(estim)
summary(estim)
print(estim)

validation <- valid(estim, numj=floor(runif(1,1,M)))

estim.method<-'paramML'
estim_param <- mixedsde.fit(times= times, X= X, model= model, random= random, 
	estim.method = 'paramML')
outputsP <- out(estim_param)
plot(estim_param)
summary(estim_param)

test1 <- pred(estim, invariant = 1)
test2 <- pred(estim_param, invariant = 1)

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
plot(times, X[1,], ylim = range(X), type = 'l'); for(i in 2:M) lines(times, X[i,])

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

pred.result <- pred(estim_Bayes)
summary(out(pred.result))
plot(pred.result)

pred.result.trajectories <- pred(estim_Bayes, trajectories = TRUE)

validbayes <- valid(estim_Bayes, numj = 1)



graphics::par(get("par.postscript", pos = 'CheckExEnv'))
cleanEx()
nameEx("mixedsde.fit")
### * mixedsde.fit

flush(stderr()); flush(stdout())

### Name: mixedsde.fit
### Title: Estimation Of The Random Effects In Mixed Stochastic
###   Differential Equations
### Aliases: mixedsde.fit
### Keywords: estimation

### ** Examples

# Frequentist estimation
# Two random effects
model = 'CIR'; M <- 200;  T <- 10 ; delta <- 0.001; N <- floor(T/delta); sigma <- 0.01 ;
random <- c(1,2); density.phi <- 'gammainvgamma2'; param<- c(1.8, 0.8, 8, 0.05);
simu <- mixedsde.sim(M=M, T=T, N=N, model=model,random=random, density.phi=density.phi,
               param=param, sigma=sigma, invariant = 1)
X <- simu$X ; phi <- simu$phi; times <- simu$times
estim.method<- 'nonparam'
estim <- mixedsde.fit(times=times, X=X, model=model, random=random, estim.method= 'nonparam')
#To stock the results of the function, use method \code{out}
#which put the outputs of the function on a list according to the frequentist or
# Bayesian approach.
outputsNP <- out(estim)
plot(estim)
# It represents the bidimensional density, the histogram of the first estimated random
# effect \eqn{\alpha} with the  marginal of \eqn{\hat{f}} from the first coordonate which
# estimates  the density of \eqn{\alpha}. And the same for the second random effect
# \eqn{\beta}. More, it plots a qq-plot for the sample of estimator of the random effects
# compared with the quantiles of a Gaussian sample with the same mean and standard deviation.

summary(estim)
print(estim)
# Validation
# If numj is fixed by the user: this function simulates Mrep =100 (by default) new
# trajectories with the value of the estimated random effect. Then it plots on the
# left graph the Mrep new trajectories \eqn{(Xnumj^{k}(t1), ... Xnumj^{k}(tN)),
# k= 1, ... Mrep} with in red the true trajectory \eqn{(Xnumj(t1), ... Xnumj(tN))}.
#The right graph is a qq-plot of the quantiles of samples
# \eqn{(Xnumj^{1}(ti), ... Xnumj^{Mrep}(ti))}
# for each time \eqn{ti} compared with the uniform quantiles. The outputs of the function
# are: a matrix \code{Xnew} dimension Mrepx N+1, vector of quantiles \code{quantiles} length
# N and the number of the trajectory for the plot \code{plotnumj= numj}
# If numj is not precised by the user, then, this function simulates Mrep =100 (by default)
# new trajectories for each estimated random effect. Then left graph is a plot of the Mrep
# new trajectories \eqn{(Xj^{k}(t1), ... Xj^{k}(tN)), k= 1, ... Mrep}
#for a randomly chosen number j with in red the true trajectory \eqn{(Xj(t1), ... Xj(tN))}.
#The right graph is a qq-plot of the quantiles of samples \eqn{(Xj^{1}(ti), ... Xj^{Mrep}(ti))},
# for the same j and for each time \eqn{ti}. The outputs of the function are: a list of
# matrices \code{Xnew} length M, matrix of quantiles \code{quantiles} dimension MxN
# and the number of the trajectory for the plot \code{plotnumj}

validation <- valid(estim,  numj=floor(runif(1,1,M)))

# Parametric estimation
estim.method<-'paramML'
estim_param <- mixedsde.fit(times= times, X= X, model= model, random= random,
           estim.method = 'paramML')
outputsP <- out(estim_param)
plot(estim_param)
summary(estim_param)

# Prediction for the frequentist approach
# This function uses the estimation of the density function to simulate a
# new sample of random effects according to this density. If \code{plot.pred =1} (default)
# is plots on the top the predictive random effects versus the estimated random effects
# from the data. On the bottom, the left graph is the true trajectories, on the right
#the predictive trajectories and the empiric prediciton intervals at level
# \code{level=0.05} (defaut). The function return on a list the prediction of phi
# \code{phipred}, the prediction of X \code{Xpred}, and the indexes of the
# corresponding true trajectories \code{indexpred}

test1 <- pred(estim,  invariant  = 1)
test2 <- pred(estim_param, invariant  = 1)

# More graph
fhat <- outputsNP$estimf
fhat_trunc <- outputsNP$estimf.trunc
fhat_param <- outputsP$estimf

gridf <- outputsNP$gridf; gridf1 <- gridf[1,]; gridf2 <- gridf[2,]

marg1 <- ((max(gridf2)-min(gridf2))/length(gridf2))*apply(fhat,1,sum)
marg1_trunc <- ((max(gridf2)-min(gridf2))/length(gridf2))*apply(fhat_trunc,1,sum)
marg2 <- ((max(gridf1)-min(gridf1))/length(gridf1))*apply(fhat,2,sum)
marg2_trunc <- ((max(gridf1)-min(gridf1))/length(gridf1))*apply(fhat_trunc,2,sum)

marg1_param <- ((max(gridf2)-min(gridf2))/length(gridf2))*apply(fhat_param,1,sum)
marg2_param <- ((max(gridf1)-min(gridf1))/length(gridf1))*apply(fhat_param,2,sum)
f1 <-  (gridf1^(param[1]-1))*exp(-gridf1/param[2])/((param[2])^param[1]*gamma(param[1]))
f2 <-  (gridf2^(-param[3]-1)) * exp(-(1/param[4])*(1/gridf2)) *
 ((1/param[4])^param[3])*(1/gamma(param[3]))
par(mfrow=c(1,2))
plot(gridf1,f1,type='l', lwd=1,  xlab='', ylab='')
lines(gridf1,marg1_trunc,col='blue', lwd=2)
lines(gridf1,marg1,col='blue', lwd=2, lty=2)
lines(gridf1,marg1_param,col='red', lwd=2, lty=2)
plot(gridf2,f2,type='l', lwd=1, xlab='', ylab='')
lines(gridf2,marg2_trunc,col='blue', lwd=2)
lines(gridf2,marg2,col='blue', lwd=2, lty=2)
lines(gridf2,marg2_param,col='red', lwd=2, lty=2)

cutoff <- outputsNP$cutoff
phihat <- outputsNP$estimphi
phihat_trunc <- outputsNP$estimphi.trunc
par(mfrow=c(1,2))
plot.ts(phi[1,], phihat[1,], xlim=c(0, 15), ylim=c(0,15), pch=18); abline(0,1)
points(phi[1,]*(1-cutoff), phihat[1,]*(1-cutoff), xlim=c(0, 20), ylim=c(0,20),pch=18, col='red');
abline(0,1)
plot.ts(phi[2,], phihat[2,], xlim=c(0, 15), ylim=c(0,15),pch=18); abline(0,1)
points(phi[2,]*(1-cutoff), phihat[2,]*(1-cutoff), xlim=c(0, 20), ylim=c(0,20),pch=18, col='red');
abline(0,1)

# one random effect:

model <-'OU'
random <- 1
M <- 80; T <- 100  ; N <- 2000
sigma <- 0.1 ; X0 <- 0
density.phi <- 'normal'
fixed <- 2 ; param <- c(1, 0.2)
#-------------------
#- simulation
simu <- mixedsde.sim(M,T=T,N=N,model=model,random=random, fixed=fixed,density.phi=density.phi,
               param=param, sigma=sigma, X0=X0)
X <- simu$X
phi <- simu$phi
times <-simu$times
plot(times, X[10,], type='l')

#- parametric estimation
estim.method<-'paramML'
estim_param <- mixedsde.fit(times, X=X, model=model, random=random, estim.fix= 1,
               estim.method=estim.method)
outputsP <- out(estim_param)
estim.fixed <- outputsP$estim.fixed #to compare with fixed
#or
estim_param2 <- mixedsde.fit(times, X=X, model=model, random=random, fixed = fixed,
             estim.method=estim.method)
outputsP2 <- out(estim_param2)
#- nonparametric estimation
estim.method <- 'nonparam'
estim <- mixedsde.fit(times, X, model=model, random=random, fixed = fixed,
           estim.method=estim.method)
outputsNP <- out(estim)

plot(estim)
print(estim)
summary(estim)

plot(estim_param)
print(estim_param)
summary(estim_param)

valid1 <- valid(estim,  numj=floor(runif(1,1,M)))
test1 <- pred(estim )
test2 <- pred(estim_param)


# Parametric Bayesian estimation
# one random effect
random <- 1; sigma <- 0.1; fixed <- 5; param <- c(3, 0.5)
sim <- mixedsde.sim(M = 50, T = 1, N = 100, model = 'OU', random = random, fixed = fixed,
       density.phi = 'normal',param= param, sigma= sigma, X0 = 0, op.plot = 1)

# here: only 100 iterations for example - should be much more!
estim_Bayes_withoutprior <- mixedsde.fit(times = sim$times, X = sim$X, model = 'OU',
             random, estim.method = 'paramBayes',  nMCMC = 100)
prior <- list( m = c(param[1], fixed), v = c(param[1], fixed), alpha.omega = 11,
            beta.omega = param[2]^2*10, alpha.sigma = 10, beta.sigma = sigma^2*9)
estim_Bayes <- mixedsde.fit(times = sim$times, X = sim$X, model = 'OU', random,
           estim.method = 'paramBayes', prior = prior, nMCMC = 100)

validation <- valid(estim_Bayes, numj = 10)
plot(estim_Bayes)
outputBayes <- out(estim_Bayes)
summary(outputBayes)
(results_Bayes <- summary(estim_Bayes))
plot(estim_Bayes, style = 'cred.int', true.phi = sim$phi)
plot(estim_Bayes_withoutprior, style = 'cred.int', true.phi = sim$phi, reduced = TRUE)

plot2compare(estim_Bayes, estim_Bayes_withoutprior, names = c('with prior', 'without prior'))

print(estim_Bayes)

pred.result <- pred(estim_Bayes)
summary(out(pred.result))
plot(pred.result)

pred.result.trajectories <- pred(estim_Bayes, trajectories = TRUE)

# second example
## Not run: 
##D random <- 2; sigma <- 0.2; fixed <- 5; param <- c(3, 0.5)
##D sim <- mixedsde.sim(M = 20, T = 1, N = 100, model = 'CIR', random = random,
##D         fixed = fixed, density.phi = 'normal',param = param, sigma = sigma, X0 = 0.1, op.plot = 1)
##D 
##D prior <- list(m = c(fixed, param[1]), v = c(fixed, param[1]), alpha.omega = 11,
##D          beta.omega = param[2]^2*10, alpha.sigma = 10, beta.sigma = sigma^2*9)
##D 
##D estim_Bayes <- mixedsde.fit(times = sim$times, X = sim$X, model = 'CIR', random = random,
##D                  estim.method = 'paramBayes', prior = prior, nMCMC = 1000)
##D plot(estim_Bayes)
##D outputBayes <- out(estim_Bayes)
##D summary(outputBayes)
##D (results_Bayes <- summary(estim_Bayes))
##D plot(estim_Bayes, style = 'cred.int', true.phi = sim$phi, reduced = TRUE)
##D 
##D print(estim_Bayes)
##D pred.result <- pred(estim_Bayes)
##D summary(out(pred.result))
##D plot(pred.result)
## End(Not run)

# for two random effects
random <- c(1,2); sigma <- 0.1; param <- c(3, 0.5, 5, 0.2)

sim <- mixedsde.sim(M = 20, T = 1, N = 100, model = 'OU', random = random,
       density.phi = 'normalnormal', param = param, sigma = sigma, X0 = 0, op.plot = 1)

# here: only 200 iterations for example - should be much more!
estim_Bayes_withoutprior <- mixedsde.fit(times = sim$times, X = sim$X, model = 'OU',
             random = random, estim.method = 'paramBayes', nMCMC = 100)
plot(estim_Bayes_withoutprior, style = 'cred.int', true.phi = sim$phi, reduced = TRUE)

prior <- list(m = param[c(1,3)], v = param[c(1,3)], alpha.omega = c(11,11),
           beta.omega = param[c(2,4)]^2*10, alpha.sigma = 10, beta.sigma = sigma^2*9)
estim_Bayes <- mixedsde.fit(times = sim$times, X = sim$X, model = 'OU', random = random,
                estim.method = 'paramBayes', prior = prior, nMCMC = 100)
outputBayes <- out(estim_Bayes)
summary(outputBayes)
summary(estim_Bayes)
plot(estim_Bayes)
plot(estim_Bayes, style = 'cred.int', true.phi = sim$phi)
print(estim_Bayes)

pred.result <- pred(estim_Bayes)


# invariant case

random <- 1; sigma <- 0.1; fixed <- 5; param <- c(3, 0.5)
sim <- mixedsde.sim(M = 50, T = 5, N = 100, model = 'OU', random = random, fixed = fixed,
           density.phi = 'normal',param = param, sigma = sigma, invariant = 1, op.plot = 1)

prior <- list(m = c(param[1], fixed), v = c(param[1], 1e-05), alpha.omega = 11,
       beta.omega = param[2]^2*10, alpha.sigma = 10, beta.sigma = sigma^2*9)
estim_Bayes <- mixedsde.fit(times = sim$times, X = sim$X, model = 'OU', random,
       estim.method = 'paramBayes', prior = prior, nMCMC = 100)
plot(estim_Bayes)

pred.result <- pred(estim_Bayes, invariant = 1)
pred.result.traj <- pred(estim_Bayes, invariant = 1, trajectories = TRUE)



graphics::par(get("par.postscript", pos = 'CheckExEnv'))
cleanEx()
nameEx("mixedsde.sim")
### * mixedsde.sim

flush(stderr()); flush(stdout())

### Name: mixedsde.sim
### Title: Simulation Of A Mixed Stochastic Differential Equation
### Aliases: mixedsde.sim

### ** Examples

#Simulation of 5 trajectories of the OU SDE with random =1, and a Gamma distribution.

simuOU <- mixedsde.sim(M=5, T=10,N=1000,model='OU', random=1,fixed=0.5,
density.phi='gamma', param=c(1.8, 0.8) , sigma=0.1,op.plot=1)
X <- simuOU$X ;
phi <- simuOU$phi
hist(phi)



cleanEx()
nameEx("mixture.sim")
### * mixture.sim

flush(stderr()); flush(stdout())

### Name: mixture.sim
### Title: Simulation Of A Mixture Of Two Normal Or Gamma Distributions
### Aliases: mixture.sim

### ** Examples

density.phi <- 'mixture.gamma'
param <- c(0.2,1.8,0.5,5.05,1); M <- 200
gridf <- seq(0, 8, length = 200)
f <- param[1] * 1/gamma(param[2]) * (gridf)^(param[2]-1) *
           exp(-(gridf) / param[3]) / param[3]^param[2] +
	(1-param[1]) * 1/gamma(param[4]) * (gridf)^(param[4]-1) *
	    exp(-(gridf) / param[5]) / param[5]^param[4]
Y <- mixture.sim(M, density.phi, param)
hist(Y)
lines(gridf, f)



cleanEx()
nameEx("neuronal.data")
### * neuronal.data

flush(stderr()); flush(stdout())

### Name: neuronal.data
### Title: Trajectories Interspike Of A Single Neuron Of A Ginea Pig
### Aliases: neuronal.data
### Keywords: data

### ** Examples

require(plot3D)
model <- "OU"
random <- c(1,2)
M <- 240     # number of trajectories, number of rows of the matrix of the data
T <- 0.3     # width of the interval of observation 
delta <- 0.00015   # step time
N <- T/delta  # number of points in the time interval 2000
# load ("data/neuronal.data.rda")
data(neuronal.data)
X <- neuronal.data[[1]]
times <- neuronal.data[[2]]

#plot(times,X[10, ], type = 'l', xlab = 'time', ylab='', col = 'blue', ylim=c(0,0.016))

random <- c(1,2)

#- nonparametric estimation
estim.method <- 'nonparam'
estim <- mixedsde.fit(times=times, X=X, model=model, random=random,  estim.method='nonparam') 

#- parametric estimation   
estim.method<-'paramML'
estim_param <- mixedsde.fit(times=times, X=X, model=model, random= random, estim.method= 'paramML')

#- implemented methods
# plot(estim); 
print(estim); #valid(estim)
print(estim_param); #plot(estim_param);  valid(estim_param)

#test1 <- pred(estim, X,  estim.method= 'nonparam',times = times)
#test2 <- pred(estim_param, X,estim.method= 'paramML', times = times) 

#- Other possible plots
par(mfrow=c(1,2))

outputsNP <-  out(estim)
outputsP <- out(estim_param)
fhat <- outputsNP$estimf
fhat_param <- outputsP$estimf 

 gridf <- outputsNP$gridf
 gridf1 <- gridf[1,]; gridf2 <- gridf[2,]
 marg1 <- ((max(gridf2)-min(gridf2))/length(gridf2))*apply(fhat,1,sum) #with cutoff
 marg2 <- ((max(gridf1)-min(gridf1))/length(gridf1))*apply(fhat,2,sum)
 marg1_param <- ((max(gridf2)-min(gridf2))/length(gridf2))*apply(fhat_param,1,sum) 
 marg2_param <- ((max(gridf1)-min(gridf1))/length(gridf1))*apply(fhat_param,2,sum)

 plot(gridf1,marg1,type='l', col='red')
 lines(gridf1,marg1_param, lwd=2, col='red')
 plot(gridf2, marg2,type='l', col='red')
 lines(gridf2,marg2_param, lwd=2, col='red')


# Bayesian
ind <- seq(1, 2000, by = 10)
estim_Bayes <- mixedsde.fit(times[ind], X[,ind], model = "OU", random = 1, 
              estim.method = "paramBayes", nMCMC = 1000) 
plot(estim_Bayes)
pred_Bayes1 <- pred(estim_Bayes)
pred_Bayes2 <- pred(estim_Bayes, trajectories = TRUE)



graphics::par(get("par.postscript", pos = 'CheckExEnv'))
### * <FOOTER>
###
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
