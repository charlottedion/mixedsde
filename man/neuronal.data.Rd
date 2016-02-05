
\name{neuronal.data}
\alias{neuronal.data}
\title{
  Trajectories Interspike Of A Single Neuron Of A Ginea Pig
}
\description{
  The \code{neuronal.data} data has 240 measurements of the membrane potential in volts for one single neuron of a pig between the spikes, along time, with 2000 points for each. The step time is \eqn{delta= 0.00015} s. 
}
\usage{
neuronal.data
}
\format{
This data frame has a list form of length 2. The first element in the matrix named \code{Xreal}. Each row is a trajectory, that one can model by a diffusion process with random effect. The realisation can be assumed independent. The second element is a vector of times of observations \code{times}
}

\examples{
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

plot(times,X[10, ], type = 'l', xlab = 'time', ylab='', col = 'blue', ylim=c(0,0.016))

random <- c(1,2)

#- nonparametric estimation
estim.method <- 'nonparam'
estim <- mixedsde.fit(times=times, X=X, model=model, random=random,  estim.method='nonparam') 

#- parametric estimation   
estim.method<-'paramML'
estim_param <- mixedsde.fit(times=times, X=X, model=model, random= random, estim.method= 'paramML')

#- implemented methods
# plot(estim); print(estim); valid(estim)
print(estim_param); plot(estim_param); # Freq.valid(estim_param)

# validation <- Freq.valid(estim, X=X, times = times)
# it produce the same method with estim_param

test1 <- pred(estim, X,  estim.method= 'nonparam',times = times)
test2 <- pred(estim_param, X,estim.method= 'paramML', times = times) 

#- Other possible plots
par(mfrow=c(1,2))

outputsNP <-  out(estim)
outputsP <- out(estim_param)
fhat <- outputsNP$estimf
fhat_param <- outputsP$estimf 

gridf <- outputsNP$gridf
gridf1 <- gridf[1,]; gridf2 <- gridf[2,]
# marg1T03OU <- ((max(gridf2)-min(gridf2))/length(gridf2))*apply(fhatT03OU,1,sum) #with cutoff
# marg2T03OU <- ((max(gridf1)-min(gridf1))/length(gridf1))*apply(fhatT03OU,2,sum)
# marg1_paramT03OU <- ((max(gridf2)-min(gridf2))/length(gridf2))*apply(fhat_paramT03OU,1,sum) 
# marg2_paramT03OU <- ((max(gridf1)-min(gridf1))/length(gridf1))*apply(fhat_paramT03OU,2,sum)

# plot(gridf1,marg1T03OU,type='l', col='red')
# lines(gridf1,marg1_paramT03OU, lwd=2, col='red')
# plot(gridf2, marg2T03OU,type='l', col='red')
# lines(gridf2,marg2_paramT03OU, lwd=2, col='red')


# Bayesian
ind <- seq(1, 2000, by=10)
estim_Bayes <- mixedsde.fit(times[ind], X[,ind], model="OU", random=1, estim.method='paramBayes', nMCMC=2000) 
pred_Bayes <- pred(estim_Bayes, trajectories=TRUE)
}

\source{
  Rune Berg, Lansky, Yu
  \emph{?}
}
\references{
 Rune Berg, Lansky, Yu
  \emph{?}
}





\keyword{datasets}
