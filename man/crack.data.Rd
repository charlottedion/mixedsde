
\name{crack.data}
\alias{crack.data}
\title{
Crack Growth
}
\description{
  68 measurement for 164 times
}
\usage{
Virkler
}
\format{
Sixty-eight replicate constant amplitude tests in aluminium alloy were carried out to investigate the fatigue crack propagation. In each of these tests, the number of cycles that leads to fixed crack lengths was observed. Against the natural assumption that something is observed at fixed times, here the time is the dependent variable and the crack length is the independent variable. Therefore, from the mathematical viewpoint the crack length will here be treated as time vector \code{t}.

The crack data comes as a list of length 2. The first element is a matrix of dimension 68 x 164, the second entry a vector of length 164. The matrix contains 68 observation series of observed times (in 100000 load cycles) up to a fixed crack length increment, stored in the second list entry. 

The original observations, i.e., the integer number of load cycles and the whole crack length, can be loaded by \code{data(Virkler)}, which contains the matrix of observed times, and \code{timesVirkler}, which contains the actual crack length beginning by 9. 

We want to thank Eric J. Tuegel for providing us the data that were collected by Prof. B. M. Hillberry, published in Virkler et al. (1979). 
}

\examples{
data(crack.data)

X <- crack.data[[1]]
times <- crack.data[[2]]

plot(times, X[1,], type='l', ylim=range(X), xlab="increment of crack length in mm", ylab="time in load cycles / 100000")
for (i in 2:nrow(X)){
 lines(times, X[i,], col=i)
}

estim_OU <- mixedsde.fit(times, X, model="OU", random=1, estim.method='paramBayes', nMCMC=2000) 
plot(estim_OU)

pred_OU <- pred(estim_OU)

}
\source{
  \emph{Eric J. Tuegel}
}
\references{
  \emph{Virkler, D. A., Hillberry, B. M. and Goel, P. K. (1979).
 The Statistical Nature of Fatigue Crack Propagation.
Journal of Engineering Materials and Technology 101, 148--153.}
}





\keyword{datasets}
