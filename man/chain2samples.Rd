% Generated by roxygen2 (4.1.0): do not edit by hand
% Please edit documentation in R/BayesianNormal.R
\name{chain2samples}
\alias{chain2samples}
\title{Removing Of Burn-in Phase And Thinning}
\usage{
chain2samples(res, burnIn, thinning)
}
\arguments{
\item{res}{Bayes.fit class object}

\item{burnIn}{number of burn-in samples}

\item{thinning}{thinning rate}
}
\description{
Transfers class object Bayes.fit from the original to the thinned chains
}

