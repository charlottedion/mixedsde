#' Computation Of The Log Likelihood In Mixed Stochastic Differential Equations
#' 
#' @description Computation of -2 loglikelihood of the mixed SDE with Normal distribution of the random effects
#'  \eqn{dXj(t)= (\alpha_j- \beta_j Xj(t))dt + \sigma a(Xj(t)) dWj(t)}.
#' @param mu current value of the mean of the normal distribution
#' @param omega current value of the standard deviation of the normal distribution
#' @param estimphi vector or matrix of estimators of the random effects 
#' @param V vector of the M sufficient statistics V (see \code{\link{UV}})
#' @param random random effects in the drift: 1 if one additive random effect, 2 if one multiplicative random effect or c(1,2) if 2 random effects.
#' @return
#' \item{L}{value of -2 x loglikelihood}
#' @references
#' Maximum likelihood estimation for stochastic differential equations with random effects, M. Delattre, V. Genon-Catalot and A. Samson, \emph{Scandinavian Journal of Statistics 2012}, Vol 40, \bold{322--343}
#' 
#' 
#' 
#' 
# 
# NE PAS EFFACER fonction pour omega diagonale
likelihoodNormal <- function(mu, omega, estimphi, V, random) { 
  if (length(random) == 1) { 
    Omega <- omega^2 
    L <- sum(log(1 + Omega * V)) + sum(V/(1 + Omega * V) * (mu- estimphi)^2)
    } 
  
  if (length(random) == 2) { 
    Omega <- matrix(c(omega[1]^2, 0, 0, omega[2]^2), 2, 2, byrow = TRUE) 
    M <- dim(estimphi)[2] 
    loglik <- vector(length = M)
    I2 <- diag(c(1, 1)) 
    for (j in 1:M) { 
        A <- (I2 + V[[j]] %*% Omega) 
        Rinv <- solve(A) %*% V[[j]] 
        b <- mu - estimphi[, j] 
        loglik[j] <- log(det(A)) + t(b) %*% Rinv %*% b 
      } 
    L <- sum(loglik) 
    } 
  return(L) 
} 


# # NE PAS EFFACER fonction pour omega NON diagonale
# likelihoodNormal = function(mu, omega, cova, estimphi, V, random) {
#     
#     
#     if (length(random) == 1) {
#       
#       Omega <- omega^2
#       L <- sum(log(1 + Omega * V)) + sum(V/(1 + Omega * V) * (mu - estimphi)^2)
#     }
#     
#     if (length(random) == 2) {
#         
#         Omega <- matrix(c(omega[1]^2, cova, cova, omega[2]^2), 2, 2, byrow = TRUE)
#        
#         M <- dim(estimphi)[2]
#         loglik <- vector(length = M)
#         I2 <- diag(c(1, 1))
#         
#         for (j in 1:M) {
#             A <- (I2 + V[[j]] %*% Omega)
#             #print(det(A))
#             #print(Omega)
#             Rinv <- solve(A) %*% V[[j]]
#             b <- mu - estimphi[ , j]
#             loglik[j] <- log(det(A)) + t(b) %*% Rinv %*% b
#         }
#         L <- sum(loglik)
#         
#     }
#     
#     return(L)
#     
# }
# 



