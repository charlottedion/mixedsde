#' Maximization Of The Log Likelihood In Mixed Stochastic Differential Equations
#' 
#' @description Maximization of the loglikelihood of the mixed SDE with Normal distribution of the random effects
#'  \eqn{dXj(t)= (\alpha_j- \beta_j Xj(t))dt + \sigma a(Xj(t)) dWj(t)}, done with \code{\link[=mixedsde]{likelihoodNormal}}
#' @param estimphi vector or matrix of estimators of the random effects
#' @param V vector of the M sufficient statistics V
#' @param random random effects in the drift: 1 if one additive random effect, 2 if one multiplicative random effect or c(1,2) if 2 random effects.
#' @return
#' \item{mu}{estimated value of the mean}
#' \item{Omega}{estimated value of the variance}
#' \item{cova}{estimated value of the covariance}
#' \item{BIChere}{BIC indicator}
#' \item{AIChere}{AIC indicator}

# NE PAS EFFACER fonction pour omega diagonale
EstParamNormal = function(estimphi, V, random) {
    
    if (length(random) == 1) {
     
        ln = function(param) {
            likelihoodNormal(param[1], param[2], estimphi, V, random)
        }
        init.mu <- mean(estimphi)
        init.omega <- sd(estimphi)
        
        res = optim(c(init.mu, init.omega), f = ln, method = "Nelder-Mead")
        mu = res$par[1]
        omega = abs(res$par[2])
        
        BIChere <- likelihoodNormal(mu, omega,  estimphi, V, random) + log(2 * pi) + log(length(estimphi)) * 2
        AIChere <- likelihoodNormal(mu, omega, estimphi, V, random) + log(2 * pi) + 2
    }
    
    if (length(random) == 2) {
        
        loglik_Omegadiag <- function(param, estimphi, V) { 
          mu <- param[1:2] 
          omega <- param[3:4]
        
          return(likelihoodNormal(mu, omega, estimphi, V, random)) 
          }
        
        init.mu <- c(mean(estimphi[1, ]), mean(estimphi[2, ]))
        init.omega <- c(sd(estimphi[1, ]), sd(estimphi[2, ]))
        
        res = optim(c(init.mu, init.omega), loglik_Omegadiag, gr = NULL, estimphi, V, method = "Nelder-Mead")
        mu = c(res$par[1], res$par[2])
        omega = abs(c(res$par[3], res$par[4]))

        BIChere <- loglik_Omegadiag(c(mu, omega), estimphi, V) + log(4 * pi) + log(dim(estimphi)[2]) * 4 # CHANGER 
        AIChere <- loglik_Omegadiag(c(mu, omega), estimphi, V) + log(4 * pi) + 4
        
    }
    return(list(mu = mu, omega = omega, BIChere = BIChere, AIChere = AIChere))
} 



# #NE PAS EFFACER fonction pour omega NON diagonale
# EstParamNormal = function(estimphi, V, random) {
#   
#   if (length(random) == 1) {
#     cova <- 0
#     ln = function(param) {
#       likelihoodNormal(param[1], param[2], cova = 0, estimphi, V, random)
#     }
#     init.mu <- mean(estimphi)
#     init.omega <- sd(estimphi)
#     
#     res = optim(c(init.mu, init.omega), f = ln, method = "Nelder-Mead")
#     mu = res$par[1]
#     omega = abs(res$par[2])
#     
#     BIChere <- likelihoodNormal(mu, omega, cova = 0, estimphi, V, random) + log(2 * pi) + log(length(estimphi)) * 2
#     AIChere <- likelihoodNormal(mu, omega, cova = 0, estimphi, V, random) + log(2 * pi) + 2
#   }
#   
#   if (length(random) == 2) {
#     
#     loglik_Omegacov = function(param, estimphi, V) {
#       mu <- param[1:2]
#       omega <- param[3:4]
#       cova <- param[5]
#       return(likelihoodNormal(mu, omega, cova, estimphi, V, random))
#     }
#     
#     init.mu <- c(mean(estimphi[1, ]), mean(estimphi[2, ]))
#     init.omega <- c(sd(estimphi[1, ]), sd(estimphi[2, ]))
#     
#     init.cov <- cov(estimphi[1, ], estimphi[2, ])
#     
#     res <- optim(c(init.mu, init.omega, init.cov), loglik_Omegacov, gr = NULL, estimphi, V, method = "Nelder-Mead")
#     
#     mu <- c(res$par[1], res$par[2])
#     omega <- abs(c(res$par[3], res$par[4]))
#     cova <- res$par[5]
#     
#     BIChere <- loglik_Omegacov(c(mu, omega, cova), estimphi, V) + log(4 * pi) + log(dim(estimphi)[2]) * 4
#     AIChere <- loglik_Omegacov(c(mu, omega, cova), estimphi, V) + log(4 * pi) + 4
#     
#   }
#   return(list(mu = mu, omega = omega, cova = cova, BIChere = BIChere, AIChere = AIChere))
# } 
