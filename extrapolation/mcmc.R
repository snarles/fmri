####
##  Use MCMC to fit the mixture of binomials
####

library(rstan)
mix_model   <- "
    data{
        int K;
        int N;
        int Y[N];
    }
    parameters{
        simplex[K] theta;
        real<lower=0,upper=1> us[K];
    }
    model{
        real ps[K];
        for (n in 1:N) {
          for (k in 1:K) {
            ps[k] <- log(theta[k]) + binomial_log(Y[n], K, us[k]);
          }
          increment_log_prob(log_sum_exp(ps));
        }
    }
"

mixdata <- list(K=3, N= 5, Y = c(1, 2, 3, 2, 1))
testfit <- stan(model_code=mix_model, data=mixdata, iter=10)
fit     <- stan(fit=testfit, data=mixdata, iter=25000, chains=5)
