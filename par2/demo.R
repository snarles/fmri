## Generate data and try to guess mu, sigma2

library(pracma)

objective_function <- function(accs, ks, mu, tau, mc.reps = 1e4) {
  xs <- ((1:mc.reps) - 0.5)/mc.repsy
  ts <- qnorm(xs, mean = mu, sd = sqrt(tau))
  Ts <- repmat(t(ts), length(ks), 1)
  pts <- repmat(t(pnorm(ts)), length(ks), 1)
  Ks <- repmat(t(t(ks-1)), 1, length(ts))
  fs <- rowMeans(pts^Ks)
  sum((accs - fs)^2)
}

objective_grad <- function(accs, ks, mu, tau, mc.reps = 1e4) {
  xs <- ((1:mc.reps) - 0.5)/mc.reps
  ts <- qnorm(xs, mean = mu, sd = sqrt(tau))
  Ts <- repmat(t(ts), length(ks), 1)
  pts <- repmat(t(pnorm(ts)), length(ks), 1)
  Ks <- repmat(t(t(ks-1)), 1, length(ts))
  fs <- rowMeans(pts^Ks)
  dmu <- rowMeans(pts^Ks * (Ts - mu)/tau)
  dtau <- rowMeans(pts^Ks * ((mu-Ts)^2 - tau)/(2 * tau^2))
  c(sum(2 * (fs - accs) * dmu),
    sum(2 * (fs - accs) * dtau))
}

objective_hess <- function(accs, ks, mu, tau, mc.reps = 1e4) {
  xs <- ((1:mc.reps) - 0.5)/mc.reps
  ts <- qnorm(xs, mean = mu, sd = sqrt(tau))
  Ts <- repmat(t(ts), length(ks), 1)
  pts <- repmat(t(pnorm(ts)), length(ks), 1)
  Ks <- repmat(t(t(ks-1)), 1, length(ts))
  fs <- rowMeans(pts^Ks)
  dmu <- rowMeans(pts^Ks * (Ts - mu)/tau)
  dtau <- rowMeans(pts^Ks * ((mu-Ts)^2 - tau)/(2 * tau^2))
  d2mu <- rowMeans(pts^Ks * (-1/tau + ((Ts - mu)/tau)^2))
  d2mu_tau <- rowMeans(pts^Ks *((mu - Ts)/(tau^2) + (Ts - mu)/tau * ((mu-Ts)^2 - tau)/(2 * tau^2)))
  d2tau <- rowMeans(pts^Ks * ((tau - 2 * (mu - Ts)^2)/(2 * tau^3) + (((mu-Ts)^2 - tau)/(2 * tau^2))^2))
  h11 <- 2 * sum(dmu^2 + (fs - accs) * d2mu)
  h12 <- 2 * sum(dmu*dtau + (fs - accs) * d2mu_tau)
  h22 <- 2 * sum(dtau^2 + (fs - accs) * d2tau)
  matrix(c(h11, h12, h12, h22), 2, 2)
}

attr(objective_function, "gradient") <- objective_grad
attr(objective_function, "hessian") <- objective_hess

par2_fit <- function(ks, accs, mc.reps = 1e4) {
    if (min(accs) == 1) {
        return(c(Inf,0))
    }
    mu.init <- par2_initialize_mu(accs, ks)
    tau.init <- 1
    res <- tryCatch({
        nlm(function(x) par2_objective_function(accs, ks, x[1], x[2]),
            c(mu.init, tau.init))},
        error = function(e) e)
    if ("error" %in% class(res)) {
        lala <- tryCatch({
            nlm(function(x) par2_objective_function(accs, ks, x[1], x[2]),
                c(mu.init, tau.init))},
            warning = function(e) e)
        if("sqrt(tau)" %in% lala) {
            res2 <- optimize(function(x) par2_objective_function(accs, ks, x, 0), interval = c(-20, 20))
            (mu.hat <- res2$minimum)
            tau.hat <- 0      
        } else {
            print(lala)
            stop("unknown error")
        }    
    } else {
        (mu.hat <- res$estimate[1])
        (tau.hat <- res$estimate[2])
    }
    return(c(mu.hat, tau.hat))
}

## data settings
(mu <- 5 * rexp(1))
(sigma2 <- exp(rnorm(1)/2))
tau <- sigma2
n.reps <- 100
ks <- 2:100
K <- max(ks)

z_stars <- rnorm(n.reps, mu, sqrt(sigma2))
zmat <- pracma::randn(n.reps, K)
zmaxs <- t(apply(zmat, 1, cummax))[, ks - 1]
dim(zmaxs)
accs <- colMeans(zmaxs < z_stars)
plot(ks, accs, type = "l")

par2_objective_function(accs, ks, mu, sigma2)
par2_objective_function(accs, ks, mu - 0.1, sigma2 + 0.1)
(fits <- par2_fit(ks, accs))
c(mu, sigma2)
par2_objective_function(accs, ks, fits[1], fits[2])


plot(ks, accs, type = "l")
lines(ks, par2_acc_k(ks, fits[1], fits[2]), col = "red")

par(mfcol = c(1,2))
plot(accs, par2_acc_k(ks, fits[1], fits[2]))
abline(0,1)

plot(accs, par2_acc_k(ks, mu, sigma2))
abline(0,1)


numDeriv::grad(function(x) objective_function(accs, ks, x[1], x[2]),
               c(mu, sigma2))

objective_grad(accs, ks, mu, sigma2)

numDeriv::hessian(function(x) objective_function(accs, ks, x[1], x[2]),
               c(mu, sigma2))
objective_hess(accs, ks, mu, sigma2)

nlm(function(x) objective_function(accs, ks, x[1], x[2]),
    c(mu, sigma2))