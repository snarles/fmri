####
##  Porting JVHW estimator
####

library(R.matlab)
library(pracma)

fingerprint <- function(samp) {
  tt <- table(table(samp))
  vs <- as.numeric(names(tt))
  ff <- numeric(max(vs))
  ff[vs] <- tt
  ff
}

h_mle <- function(samp) {
  freqs <- table(samp)
  ps <- freqs/sum(freqs)
  -sum(ps * log(ps, base=2))
}

## block 0
# est_entro_MLE(samp) % 5.3105
# est_entro_JVHW(samp)
# sum(samp) % 58421
# sum(samp.^2) % 4248073


source("entropy//temp_samp.R")
h_mle(samp)
sum(samp)
sum(samp^2)

## block 1
# [n, wid] = size(samp)
# order = 4 + ceil(1.2*log(n))
# V = [0.3303  -0.3295  0.4679]

(n <- length(samp))
(order <- 4 + ceiling(1.2*log(n)))
V <- c(0.3303, -0.3295, 0.4679)

## block2
# load poly_coeff_r.mat poly_coeff_r;
# coeff = poly_coeff_r{order}
# num2cell(coeff)

poly_coeff_r <- readMat('entropy//poly_coeff_r.mat')
(coeff <- poly_coeff_r[[1]][order][[1]][[1]])

## block 3
# f = find([diff(sort(samp)); ones(1,wid,class(samp))]) % f: linear index of the last occurrence in each set of repeated values along the columns of samp
# f = accumarray({[f(1);diff(f)],ceil(f/n)},1)   % f: fingerprint
# prob = (1:size(f,1))/n

(ff <- fingerprint(samp))
(prob <- 1:length(ff)/n)

## block 4
# f1nonzero = find(f(1,:) > 0)
# lenf1nonzero = length(f1nonzero)
# c_1 = zeros(1, wid)
# if n > 15 && lenf1nonzero >0
# [ log(n) * ones(1,lenf1nonzero); log(f(1,f1nonzero)); ones(1,lenf1nonzero)]
# c_1(f1nonzero) = V * [ log(n) * ones(1,lenf1nonzero); log(f(1,f1nonzero)); ones(1,lenf1nonzero)];
# c_1 = max(c_1, 1/(1.9*log(n))); % make sure threshold is higher than 1/n
# end
# c_1

c_1 <- 0
if (n > 15 && ff[1] > 0) {
  c_1 <- sum(V * c(log(n), log(ff[1]), 1))
  c_1 <- pmax(c_1, 1/(1.9*log(n)))
}
c_1

## block 5
# x = prob;
# g_coeff = coeff;

x <- prob
g_coeff <- coeff

## block 6
# K = length(g_coeff) - 1   % K: the order of best polynomial approximation, g_coeff = {g0, g1, g2, ..., g_K}
# thres = 4*c_1*log(n)/n
# output = zeros(length(x), length(c_1));
# [thres, x] = meshgrid(thres,x)   

K <- length(g_coeff) - 1
(thres = 4*c_1*log(n)/n)
output = zeros(length(x), length(c_1));

## block 7
# region_large = x>thres
# region_nonlarge = ~region_large;
# region_mid = x>thres/2 & region_nonlarge
# output(region_large) = -x(region_large) .* log(x(region_large)) + 1/(2*n);   

region_large = x>thres
region_nonlarge = !region_large
region_mid = x>thres/2 & region_nonlarge
output[region_large] = -x[region_large] * log(x[region_large]) + 1/(2*n);   
t(output)

## block 8
# x1(:,1) = x(region_nonlarge)  
# thres1(:,1) = thres(region_nonlarge) % Ensure x1 and thres1 are column vectors
# q = 0:K-1
# bsxfun(@minus, n*x1, q)
# bsxfun(@times, thres1, n-q)
# [thres1, bsxfun(@minus, n*x1, q)./bsxfun(@times, thres1, n-q)]
# temp = cumprod([thres1, bsxfun(@minus, n*x1, q)./bsxfun(@times, thres1, n-q)],2)
# temp * g_coeff.'
#         temp * g_coeff.'- x1.*log(thres1)
# output(region_nonlarge) = cumprod([thres1, bsxfun(@minus, n*x1, q)./bsxfun(@times, thres1, n-q)],2)*g_coeff.'- x1.*log(thres1);

x1 <- x[region_nonlarge]
q = 0:(K-1)
n * x1
repmat(n * t(t(x1)), 1, K) - repmat(q, length(x1), 1)
thres * repmat(n - q, length(x1), 1)
temp0 <- cbind(thres, (repmat(n * t(t(x1)), 1, K) - repmat(q, length(x1), 1))/(thres * repmat(n - q, length(x1), 1)))
temp <- t(apply(temp0, 1, cumprod))
output[region_nonlarge] <- temp %*% t(g_coeff) - x1 * log(thres)


## block 9
#     ratio = 2*x(region_mid)./thres(region_mid) - 1
#     output(region_mid) = ratio.*(-x(region_mid) .* log(x(region_mid)) + 1/(2*n)) + (1-ratio).*output(region_mid); 
#     output = max(output,0);


ratio = 2*x[region_mid]/thres - 1
output[region_mid] = ratio*(-x[region_mid] * log(x[region_mid]) + 1/(2*n)) + 
  (1-ratio)*output[region_mid]

t(output)

## block 10

(est = sum(ff*output)/log(2))
