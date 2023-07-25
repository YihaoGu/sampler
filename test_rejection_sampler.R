#### THIS CODE IS MESSY. BE CAREFUL. ####

### Test rejection sampler
source('rejection_sampler_negative.R')

## Visually compare the empirical distribution to the target.
compute_target_pdf <- function(x, beta, alpha, normalized = TRUE) {
  log_prior <- - log(1 + x)
  loglik <- - beta * x + alpha * sqrt(x)
  logp <- loglik + log_prior
  logp <- logp - max(logp)
  # Avoid numerical under-flow when exponentiating.
  pdf <- exp(logp)
  if (normalized) {
    pdf <- pdf / trapz(pdf, x)
  }
  return(pdf)
}

trapz <- function(f, x) {
  dx <- tail(x, -1) - head(x, -1)
  return(
    sum((head(f, -1) + tail(f, -1)) * dx / 2)
  )
}


## Sample the horseshoe local scale conditionally on fixed beta and global scale.
alpha <- -0.99
# m1 = (1 - sqrt(1 - alpha^2)) / 2
# beta <- m1 / 100
beta <- 0.1

n_samples <- 10^4
lscale_samples <-
  sapply(1:n_samples, function(i) rejection_sampler(beta, alpha))
1 / mean(lscale_samples)

# Restrict the plot range; otherwise, the empirical 
# distribution of a heavy-tailed target is unstable.
max_quantile <- .99
upper_xlim <- quantile(lscale_samples, max_quantile)
breaks <- seq(0, upper_xlim, length.out = 100)

x <- seq(0, upper_xlim, length.out = 10^5)[-1]
pdf <- compute_target_pdf(x, beta, alpha)
max_pdf = max(pdf)

hist(
  lscale_samples[lscale_samples < upper_xlim], 
  breaks = breaks, col='blue',
  # xlab = bquote(lambda ~ "|" ~ beta / tau ~ "=" ~ .(beta / gscale)),
  xlab = 'local scale sample',
  # ylim = c(0,max_pdf*1.1),
  main=NULL,
  prob=TRUE
)

lines(x, pdf, col = 'orange', lw=2)

legend(
  "topright", 
  legend = c('empirical dist', 'target density'),
  lty = c(2, 1), lwd=c(10, 2), bty='n',
  col = c('blue', 'orange')
)

fixed_beta_ar = function(beta){
  helper = function(alpha){
    lscale_samples <-
      sapply(1:n_samples, function(i) rejection_sampler(beta, alpha))
    1 / mean(lscale_samples)
  }
  return(helper)
}

fixed_alpha_ar = function(alpha){
  helper = function(beta){
    lscale_samples <-
      sapply(1:n_samples, function(i) rejection_sampler(beta, alpha))
    1 / mean(lscale_samples)
  }
  return(helper)
}

beta = 1
breakpoints = seq(-1, -10^(-4), length.out=10)
breakpoints = matrix(breakpoints, nrow = 1)
ar = apply(breakpoints, 2, fixed_beta_ar(beta))
plot(breakpoints, ar, type = 'l', ylim = c(0,1))

alpha = -1
# m1 = (1 - sqrt(1 - alpha^2)) / 2
# m2 = (1 + sqrt(1 - alpha^2)) / 2
breakpoints = seq(1, 10, length.out=100)
breakpoints = matrix(breakpoints, nrow = 1)
ar = apply(breakpoints, 2, fixed_alpha_ar(alpha))
plot(breakpoints, ar, type = 'l', ylim = c(0,1), 
     ylab = 'Acceptance Rate', xlab = expression(beta))
mtext(expression(paste(alpha, " = - 1")), side = 3)

# breakpoints2 = seq(1, 10, length.out=30)
# breakpoints2 = matrix(breakpoints2, nrow = 1)
# ar2 = apply(breakpoints2, 2, fixed_alpha_ar(alpha))
# plot(breakpoints2, ar2, type = 'l', ylim = c(0,1), 
#      ylab = 'Acceptance Rate', xlab = expression(beta))
# mtext(expression(paste(alpha, " = - 0.3")), side = 3)
# 
# breakpoints = cbind(breakpoints, breakpoints2)
# ar = cbind(ar, ar2)
# plot(breakpoints, ar, type = 'l', ylim = c(0,1), 
#      ylab = 'Acceptance Rate', xlab = expression(beta))
# mtext(expression(paste(alpha, " = - 0.9")), side = 3)

alpha = -0.1
a2 = (2 - alpha^2 + 2 * sqrt(1-alpha^2)) / alpha^2
m1 = (1 - sqrt(1 - alpha^2)) / 2
beta = m1 / 100
b = alpha^2 / beta^2 / 4
rtrunc(1, spec = 'gamma', a = a2, b = b, shape = 2, rate = -alpha)