### Helper Functions ###

## Helper functions for the rejection sampler
# f(x)
f = function(x, beta, alpha){
  exp(- beta * x + alpha * sqrt(x)) / (1 + x)
}

# h(x)
h = function(x, beta, alpha){
  -beta * x + alpha * sqrt(x) - log(1 + x)
}

# Compute \int_0^m exp(alpha * sqrt(x)) dx
exp_sqrt_term = function(alpha, m){
  x = alpha * sqrt(m)
  x = exp(x) * (x - 1) + 1
  x = x * 2 / (alpha^2)
  return(x)
}

# Compute the slope of the line segment between (a, h(a)) and (b, h(b))
slope = function(beta, alpha, a, b){
  ha = h(a, beta, alpha)
  hb = h(b, beta, alpha)
  s = (hb - ha) / (b - a)
  return(s)
}

# Compute \int_a^b exp(h(a) + (x-a) * s) dx, where s = (h(b) - h(a)) / (b - a)
secant_term = function(beta, alpha, a, b){
  if (a == b){
    return(0)
  }
  fa = f(a, beta, alpha)
  fb = f(b, beta, alpha)
  s = slope(beta, alpha, a, b)
  if (s != 0){
    s = (fb - fa) / s
  }
  return(s)
}

# Compute \int_a^b exp(-beta * a + alpha * sqrt(a))/ (1+x) dx
log_term = function(beta, alpha, a, b){
  exp(-beta * a + alpha * sqrt(a)) * log((1 + b) / (1 + a))
}

# Compute \int_a^{+\infty} exp(h(a) - beta * (x - a)) dx
exp_term = function(beta, alpha, a){
  fa = f(a, beta, alpha)
  fa = fa / beta
  return(fa)
}

################################################
## Helper functions for the test file
# Compute the target pdf
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

# Trapz method
trapz <- function(f, x) {
  dx <- tail(x, -1) - head(x, -1)
  return(
    sum((head(f, -1) + tail(f, -1)) * dx / 2)
  )
}

# Visually compare the empirical distribution and the target
plot_lscale_sample <- function(lscale_sample, beta, alpha){
  # Restrict the plot range; otherwise, the empirical 
  # distribution of a heavy-tailed target is unstable.
  max_quantile <- .99
  upper_xlim <- quantile(lscale_sample, max_quantile)
  breaks <- seq(0, upper_xlim, length.out = 100)
  
  x <- seq(0, upper_xlim, length.out = 10^5)[-1]
  pdf <- compute_target_pdf(x, beta, alpha)
  max_pdf = max(pdf)
  
  hist(
    lscale_sample[lscale_sample < upper_xlim], 
    breaks = breaks, col='blue',
    xlab = bquote(x ~ "|" ~ beta ~ "=" ~ .(beta) ~ "," ~ alpha ~ "=" ~ .(alpha)),
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
}

# Fix alpha and compute the acceptance rate for different betas
fixed_alpha_ar = function(alpha){
  helper = function(beta){
    res = rejection_sampler_negative(beta, alpha, Nsim = 10^4)
    return(res$ar)
  }
  return(helper)
}

# Fix beta and compute the acceptance rate for different alphas
fixed_beta_ar = function(beta){
  helper = function(alpha){
    res = rejection_sampler_negative(beta, alpha, Nsim = 10^4)
    return(res$ar)
  }
  return(helper)
}

