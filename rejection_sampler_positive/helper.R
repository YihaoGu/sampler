### Helper Functions ###

## Helper functions for the test file
# Compute the target pdf
compute_target_pdf <- function(x, a, c, normalized = TRUE) {
  log_prior <- log(x) - log(1 + x^2)
  loglik <- - a^2 * (x- c)^2 
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
plot_lscale_sample <- function(lscale_sample, a, c){
  # Restrict the plot range; otherwise, the empirical 
  # distribution of a heavy-tailed target is unstable.
  max_quantile <- .99
  upper_xlim <- quantile(lscale_sample, max_quantile)
  breaks <- seq(0, upper_xlim, length.out = 100)
  
  x <- seq(0, upper_xlim, length.out = 10^5)[-1]
  pdf <- compute_target_pdf(x, a, c)
  max_pdf = max(pdf)
  
  hist(
    lscale_sample[lscale_sample < upper_xlim], 
    breaks = breaks, col='blue',
    xlab = bquote(x ~ "|" ~ a ~ "=" ~ .(a) ~ "," ~ c ~ "=" ~ .(c)),
    ylim = c(0, 1.1 * max_pdf), 
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

# Fix a and compute the acceptance rate for different c
fixed_a_ar = function(a){
  helper = function(c){
    res = rejection_sampler_positive(a, c, Nsim = 10^4)
    return(res$ar)
  }
  return(helper)
}

# Fix c and compute the acceptance rate for different a
fixed_c_ar = function(c){
  helper = function(a){
    res = rejection_sampler_positive(a, c, Nsim = 10^4)
    return(res$ar)
  }
  return(helper)
}

