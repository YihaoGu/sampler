## Sample from a density proportional to x / (1 + x^2) * exp(-a^2 * (x-c)^2)
# c is positive

rejection_sampler_positive = function(a, c, Nsim){
  RV = rep(0, Nsim)
  AR = rep(0, Nsim)
  for (idx in 1 : Nsim){
    accepted = FALSE
    i = 0
    while (!accepted) {
      i = i + 1
      rv = sample_from_bounding_dist(a, c)
      target_logp = log(rv) - log(1 + rv^2) - a^2 * (rv-c)^2
      bound_logp = compute_bounding_logp(rv, a, c)
      accept_prob = exp(target_logp - bound_logp)
      accepted = accept_prob > runif(1)
    }
    RV[idx] = rv
    AR[idx] = i
  }
  res = list('sample' = RV, 'ar' = 1 / mean(AR))
  return(res)
}

# Sample from bounding distribution
sample_from_bounding_dist = function(a, c){
  rv = sample_from_a_small_c_small(a, c)
  return(rv)
}

# Compute the bounding logp
compute_bounding_logp = function(x, a, c){
  logp = a_small_c_small_logp(x, a, c)
  return(logp)
}

# Case 1: o < a < sqrt(pi) / 2 and 0 < c < 1 / a
sample_from_a_small_c_small = function(a, c){
  m = (sqrt(pi) + sqrt(pi - 4 * a^2)) / 2 / a
  sd = sqrt(0.5) / a
  mass_to_left = 0.5 * log(1 + m^2)
  mass_to_right = 1 - pnorm(m, mean = c, sd = sd)
  prob_to_left = mass_to_left / (mass_to_left + mass_to_right)
  if (prob_to_left > runif(1)){
    v = runif(1, min = 0, max = mass_to_left)
    rv = sqrt(exp(2*v) - 1)
  }
  else{
    rv = rtruncnorm(1, a = m, b = Inf, mean = c, sd = sd)
  }
  return(rv)
}

a_small_c_small_logp = function(x, a, c){
  m = (sqrt(pi) + sqrt(pi - 4 * a^2)) / 2 / a
  if (x < m){
    logp = log(x) - log(1 + x^2)
  }
  else{
    logp = log(a / sqrt(pi)) - a^2 * (x - c)^2
  }
  return(logp)
}