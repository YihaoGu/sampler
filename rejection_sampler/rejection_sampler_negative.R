## to do:
  # avoid repetitive computing
  # ...

## Sample from a density proportional to 1 / (1 + x) * exp(-beta * x + alpha * sqrt(x)) 
# alpha is negative

# The  rejection sampler
rejection_sampler_negative = function(beta, alpha, Nsim){
  RV = rep(0, Nsim)
  AR = rep(0, Nsim)
  for (idx in 1 : Nsim){
    accepted = FALSE
    i = 0
    while (!accepted) {
      i = i + 1
      rv = sample_from_bounding_dist(beta, alpha)
      bound_logp = compute_bounding_logp(rv, beta, alpha)
      target_logp = -log(1 + rv) - beta * rv + alpha * sqrt(rv)
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
sample_from_bounding_dist = function(beta, alpha){
  if ((alpha <= -1) || (beta >= ((1 + sqrt(1 - alpha^2)) / 2))){
    # cat('Easy Case! \n')
    rv = sample_from_easy_case(beta, alpha)
  }
  else if (beta <= ((1 - sqrt(1 - alpha^2)) / 2)){
    # cat('alpha > -1 & beta small \n')
    rv = sample_from_alpha_small_beta_small(beta, alpha)
  }
  else{
    # cat('alpha > -1 & beta medium \n')
    rv = sample_from_alpha_small_beta_medium(beta, alpha)
  }
  return(rv)
}

# Compute the bounding logp
compute_bounding_logp = function(x, beta, alpha){
  if ((alpha <= -1) || (beta >= ((1 + sqrt(1 - alpha^2)) / 2))){
    logp = easy_case_logp(x, beta, alpha)
  }
  else if (beta <= ((1 - sqrt(1 - alpha^2)) / 2)){
    logp = alpha_small_beta_small_logp(x, beta, alpha)
  }
  else{
    logp = alpha_small_beta_medium_logp(x, beta, alpha)
  }
  return(logp)
}

# Case 1: alpha <= -1 or beta => (1 + sqrt(1-alpha^2)) / 2
sample_from_easy_case = function(beta, alpha){
  m = alpha^2 / beta^2 / 4
  mass_to_left = exp_sqrt_term(alpha, m)
  mass_to_right = exp_term(beta, alpha, m)
  prob_to_left = mass_to_left / (mass_to_left + mass_to_right)
  sampled_from_left = (prob_to_left > runif(1))
  if (sampled_from_left){
    rv = rtrunc(1, spec = 'gamma', a = 0, b = - alpha / beta / 2, shape = 2, rate = - alpha)
    rv = rv^2
  }
  else{
    rv = rexp(1, rate = beta) + m
  }
  return(rv)
}

easy_case_logp = function(x, beta, alpha){
  m = alpha^2 / beta^2 / 4
  if (x < m){
    logp = alpha * sqrt(x)
  }
  else{
    logp = - log(1 + m) + alpha * sqrt(m) - beta * x
  }
  return(logp)
}

# Case 2: -1 < alpha < 0 and beta <= (1 - sqrt(1-alpha^2)) / 2
sample_from_alpha_small_beta_small = function(beta, alpha){
  r = 1 / alpha^2
  t = 2 * sqrt(r^2 - r)
  r = 2 * r
  a1 = r - t - 1
  a2 = r + t - 1
  b = alpha^2 / beta^2 / 4
  int1 = secant_term(beta, alpha, 0, a1)
  int2 = log_term(beta, alpha, a1, a2)
  int3 = secant_term(beta, alpha, a2, b)
  int4 = exp_term(beta, alpha, b)
  int_all = sum(int1, int2, int3, int4)
  u = runif(1)
  if (u < (int1 / int_all)){
    s = slope(beta, alpha, 0, a1)
    rv = rtrunc(1, spec = 'exp', a = 0 , b = a1, rate = - s)
  }
  else if (u < (sum(int1, int2) / int_all)){
    v = runif(1, min = log(1+a1), max = log(1+a2))
    rv = exp(v) - 1
  }
  else if (u < (sum(int1, int2, int3) / int_all)){
    s = slope(beta, alpha, a2, b)
    rv = rtrunc(1, spec = 'exp', a = 0, b = b - a2, rate = - s) + a2
  }
  else{
    rv = rexp(1, rate = beta) + b
  }
  return(rv)
}

alpha_small_beta_small_logp = function(x, beta, alpha){
  r = 1 / alpha^2
  t = 2 * sqrt(r^2 - r)
  r = 2 * r
  a1 = r - t - 1
  a2 = r + t - 1
  b = alpha^2 / beta^2 / 4
  if (x < a1){
    logp = slope(beta, alpha, 0, a1) * x
  }
  else if (x < a2){
    logp = - beta * a1 + alpha * sqrt(a1) - log(1 + x)
  }
  else if (x < b){
    logp = h(a2, beta, alpha) + slope(beta, alpha, a2, b) * (x - a2)
  }
  else{
    logp = - log(1 + b) + alpha * sqrt(b) - beta * x
  }
  return(logp)
}

#  Case 3: -1 < alpha < 0 and (1 - sqrt(1-alpha^2)) / 2 < beta < (1 + sqrt(1-alpha^2)) / 2
sample_from_alpha_small_beta_medium = function(beta, alpha){
  r = 1 / alpha^2
  t = 2 * sqrt(r^2 - r)
  r = 2 * r
  a = r - t - 1
  b = 1 / beta - 1
  int1 = secant_term(beta, alpha, 0, a)
  int2 = log_term(beta, alpha, a, b)
  int3 = exp_term(beta, alpha, b)
  int_all = sum(int1, int2, int3)
  u = runif(1)
  if (u < (int1 / int_all)){
    s = slope(beta, alpha, 0, a)
    rv = rtrunc(1, spec = 'exp', a = 0 , b = a, rate = - s)
  }
  else if (u < (sum(int1, int2) / int_all)){
    v = runif(1, min = log(1+a), max = log(1+b))
    rv = exp(v) - 1
  }
  else{
    rv = rexp(1, rate = beta) + b
  }
  return(rv)
}

alpha_small_beta_medium_logp = function(x, beta, alpha){
  r = 1 / alpha^2
  t = 2 * sqrt(r^2 - r)
  r = 2 * r
  a = r - t - 1
  b = 1 / beta - 1
  if (x < a){
    logp = slope(beta, alpha, 0, a) * x
  }
  else if (x < b){
    logp = - beta * a + alpha * sqrt(a)- log(1 + x) 
  }
  else{
    logp = - log(1 + b) + alpha * sqrt(b) - beta * x 
  }
  return(logp)
}

