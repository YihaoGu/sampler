## to do:
  # avoid repetitive computing
  # ...

# Compute \int_0^m exp(alpha * sqrt(x)) dx
exp_sqrt_term = function(alpha, m){
  x = alpha * sqrt(m)
  x = exp(x) * (x - 1) + 1
  x = x * 2 / (alpha^2)
  # if (is.na(x)){
  #   x = m
  # }
  # cat('exp_sqrt_term:', x, '\n')
  return(x)
}

# Compute \int_a^b exp(h(a) + (x-a) * s) dx, where s = (h(b) - h(a)) / (b - a)
secant_term = function(beta, alpha, a, b){
  if (a == b){
    return(0)
  }
  fa = exp(- beta * a + alpha * sqrt(a)) / (1 + a)
  fb = exp(- beta * b + alpha * sqrt(b)) / (1 + b)
  ha = -log(1 + a) - beta * a + alpha * sqrt(a)
  hb = -log(1 + b) - beta * b + alpha * sqrt(b)
  s = (hb - ha) / (b - a)
  if (s != 0){
    s = (fb - fa) / s
  }
  return(s)
}

slope = function(beta, alpha, a, b){
  ha = -log(1 + a) - beta * a + alpha * sqrt(a)
  hb = -log(1 + b) - beta * b + alpha * sqrt(b)
  s = (hb - ha) / (b - a)
  return(s)
}

# Compute \int_a^{+\infty} exp(h(a) - beta * (x - a)) dx
exp_term = function(beta, alpha, a){
  fa = exp(- beta * a + alpha * sqrt(a)) / (1 + a)
  fa = fa / beta
  # cat('exp_term:', fa, '\n')
  return(fa)
}
#############################################
# Sample from a density propotional to 1 / (1 + x) * exp(-beta * x + alpha * sqrt(x)) 
# alpha is negative

rejection_sampler = function(beta, alpha){
  accepted = FALSE
  i = 0
  while (!accepted) {
    i = i + 1
    rv = sample_from_bounding_dist(beta, alpha)
    target_logp = -log(1 + rv) - beta * rv + alpha * sqrt(rv)
    bound_logp = compute_bounding_logp(rv, beta, alpha)
    accept_prob = exp(target_logp - bound_logp)
    accepted = accept_prob > runif(1)
  }
  return(i)
}

sample_from_bounding_dist = function(beta, alpha){
  if ((alpha <= -1) || (beta >= ((1 + sqrt(1 - alpha^2)) / 2))){
  # if ((alpha <= -1) || (beta >= 1)){
    #cat('Easy Case! \n')
    rv = sample_from_easy_case(beta, alpha)
  }
  else if (beta <= ((1 - sqrt(1 - alpha^2)) / 2)){
  # else if (beta <= (- alpha / 4)){
    #cat('alpha > -1 & beta small \n')
    rv = sample_from_alpha_small_beta_small(beta, alpha)
  }
  else{
    #cat('alpha > -1 & beta medium \n')
    rv = sample_from_alpha_small_beta_medium(beta, alpha)
  }
  return(rv)
}

# if alpha <= -1 or beta => (1 + sqrt(1-alpha^2)) / 2
sample_from_easy_case = function(beta, alpha){
  alpha_sq = alpha ^ 2
  t = alpha_sq / beta
  mass_to_left = (exp(- t / 2) * (- t / 2 - 1) + 1) * 2 / alpha_sq
  mass_to_right = exp(- 3 / 4 * t) / (beta + t / 4)
  prob_to_left = mass_to_left / (mass_to_left + mass_to_right)
  sampled_from_left = (prob_to_left > runif(1))
  if (sampled_from_left){
    rv = rtrunc(1, spec = 'gamma', a = 0, b = - alpha / beta / 2, shape = 2, rate = - alpha)
    rv = rv^2
  }
  else{
    rv = rexp(1, rate = beta) + alpha_sq / beta^2 / 4
  }
  return(rv)
}

# if -1 <= alpha < 0 and 0 < beta <= (1 - sqrt(1-alpha^2)) / 2
sample_from_alpha_small_beta_small = function(beta, alpha){
  r = 1 / alpha^2
  t = 2 * sqrt(r^2 - r)
  r = 2 * r
  a1 = r - t - 1
  a2 = r + t - 1
  b = alpha^2 / beta^2 / 4
  # int1 = exp_sqrt_term(alpha, a1)
  int1 = secant_term(beta, alpha, 0, a1)
  int2 = exp(-beta * a1 + alpha * sqrt(a1)) * log((r+t)/(r-t))
  int3 = secant_term(beta, alpha, a2, b)
  int4 = exp_term(beta, alpha, b)
  int_all = sum(int1, int2, int3, int4)
  u = runif(1)
  if (u < (int1 / int_all)){
    # rv = rtrunc(1, spec = 'gamma', a = 0, b = sqrt(a1), shape = 2, rate = - alpha)
    # rv = rv^2
    s = slope(beta, alpha, 0, a1)
    rv = rtrunc(1, spec = 'exp', a = 0 , b = a1, rate = - s)
  }
  else if (u < (sum(int1, int2) / int_all)){
    # s = slope(beta, alpha, a1, a2)
    # rv = rtrunc(1, spec = 'exp', a = 0, b = a2 - a1, rate = - s) + a1
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

# sample_from_alpha_small_beta_small = function(beta, alpha){
#   a = 4
#   b = alpha^2 / beta^2 / 4
#   int1 = secant_term(beta, alpha, 0, a)
#   int2 = (exp_sqrt_term(alpha, b) - exp_sqrt_term(alpha, a)) / (1 + a) / exp(beta * a)
#   int3 = exp_term(beta, alpha, b)
#   int_all = sum(int1, int2, int3)
#   u = runif(1)
#   if (u < (int1 / int_all)){
#     s = slope(beta, alpha, 0, a)
#     rv = rtrunc(1, spec = 'exp', a = 0, b = a, rate = -s)
#   }
#   else if (u < (1 - int3 / int_all)){
#     rv = rtrunc(1, spec = 'gamma', a = a, b = b, shape = 2, rate = -alpha)
#     rv = rv^2
#   }
#   else{
#     rv = rexp(1, rate = beta) + b
#   }
#   return(rv)
# }

sample_from_alpha_small_beta_medium = function(beta, alpha){
  r = 1 / alpha^2
  t = 2 * sqrt(r^2 - r)
  r = 2 * r
  a = r - t - 1
  b = 1 / beta - 1
  int1 = secant_term(beta, alpha, 0, a)
  int2 = exp(-beta * a + alpha * sqrt(a)) * log((1+b)/(1+a))
  int3 = exp_term(beta, alpha, b)
  int_all = sum(int1, int2, int3)
  u = runif(1)
  if (u < (int1 / int_all)){
    s = slope(beta, alpha, 0, a)
    rv = rtrunc(1, spec = 'exp', a = 0 , b = a, rate = - s)
  }
  else if (u < (sum(int1, int2) / int_all)){
    u = runif(1, min = log(1+a), max = log(1+b))
    rv = exp(u) - 1
  }
  else{
    rv = rexp(1, rate = beta) + b
  }
  return(rv)
}
  
# sample_from_alpha_small_beta_medium = function(beta, alpha){
#   thresh = alpha^2 / beta^2 / 4
#   int1 = secant_term(beta, alpha, 0, thresh)
#   int2 = exp_term(beta, alpha, thresh)
#   int_all = sum(int1, int2)
#   u = runif(1)
#   if (u < (int1 / int_all)){
#     s = slope(beta, alpha, 0, thresh)
#     rv = rtrunc(1, spec = 'exp', a = 0, b = thresh, rate = - s)
#   }
#   else{
#     rv = rexp(1, rate = beta) + thresh
#   }
#   return(rv)
# }

compute_bounding_logp = function(x, beta, alpha){
  # if ((alpha <= -1) || (beta >= ((1 + sqrt(1 - alpha^2)) / 2))){
  if ((alpha <= -1) || (beta >= 1)){
    logp = easy_case_logp(x, beta, alpha)
  }
  # else if (beta <= ((1 - sqrt(1 - alpha^2)) / 2)){
  else if (beta <= (- alpha / 4)){
    logp = alpha_small_beta_small_logp(x, beta, alpha)
  }
  else{
    logp = alpha_small_beta_medium_logp(x, beta, alpha)
  }
  return(logp)
}

easy_case_logp = function(x, beta, alpha){
  thresh = alpha^2 / beta^2 / 4
  if (x < thresh){
    logp = alpha * sqrt(x)
  }
  else{
    logp = - log(1 + thresh) - beta * thresh + alpha * sqrt(thresh) - beta * x
  }
  return(logp)
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
    logp = - log(1 + a2) - beta * a2 + alpha * sqrt(a2) + slope(beta, alpha, a2, b) * (x - a2)
  }
  else{
    logp = - log(1 + b) - beta * b + alpha * sqrt(b) - beta * x
  }
  return(logp)
}

# alpha_small_beta_small_logp = function(x, beta, alpha){
#   a = 4
#   b = alpha^2 / beta^2 / 4
#   if (x < a){
#     logp = slope(beta, alpha, 0, a) * x
#   }
#   else if (x < b){
#     logp = -log(1 + a) - beta * a + alpha * sqrt(x)
#   }
#   else{
#     logp = - log(1 + b) - beta * b + alpha * sqrt(b) - beta * x
#   }
#   return(logp)
# }

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
    logp = - log(1 + x) - beta * a + alpha * sqrt(a)
  }
  else{
    logp = - log(1 + b) - beta * b + alpha * sqrt(b) - beta * x
  }
  return(logp)
}

# alpha_small_beta_medium_logp = function(x, beta, alpha){
#   thresh = alpha^2 / beta^2 / 4
#   if (x < thresh){
#     logp = slope(beta, alpha, 0, thresh) * x
#   }
#   else{
#     logp = - log(1 + thresh) - beta * thresh + alpha * sqrt(thresh) - beta * x
#   }
#   return(logp)
# }

