## Sample from a density proportional to x / (1 + x^2) * exp(-a^2 * (x-c)^2)
# a > 0, c \in \mathbb{R}

rejection_sampler = function(a, c, Nsim){
  a = abs(a)
  case = classifier(a, c)
  cat('Case', case, '\n')
  sampler = choose_sampler(case)
  
  RV = rep(0, Nsim)
  AR = rep(0, Nsim)
  for (idx in 1 : Nsim){
    accepted = FALSE
    i = 0
    while (!accepted) {
      i = i + 1
      res = sampler(a, c)
      rv = res$rv
      bound_logp = res$logp
      target_logp = compute_target_logp(rv, a, c, case)
      accept_prob = exp(target_logp - bound_logp)
      accepted = accept_prob > runif(1)
    }
    RV[idx] = transform_rv(rv, case)
    AR[idx] = i
  }
  
  ret = list('sample' = RV, 'ar' = 1 / mean(AR))
  return(ret)
}

# Categorize inputs a and c into one of the eight cases 
classifier = function(a, c){
  if (a == 0){
    stop('a cannot be 0!')
  }
  else if ((a < 1) && (c >= 0) && (c < 1 / a)){
    case = 1
  }
  else if ((a < 0.25) && (c >= 1 / a)){
    case = 2
  }
  else if ((a >= 0.25) && (c >= max(1 / a, 1))){
    case = 3
  }
  else if ((a >= 1) && (c >= 0.5)){
    case = 4
  }
  else if ((a >= 1) && (c >= 0) && (c < 0.5)){
    case = 5
  }
  else if ((a < 1) && (c >= -0.5/a^2)){
    case = 6
  }
  else {
    case = 7
  }
  return(case)
}

# Choose the sampler for different cases
choose_sampler = function(case){
  sampler_set = c(sample_from_case_1, 
                  sample_from_case_2,
                  sample_from_case_3,
                  sample_from_case_4,
                  sample_from_case_5,
                  sample_from_case_6,
                  sample_from_case_7)
  sampler = sampler_set[[case]]
  return(sampler)
}

# Compute log probability of the target density
compute_target_logp = function(x, a, c, case){
  if (case == 7){
    logp = - log(1 + x) - a^2 * x + 2 * a^2 * c * sqrt(x)
  }
  else{
    logp = log(x) - log(1 + x^2) - a^2 * (x - c)^2
  }
  return(logp)
}

# Transform the random variable if necessary
transform_rv = function(x, case){
  if (case == 7){
    x = sqrt(x)
  }
  return(x)
}


# Case 1: 0 < a <= 1 and 0 <= c <= 1 / a
# g(x) = x / (1 + x^2)                    if 0 <= x < m,
#        0.5 * a * exp(-a^2 * (x - c)^2)  if x >= m,
# where m = (1 + sqrt(1 - a^2)) / a.

sample_from_case_1 = function(a, c){
  m = (1 + sqrt(1 - a^2)) / a
  sd = sqrt(0.5) / a
  mass_to_left = 0.5 * log(1 + m^2)
  mass_to_right = pnorm(m, mean = c, sd = sd, lower.tail = FALSE) * sqrt(pi) / 2
  
  prob_to_left = mass_to_left / (mass_to_left + mass_to_right)
  if (prob_to_left > runif(1)){
    v = runif(1, min = 0, max = mass_to_left)
    rv = sqrt(exp(2*v) - 1)
    logp = log(rv) - log(1 + rv^2)
  }
  else{
    rv = rtrunc(1, spec = 'norm', a = m, b = Inf, mean = c, sd = sd)
    logp = - log(2) + log(a) - a^2 * (rv - c)^2
  }
  res = list('rv' = rv, 'logp' = logp)
  return(res)
}

# Case 2: 0 < a < 0.25, c >= 1 / a
# g(x) = x / (1 + x^2) * exp(-9 / 16 * a^2 * c^2) if 0 <=x < m1,
#        m1 / (1 + m1^2) * exp(-a^2 * (x - c)^2)  if m1 <= x < m2,
#        m2 / (1 + m2^2) * exp(-a^2 * (x - c)^2)  if m2 <= x < m3,
#        m3 / (1 + m3^2) * exp(-a^2 * (x - c)^2)  if x >= m3,
# where m1 = c / 4, m2 = c / 2, m3 = c.

sample_from_case_2 = function(a, c){
  m1 = c / 4
  m2 = c / 2
  m3 = c
  sd = sqrt(0.5) / a
  # int1 = 0.5 * log(1 + m1^2)* exp(-9 / 16 * a^2 * c^2)
  logint1 = -log(2) + log(log(1 + m1^2)) - 9 / 16 * a^2 * c^2
  int1 = exp(logint1)
  int2 = (pnorm(m2, mean = c, sd = sd) - pnorm(m1, mean = c, sd = sd)) * sqrt(pi) / a * m1 / (1 + m1^2)
  int3 = (pnorm(m3, mean = c, sd = sd) - pnorm(m2, mean = c, sd = sd)) * sqrt(pi) / a * m2 / (1 + m2^2)
  int4 = 0.5 * sqrt(pi) / a * m3 / (1 + m3^2)
  ints = c(int1, int2, int3, int4)
  probs = cumsum(ints / sum(ints))
  u = runif(1)
  
  if (u < probs[1]){
    v = runif(1, min = 0, max = 0.5 * log(1 + m1^2))
    rv = sqrt(exp(2*v) - 1)
    logp = log(rv) - log(1 + rv^2) - 9 / 16 * a^2 * c^2
  }
  else if (u < probs[2]){
    rv = rtrunc(1, spec = 'norm', a = m1, b = m2, mean = c, sd = sd)
    logp = log(m1) - log(1 + m1^2) - a^2 * (rv - c)^2
  }
  else if (u < probs[3]){
    rv = rtrunc(1, spec = 'norm', a = m2, b = m3, mean = c, sd = sd)
    logp = log(m2) - log(1 + m2^2) - a^2 * (rv - c)^2
  }
  else{
    rv = rtrunc(1, spec = 'norm', a = m3, b = Inf, mean = c, sd = sd)
    logp = log(m3) - log(1 + m3^2) - a^2 * (rv - c)^2
  }
  res = list('rv' = rv, 'logp' = logp)
  return(res)
}

# Case 3: 0.25 <= a < 1, c >= 1 / a or a >= 1, c >= 1
# g(x) = 0.5 * exp(-a^2 * (x - c)^2)             if 0 <= x < m,
#        m / (1 + m^2) * exp(-a^2 * (x - c)^2)   if x >= m,
# where m = max{1, c - log(c) / sqrt(2 * a)}.

sample_from_case_3 = function(a, c){
  sd = sqrt(0.5) / a
  m = max(1, c - log(c) * sd)
  int1 = (pnorm(m, mean = c, sd = sd) - pnorm(0, mean = c, sd = sd)) * sqrt(pi) / a / 2
  int2 = pnorm(m, mean = c, sd = sd, lower.tail = FALSE) * sqrt(pi) / a * m / (1 + m^2)
  prob = int1 / (int1 + int2)
  u = runif(1)
  
  if (u < prob){
    rv = rtrunc(1, spec = 'norm', a = 0, b = m, mean = c, sd = sd)
    logp = -log(2) - a^2 * (rv - c)^2
  }
  else{
    rv = rtrunc(1, spec = 'norm', a = m, b = Inf, mean = c, sd = sd)
    logp = log(m) - log(1 + m^2) - a^2 * (rv - c)^2
  }
  res = list('rv' = rv, 'logp' = logp)
  return(res)
}

# Case 4: a >= 1, 0.5 <= c < 1
# g(x) = 0.5 * exp(- a^2 * (x - c)^2)

sample_from_case_4 = function(a, c){
  sd = sqrt(0.5) / a
  rv = rtrunc(1, spec = 'norm', a = 0, b = Inf, mean = c, sd = sd)
  logp = -log(2) - a^2 * (rv - c)^2
  res = list('rv' = rv, 'logp' = logp)
  return(res)
}

# Case 5: a >= 1, 0 <= c < 0.5
# g(x) = c / (1 + c^2) * exp(- a^2 * (x - c)^2)   if 0 <= x < c, 
#        x * exp(-a^2 * (x - c)^2),               if x >= c.

sample_from_case_5 = function(a, c){
   sd = sqrt(0.5) / a
   int1 = (0.5 - pnorm(0, mean = c, sd = sd)) * sqrt(pi) / a * c / (1 + c^2)
   int2 = 0.5 / a^2 + 0.5 * sqrt(pi) / a * c
   prob = int1 / (int1 + int2)
   u = runif(1)
   
   if (u <= prob){
     rv = rtrunc(1, spec = 'norm', a = 0, b = c, mean = c, sd = sd)
     logp = log(c) - log(1 + c^2) - a^2 * (rv - c)^2
   }
   else{
     v = runif(1)
     if (v <= (1 / (1 + a * c * sqrt(pi)))){
       w = runif(1)
       rv = sqrt(-2 * log(1 - w)) * sd + c
     }
     else{
       rv = rtrunc(1, spec = 'norm', a = c, b = Inf, mean = c, sd = sd)
     }
     logp = log(rv) - a^2 * (rv - c)^2
   }
   
   res = list('rv' = rv, 'logp' = logp)
   return(res)
}

# Case 6: 0 < a < 0.5, -0.5 / a^2 <= c < -1 / a
# g(x) = x / (1 + x^2) * exp(-a^2 * c^2)    if x < m1,
#        1 / (1 + x^2) * m * exp(-a^2 * (m - c)^2)  if m1 <= x < m2,
#        m2 / (1 + m2^2) * exp(-a^2 * (x - c)^2)    if x >= m2,
# where m = (a * c + sqrt(a^2 * c^2 + 2)) / a / 2, 
#       m1 = m * exp(- a^2 * (m^2 - 2 * m * c)), 
#       m2 = max(c + 0.5 / a^2, 1).

sample_from_case_6 = function(a, c){
  
  m = (a * c + sqrt(a^2 * c^2 + 2)) / a / 2
  m1 = m * exp(- a^2 * (m^2 - 2 * m * c))
  m2 = max(c + 0.5 / a^2, 1)
  sd = sqrt(0.5) / a
  
  logint1 = log(0.5) + log(log(1 + m1^2))
  logint2 = log((atan(m2) - atan(m1))) + log(m1)
  logint3 = pnorm(m2, mean = c, sd = sd, log = TRUE, lower.tail = FALSE) + 0.5* log(pi) - log(a) + log(m2) - log(1+m2^2) + a^2 * c^2
  logint = c(logint1, logint2, logint3) - max(logint1, logint2, logint3)
  ints = exp(logint)
  prob = ints / sum(ints)
  u1 = prob[1]
  u2 = prob[1] + prob[2]
  u = runif(1)
  
  if (u < u1){
    v = runif(1, min = 0, max = log(1 + m1^2))
    rv = sqrt(exp(v) - 1)
    logp = log(rv) - log(1 + rv^2) - a^2 * c^2
  }
  else if (u < u2){
    v = runif(1, min = atan(m1), max = atan(m2))
    rv = tan(v)
    logp = - log(1 + rv^2) + log(m1) - a^2 * c^2
  }
  else{
    rv = rtrunc(1, spec = 'norm', a = m2, b = Inf, mean = c, sd = sd)
    logp = log(m2) - log(1 + m2^2) - a^2 * (rv - c)^2
  }
  
  res = list('rv' = rv, 'logp' = logp)
  return(res)
}
 
# Case 7: 0 < a <= 1, c < -0.5 / a^2, or a >= 1
# beta = alpha^2, alpha = 2 * a^2 * c, alpha^2 / beta ^2 / 4 = c^2
# transformation: x <- x^2
# f(x) = 1 / (1 + x) * exp(- beta * x + alpha * sqrt(x))
# g(x) = exp(alpha * sqrt(x))           if x < c^2
#        f(c^2) * exp(-beta*(x - c^2))  if x >= c^2
 
sample_from_case_7 = function(a, c){
  beta = a^2
  alpha = 2 * beta * c
  m = c^2
  sqrt_m = abs(c)
  mass_to_left = 2 / alpha^2 * (exp(alpha * sqrt_m) * (alpha * sqrt_m - 1) + 1)
  mass_to_right = 1 / (1 + m) * exp(- beta * m + alpha * sqrt_m) / beta
  prob = mass_to_left / (mass_to_left + mass_to_right)
  u = runif(1)
 
  if (u < prob){
    rv = rtrunc(1, spec = 'gamma', a = 0, b = - alpha / beta / 2, shape = 2, rate = - alpha)
    rv = rv^2
    logp = alpha * sqrt(rv)
  }
  else{
    rv = rexp(1, rate = beta) + m
    logp = - log(1 + m) + alpha * sqrt_m - beta * rv
  }
  res = list('rv' = rv, 'logp' = logp)
  return(res)
}

