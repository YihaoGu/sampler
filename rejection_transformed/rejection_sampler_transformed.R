## Sample from a density proportional to exp(-a^2 * (sqrt(e^(2x)-1) - c)^2)
# a > 0, c \in \mathbb{R}

rejection_sampler_transformed = function(a, c, Nsim){
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
      if (c >= 0){
        target_logp = -a^2 * (sqrt(exp(2 * rv) - 1) - c)^2
      }
      else{
        target_logp = -a^2 * (exp(2 * rv) - 1) + 2 * a^2 * c * sqrt(exp(2 * rv) - 1)
      }
      accept_prob = exp(target_logp - bound_logp)
      # if (accept_prob > 1){
      #   cat(accept_prob, '\n')
      # }
      accepted = accept_prob > runif(1)
    }
    RV[idx] = rv
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
  if (c >= 0){
    if (a^2 * c^2 > log(4)){
      if (c <= 3 * sqrt(3)){
        case = 1
      }
      else{
        case = 2
      }
    }
    else if (a^2 * c^2 > log(2)){
      case = 3
    }
    else{
      case = 4
    }
  }
  else{
    case = 5
  }
  return(case)
}

# Choose the sampler for different cases
choose_sampler = function(case){
  sampler_set = c(sample_from_case_1, 
                  sample_from_case_2,
                  sample_from_case_3,
                  sample_from_case_4,
                  sample_from_case_5)
  sampler = sampler_set[[case]]
  return(sampler)
}

# Case 1: a^2 * c^2 > log(4) && c <= 3 * sqrt(3)

sample_from_case_1 = function(a, c){
  m0 = c(c - sqrt(log(4)) / a, c - sqrt(log(2)) / a, c + sqrt(log(2)) / a, c + sqrt(log(4)) / a)
  m = log(1 + m0^2) / 2
  d1 = 2 * a * sqrt(log(4)) * (1 + m0[1]^2) / m0[1]
  int1 = 0.25 / d1 * (1 - exp(-d1 * m[1]))
  int2 = (m[2] - m[1]) / 2
  int3 = m[3] - m[2]
  int4 = (m[4] - m[3]) / 2
  d2 = 2 * a * sqrt(log(4)) * (1 + m0[4]^2) / m0[4]
  int5 = 0.25 / d2
  ints = c(int1, int2, int3, int4, int5)
  probs = cumsum(ints / sum(ints))
  u = runif(1)
  
  if (u < probs[1]){
    rv = runif(1, exp(-d1 * m[1]), 1)
    rv = m[1] + log(rv) / d1
    logp = -log(4) + d1 * (rv - m[1])
  }
  else if (u < probs[2]){
    rv = runif(1, m[1], m[2])
    logp = -log(2)
  }
  else if (u < probs[3]){
    rv = runif(1, m[2], m[3])
    logp = 0
  }
  else if (u < probs[4]){
    rv = runif(1, m[3], m[4])
    logp = -log(2)
  }
  else{
    rv = m[4] + rexp(1, rate = d2)
    logp = -log(4) - d2 * (rv - m[4])
  }
  res = list('rv' = rv, 'logp' = logp)
  return(res)
}

# Case 2: a^2 * c^2 > log(4) && c > 3 * sqrt(3)
sample_from_case_2 = function(a, c){
  r = compute_inflection_point(c)
  m0 = c(c - sqrt(log(4)) / a, c - sqrt(log(2)) / a, c + sqrt(log(2)) / a, c + sqrt(log(4)) / a)
  if (r[2] > m0[1]){
    r[2] = m0[1]
  }
  m0 = c(r, m0)
  m = log(1 + m0^2) / 2
  
  d1 = 2 * a * sqrt(log(4)) * (1 + m0[1]^2) / m0[1]
  int1 = exp(-a^2 * (m0[1] - c)^2) / d1 * (1 - exp(-d1 * m[1]))
  
  d2 = - a^2 * ((m0[2] - c)^2 - (m0[1] - c)^2)/ (m[2] - m[1])
  int2 = exp(-a^2 * (m0[2] - c)^2) / d2 * (1 - exp(d2 * (m[1] - m[2])))
  
  d3 = 2 * a * sqrt(log(4)) * (1 + m0[3]^2) / m0[3]
  int3 = 0.25 / d3 * (1 - exp(d3 * (m[2] - m[3])))
  
  int4 = (m[4] - m[3]) / 2
  int5 = m[5] - m[4]
  int6 = (m[6] - m[5]) / 2
  
  d4 = 2 * a * sqrt(log(4)) * (1 + m0[6]^2) / m0[6]
  int7 = 0.25 / d4
  
  ints = c(int1, int2, int3, int4, int5, int6, int7)
  probs = cumsum(ints / sum(ints))
  u = runif(1)
  
  if (u < probs[1]){
    rv = runif(1, exp(-d1 * m[1]), 1)
    rv = m[1] + log(rv) / d1
    logp = -a^2 * (m0[1] - c)^2 + d1 * (rv - m[1])
  }
  else if (u < probs[2]){
    rv = runif(1, exp(d2 * (m[1] - m[2])), 1)
    rv = m[2] + log(rv) / d2
    logp = -a^2 * (m0[2] - c)^2 + d2 * (rv - m[2])
  }
  else if (u < probs[3]){
    # cat('piece3\n')
    rv = runif(1, exp(d3 * (m[2] - m[3])), 1)
    rv = m[3] + log(rv) / d3
    logp = -log(4) + d3 * (rv - m[3])
  }
  else if (u < probs[4]){
    rv = runif(1, m[3], m[4])
    logp = -log(2)
  }
  else if (u < probs[5]){
    rv = runif(1, m[4], m[5])
    logp = 0
  }
  else if (u < probs[6]){
    rv = runif(1, m[5], m[6])
    logp = -log(2)
  }
  else{
    rv = m[6] + rexp(1, rate = d4)
    logp = -log(4) - d4 * (rv - m[6])
  }
  res = list('rv' = rv, 'logp' = logp)
  return(res)
}

compute_inflection_point = function(c){
  if (c > 3 * sqrt(3)){
    theta = acos(3 * sqrt(3) / c)
    t1 = cos(theta / 3) / sqrt(3)
    t2 = sin(theta / 3)
    r = c(t1 + t2, t1 - t2)
  }
  else if (c < - 3 * sqrt(3)){
    theta = acos(3 * sqrt(3) / c)
    t1 = cos(theta / 3) / sqrt(3)
    t2 = sin(theta / 3)
    r = t1 + t2
  }
  else{
    c_reciprocal =  1 / c
    delta = sqrt(c_reciprocal^2 - 1 / 27)
    r = (- c_reciprocal + delta)^(1/3) + (- c_reciprocal - delta)^(1/3)
  }
  return(1 / r)
}

# Case 3: log(2) < a^2 * c^2 <= log(4)

sample_from_case_3 = function(a, c){
  m0 = c(c - sqrt(log(2)) / a, c + sqrt(log(2)) / a, c + sqrt(log(4)) / a)
  m = log(1 + m0^2) / 2
  int1 = m[1] / 2
  int2 = m[2] - m[1]
  int3 = (m[3] - m[2]) / 2
  d = 2 * a * sqrt(log(4)) * (1 + m0[3]^2) / m0[3]
  int4 = 0.25 / d
  ints = c(int1, int2, int3, int4)
  probs = cumsum(ints / sum(ints))
  u = runif(1)
  
  if (u < probs[1]){
    rv = runif(1, 0, m[1])
    logp = -log(2)
  }
  else if (u < probs[2]){
    rv = runif(1, m[1], m[2])
    logp = 0
  }
  else if (u < probs[3]){
    rv = runif(1, m[2], m[3])
    logp = -log(2)
  }
  else{
    rv = m[3] + rexp(1, rate = d)
    logp = -log(4) - d * (rv - m[3])
  }
  res = list('rv' = rv, 'logp' = logp)
  return(res)
}

# Case 4: a^2 * c^2 <= log(2)

sample_from_case_4 = function(a, c){
  m0 = c(c + sqrt(log(2)) / a, c + sqrt(log(4)) / a)
  m = log(1 + m0^2) / 2
  int1 = m[1]
  int2 = (m[2] - m[1]) / 2
  d = 2 * a * sqrt(log(4)) * (1 + m0[2]^2) / m0[2]
  int3 = 0.25 / d
  ints = c(int1, int2, int3)
  probs = cumsum(ints / sum(ints))
  u = runif(1)
  
  if (u < probs[1]){
    rv = runif(1, 0, m[1])
    logp = 0
  }
  else if (u < probs[2]){
    rv = runif(1, m[1], m[2])
    logp = -log(2)
  }
  else{
    rv = m[2] + rexp(1, rate = d)
    logp = -log(4) - d * (rv - m[2])
  }
  res = list('rv' = rv, 'logp' = logp)
  return(res)
}

# Case 5
sample_from_case_5 = function(a, c){
  m0 = c(c + sqrt(a^2 * c^2 + log(2)) / a, c + sqrt(a^2 * c^2 + log(4)) / a)
  r = compute_inflection_point(c)
  m0[3] = max(m0[2], r)
  m = log(1 + m0^2) / 2
  
  int1 = m[1]
  int2 = (m[2] - m[1]) / 2
  if (r <= m0[2]){
    int3 = 0
    d1 = 2 * a * sqrt(a^2 * c^2 + log(4)) * (1 + m0[2]^2) / m0[2]
    int4 = 0.25 / d1
  }
  else{
    d0 = (a^2 * (m0[3]^2 - m0[2]^2) - 2 * a^2 * c * (m0[3] - m0[2])) / (m[3] - m[2])
    int3 = 0.25 / d0 * (1 - exp(d0 * (m[2] - m[3])))
    d1 = 2 * a^2 * (m0[3] - c) * (1 + m0[3]^2) / m0[3]
    int4 = exp(-a^2 * m0[3]^2 + 2 * a^2 * c * m0[3]) / d1
  }
  ints = c(int1, int2, int3, int4)
  probs = cumsum(ints / sum(ints))
  u = runif(1)
  
  if (u < probs[1]){
    rv = runif(1, 0, m[1])
    logp = 0
  }
  else if (u < probs[2]){
    rv = runif(1, m[1], m[2])
    logp = -log(2)
  }
  else if (u < probs[3]){
    rv = runif(1, 0, 1 - exp(d0 * (m[2] - m[3])))
    rv = m[2] - log(1 - rv) / d0
    logp = -log(4) - d0 * (rv - m[2])
  }
  else{
    rv = m[3] + rexp(1, rate = d1)
    if (m0[3] <= m0[2]){
      logp = -log(4) - d1 * (rv - m[3])
    }
    else{
      logp = -a^2 * m0[3]^2 + 2 * a^2 * c * m0[3] - d1 * (rv - m[3])
    }
  }
  res = list('rv' = rv, 'logp' = logp)
  return(res)
}