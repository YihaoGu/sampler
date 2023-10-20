## Sample from a density proportional to exp(-a^2 * (sqrt(e^(2x)-1) - c)^2).

# Input Arguments:
  # a > 0, c \in \mathbb{R}
  # Nsim: size of the sample
  # q: min f(x) / max f(x) on each constant piece
  # k1: use 1, q, q^2, ..., q^(k1-1) as constant proposal
  # k2: divide the convex segment into k2 pieces

# Output Arguments: a list which contains:
  # sample: the sample
  # ar: the acceptance rate obtained via simulation

rejection_sampler_transformed = function(a, c, Nsim, q = 0.5, k1 = 2, k2 = 1){
  a = abs(a)

  RV = rep(0, Nsim)
  AR = rep(0, Nsim)
  for (idx in 1 : Nsim){
    accepted = FALSE
    i = 0
    while (!accepted) {
      i = i + 1
      if (c >= 0){
        res = sample_from_c_positive(a, c, q, k1, k2)
      }
      else{
        res = sample_from_c_negative(a, c, q, k1, k2)
      }
      rv = res$rv
      bound_logp = res$logp
      if (c >= 0){
        target_logp = -a^2 * (sqrt(exp(2 * rv) - 1) - c)^2
      }
      else{
        target_logp = -a^2 * (exp(2 * rv) - 1) + 2 * a^2 * c * sqrt(exp(2 * rv) - 1)
      }
      accept_prob = exp(target_logp - bound_logp)
      accepted = accept_prob > runif(1)
    }
    RV[idx] = rv
    AR[idx] = i
  }
  
  ret = list('sample' = RV, 'ar' = 1 / mean(AR))
  return(ret)
}

# Case 1: c positive
## c positive case
# target density: exp(-a^2 * (sqrt(exp(2*x) - 1) - c)^2)
# q: min f(x) / max f(x) on each constant piece
# k1: use 1, q, q^2, ..., q^(k1-1) as constant proposal
# k2: divide the convex segment into k2 pieces

sample_from_c_positive = function(a, c, q, k1, k2){
  
  # Compute and store repeatedly used constants.
  log_q = log(q)
  sqrt_logq = sqrt(-log_q)
  a_sq = a^2
  
  # Compute how many constant pieces we can use on the left side of the mode.
  # If f(0) >= q^k1, we can use piecewise constant proposal to upper bound the 
  # left side of the mode with guaranteed acceptance rate; otherwise, we need
  # to compute where the target density is concave/convex and use exponential proposal.
  k0 = - a^2 * c^2 / log_q
  if (k0 <= k1){
    r_computed = FALSE
    k0 = floor(k0) + 1
  }
  else{
    r_computed = TRUE
    k0 = k1
  }
  
  # Compute the breakpoints determined by constant proposal densities: 
  # Suppose m are the breakpoints,
  # let m_transformed = sqrt(exp(2 * m) - 1).
  m_transformed = c + sqrt_logq / a * c(- sqrt(k0:0), sqrt(1:k1))
  if (! r_computed){
    m_transformed[1] = 0
  }
  m = log(1 + m_transformed^2) / 2
  
  # Compute the integrals of these constant proposal densities.
  log_densities_constant = log_q * abs(-k0:k1)
  densities_constant = exp(log_densities_constant)
  int_constant = c(densities_constant[2:(k0+1)] * (m[2:(k0+1)] - m[1:k0]),
                   densities_constant[(k0+1):(k0+k1)] * (m[(k0+2):(k0+k1+1)] - m[(k0+1):(k0+k1)]))
  
  # We first deal with the right tail, which is much more tractable.
  # In the right tail, the target density is always log-concave, so we use an
  # exponential proposal to upper bound it.
  d_concave_opposite = 2 * a_sq * (1 + m_transformed[k0+k1+1]^2) * (1 - c / m_transformed[k0+k1+1])
  int_concave_right = densities_constant[k0+k1+1] / d_concave_opposite
  # cat(d_concave_opposite,'\n')
  
  # Now we start to deal with the left "tail" if r_computed == TRUE.
  int_concave_left = 0
  int_convex = rep(0, k2)
  int_concave_middle = 0
  
  if (r_computed){
    r_transformed = compute_inflection_point(c)
    r = log(1 + r_transformed^2) / 2
    
    if (length(r) == 0){
      d_concave_left = 2 * a * sqrt(abs(log_densities_constant[1])) * (1 / m_transformed[1] + m_transformed[1])
      coef_concave_left = 1 - exp(- d_concave_left * m[1])
      int_concave_left = densities_constant[1] / d_concave_left * coef_concave_left
    }
    else{
      
      if (r[1] >= m[1]){
        r[1] = m[1]
        r_transformed[1] = m_transformed[1]
      }
      
      # concave part
      tmp = c - r_transformed[1]
      d_concave_left = 2 * a_sq * tmp * (1 / r_transformed[1] + r_transformed[1])
      log_density_concave_left = - a_sq * tmp^2
      density_concave_left = exp(log_density_concave_left)
      coef_concave_left = 1 - exp(- d_concave_left * r[1])
      int_concave_left = density_concave_left / d_concave_left * coef_concave_left
      
      # convex part
      if (min(r[2], m[1]) != r[1]){
        m_convex = seq(r[1], min(r[2], m[1]), length.out = k2 + 1)
        tmp = sqrt(exp(2 * m_convex) - 1) - c
        log_densities_convex = -a_sq * tmp^2
        densities_convex = exp(log_densities_convex)
        d_convex = (log_densities_convex[2:(k2+1)] - log_densities_convex[1:k2]) / (m_convex[2:(k2+1)] - m_convex[1:k2])
        coef = 1 - exp(d_convex * (m_convex[1:k2] - m_convex[2:(k2+1)]))
        int_convex = densities_convex[2:(k2+1)] / d_convex * coef
      }
      
      # concave part
      if (r[2] < m[1]){
        d_concave_middle = 2 * a * sqrt(abs(log_densities_constant[1])) * (1 / m_transformed[1] + m_transformed[1])
        coef_concave_middle = 1 - exp(d_concave_middle * (r[2] - m[1]))
        int_concave_middle = densities_constant[1] / d_concave_middle * coef_concave_middle
      }
    }
  }
  
  # Concatenate all the pieces.
  ints = c(int_concave_left, int_convex, int_concave_middle, int_constant, int_concave_right)
  probs = ints / sum(ints)
  # cat(int_concave_left,'\n', int_convex,'\n',int_concave_middle,'\n',int_constant,'\n',int_concave_right,'\n')
  
  # Determine from which piece we sample from.
  piece = sample(length(probs), 1, prob = probs) # length(probs) = 1 + k2 + 1 + k0 + k1 + 1
  
  # Sample from the chosen piece and compute the corresponding log density.
  if (piece == 1){
    rv = runif(1, 0, coef_concave_left)
    if (length(r) == 0){
      rv = m[1] + log(1 - rv) / d_concave_left
      logp = log_densities_constant[1] + d_concave_left * (rv - m[1])
    }
    else{
      rv = r[1] + log(1 - rv) / d_concave_left
      logp = log_density_concave_left + d_concave_left * (rv - r[1])
    }
  }
  else if (piece <= (1 + k2)){
    piece = piece - 1
    rv = runif(1, 0, coef[piece])
    rv = m_convex[piece+1] + log(1 - rv) / d_convex[piece]
    logp = log_densities_convex[piece+1] + d_convex[piece] * (rv - m_convex[piece+1])
  }
  else if (piece <= (2 + k2)){
    rv = runif(1, 0, coef_concave_middle)
    rv = m[1] + log(1 - rv) / d_concave_middle
    logp = log_densities_constant[1] + d_concave_middle * (rv - m[1])
  }
  else if (piece <= (2 + k2 + k0 + k1)){
    piece = piece - 2 - k2
    rv = runif(1, m[piece], m[piece + 1])
    if (piece <= k0){
      logp = log_densities_constant[piece + 1]
    }
    else{
      logp = log_densities_constant[piece]
    }
  }
  else{
    rv = m[k0+k1+1] + rexp(1, rate = d_concave_opposite)
    logp = log_densities_constant[k0+k1+1] - d_concave_opposite * (rv - m[k0+k1+1])
  }
  
  # Return the results.
  res = list('rv' = rv, 'logp' = logp)
  return(res)
}

# Case 2: c negative
## c negative case
# target density: exp(-a^2 * (exp(2*x) - 1) + 2 * a^2 * c * sqrt(exp(2*x) - 1))
# q: min f(x) / max f(x) on each constant piece
# k1: use 1, q, q^2, ..., q^(k1-1) as constant proposal
# k2: divide the convex segment into k2 pieces

sample_from_c_negative = function(a, c, q, k1, k2){
  
  # Compute and store repeatedly used constants.
  log_q = log(q)
  a_sq = a^2
  
  # Compute the breakpoints determined by constant proposal densities: 
  # Suppose m are the breakpoints,
  # let m_transformed = sqrt(exp(2 * m) - 1), and we have
  # exp(-a^2 * m_transformed^2 + 2 * a^2 * c * m_transformed) = 1 / 2 ^ (0:k1).
  m_transformed = sqrt(c^2 - log_q / a_sq * (0:k1)) + c
  # if (m_transformed[k1+1] == 0){
  #   m_transformed = - (0:k1) * neg_log_q / 2 / a_sq / c
  # }
  m = log(1 + m_transformed^2) / 2
  
  # Now we compute the inflection points where the target density
  # changes from being convex to being concave
  r_transformed = compute_inflection_point(c)
  
  # If r_transformed > m_transformed[k1+1], we keep it;
  # otherwise, it is not used so we do not need to compute r.
  r_computed = FALSE
  if (r_transformed > m_transformed[k1+1]){
    r_computed = TRUE
    r_transformed_sq = r_transformed^2
    r = log(1 + r_transformed_sq) / 2
  }
  
  # Compute the integrals of these constant proposal densities.
  log_densities = log_q * (0:k1)
  densities = exp(log_densities)
  int_constant = densities[1:k1] * (m[2:(k1+1)] - m[1:k1])
  
  # If r is not computed, we only need to compute the last exponential piece;
  # otherwise, we need to initialize on both the convex and concave intervals.
  # Note that we use the opposite values of the derivatives here.
  if (! r_computed){
    # convex part
    int_convex = rep(0, k2)
    
    # concave part
    d_concave_opposite = 2 * a_sq * (1 + m_transformed[k1+1]^2) * (1 - c / m_transformed[k1+1])
    int_concave = densities[k1+1] / d_concave_opposite
  }
  else{
    # convex part
    m_convex = seq(m[k1+1], r, length.out = k2 + 1)
    tmp = exp(2 * m_convex) - 1
    log_densities_convex = - a_sq * tmp + 2 * a_sq * c * sqrt(tmp)
    densities_convex = exp(log_densities_convex)
    d_convex_opposite = - (log_densities_convex[2:(k2+1)] - log_densities_convex[1:k2]) / (m_convex[2:(k2+1)] - m_convex[1:k2])
    coef = 1 - exp(d_convex_opposite * (m_convex[1:k2] - m_convex[2:(k2+1)]))
    int_convex = densities_convex[1:k2] / d_convex_opposite * coef
    
    # concave part
    d_concave_opposite = 2 * a_sq * (1 + r_transformed_sq) * (1 - c / r_transformed)
    int_concave = densities_convex[k2+1] / d_concave_opposite
  }
  
  # Concatenate all the pieces.
  ints = c(int_constant, int_convex, int_concave)
  probs = ints / sum(ints)
  
  # Determine from which piece we sample from.
  piece = sample(k1 + k2 + 1, 1, prob = probs)
  
  # Sample from the chosen piece and compute the corresponding log density.
  if (piece <= k1){
    rv = runif(1, m[piece], m[piece + 1])
    logp = log_densities[piece]
  }
  else if (piece <= (k1 + k2)){
    # Use the inverse transform method
    piece = piece - k1
    rv = runif(1, 0, coef[piece])
    rv = m_convex[piece] - log(1 - rv) / d_convex_opposite[piece]
    logp = log_densities_convex[piece] - d_convex_opposite[piece] * (rv - m_convex[piece])
  }
  else{
    if (! r_computed){
      rv = m[k1+1] + rexp(1, rate = d_concave_opposite)
      logp = log_densities[k1+1] - d_concave_opposite * (rv - m[k1+1])
    }
    else{
      rv = r + rexp(1, rate = d_concave_opposite)
      logp = log_densities_convex[k2+1] - d_concave_opposite * (rv - r)
    }
  }
  
  # Return the results
  res = list('rv' = rv, 'logp' = logp)
  return(res)
}

# Compute the (transformed) points where the target density changes from 
# being concave/convex to being convex/concave
compute_inflection_point = function(c){
  if ((c >= 0) && (c <= 3 * sqrt(3))){
    r = NULL
  }
  else if (c > (3 * sqrt(3))){
    theta = acos(3 * sqrt(3) / c)
    t1 = cos(theta / 3) / sqrt(3)
    t2 = sin(theta / 3)
    r = c(t1 + t2, t1 - t2)
    r = 1 / r
  }
  else if (c < - 3 * sqrt(3)){
    theta = acos(3 * sqrt(3) / c)
    t1 = cos(theta / 3) / sqrt(3)
    t2 = sin(theta / 3)
    r = t1 + t2
    r = 1 / r
  }
  else{
    c_reciprocal =  1 / c
    delta = sqrt(c_reciprocal^2 - 1 / 27)
    r = (- c_reciprocal + delta)^(1/3) + (- c_reciprocal - delta)^(1/3)
    r = 1 / r
  }
  return(r)
}