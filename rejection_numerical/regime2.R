# library('truncnorm')
# 
# # Monte Carlo Integration for f(x) on interval [m1, m2]
# mc_integral = function(a, b, m1, m2){
#   sd = sqrt(1 / 2) / a
#   X = rtruncnorm(10^4, a = m1, b = m2, mean=b, sd=sd)
#   Y = X / (1 + X^2) * sqrt(pi) / a * (pnorm(m2, mean=b, sd=sd) - pnorm(m1, mean=b, sd=sd))
#   mean(Y)
# }

get_acceptance_rate = function(a, b) {
  target = function(x) {
    x / (1 + x ^ 2) * exp(-a ^ 2 * (x - b) ^ 2)
  }
  
  mu = b
  std = 1 / (a * sqrt(2))
  
  m1 = max(1, b / 4)
  m2 = max(1, b / 2)
  m3 = b
  
  # breaking integrals into more parts helps numerical stability,
  # also makes it easy to check which region contributes most mass
  f_area1 = integrate(target, lower = 0, upper = m1)$value
  f_area2 = integrate(target, lower = m1, upper = m2)$value
  f_area3 = integrate(target, lower = m2, upper = m3)$value
  f_area4 = integrate(target, lower = m3, upper = Inf)$value
  
  # use monte carlo integration
  # f_area1 = mc_integral(a, b, 0, m1)
  # f_area2 = mc_integral(a, b, m1, m2)
  # f_area3 = mc_integral(a, b, m2, m3)
  # f_area4 = mc_integral(a, b, m3, Inf)
  
  f_area = f_area1 + f_area2 + f_area3 + f_area4
  stopifnot(f_area > 1e-10) # warning for underflow
  
  zero_to_m1 = 0.5 * log(1 + m1^2) * exp(-a^2 * (0.75 * b)^2)
  m1_to_m2 = (pnorm(m2, mean = mu, sd = std) -  pnorm(m1, mean = mu, sd = std)) / a * (sqrt(pi)) * m1 / (1 + m1 ^ 2)
  m2_to_m3 = (pnorm(m3, mean = mu, sd = std) - pnorm(m2, mean = mu, sd = std)) / a * (sqrt(pi)) * m2 / (1 + m2 ^ 2)
  m3_to_inf = (1 - pnorm(m3, mean = mu, sd = std)) / a * (sqrt(pi)) * m3 / (1 + m3 ^ 2)
  g_area =  zero_to_m1 + m1_to_m2 + m2_to_m3 + m3_to_inf

  # cat('Mass:', zero_to_m1, m1_to_m2, m2_to_m3, m3_to_inf, '\n')
  # cat('Acceptance rate:', f_area1 / zero_to_m1, f_area2 / m1_to_m2, f_area3 / m2_to_m3, f_area4 / m3_to_inf, '\n')
  return(f_area / g_area)
}

# a = 0.25
# b = c(1/a, 2/a, 10/a, 100/a, 1000/a)
# sapply(b, get_acceptance_rate, a = a)
# 
# a = 0.1
# b = c(1/a, 2/a, 10/a, 100/a, 1000/a)
# sapply(b, get_acceptance_rate, a = a)
# 
# a = 0.01
# b = c(1/a, 2/a, 10/a, 100/a, 1000/a)
# sapply(b, get_acceptance_rate, a = a)
# 
# a = 0.001
# b = c(1/a, 2/a, 10/a, 100/a, 1000/a)
# sapply(b, get_acceptance_rate, a = a)
# 
# a = 0.01
# b = 2 / a
# get_acceptance_rate(a, b)

###################################
library(ggplot2)
library(dplyr)
library(purrr)

b <- seq(1, 5, length.out = 100)
as <- seq(1e-3, 0.25, length.out = 50)

# Create a dataframe
df <- map_df(as, ~{
  data.frame(
    a = .x,
    b = b,
    acceptance_rate = sapply(b/.x, get_acceptance_rate, a = .x)
  )
})

# Plot
ggplot(df, aes(x=b, y=acceptance_rate, color=a, group=a)) + 
  geom_line() +
  labs(x="a*c", y="Acceptance rate", color="a value") +
  theme_minimal() +
  scale_y_continuous(limits = c(0, 1)) +
  scale_color_gradient(low="#FCCF31", high="#F55555")




