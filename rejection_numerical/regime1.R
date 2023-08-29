get_acceptance_rate = function(a, b) {
  target = function(x) {
    x / (1 + x ^ 2) * exp(-a ^ 2 * (x - b) ^ 2)
  }
  
  m = (1 + sqrt(1 - a^2)) / a
  mu = b
  std = 1 / (a * sqrt(2))
  
  # breaking integrals into more parts helps numerical stability,
  # also makes it easy to check which region contributes most mass
  f_area1 = integrate(target, lower = 0, upper = m)$value
  f_area2 = integrate(target, lower = m, upper = Inf)$value
  f_area = f_area1 + f_area2
  stopifnot(f_area > 1e-10) # warning for underflow
  
  g_area1 = 0.5 * log(1 + m^2)
  g_area2 = (1 -  pnorm(m, mean = mu, sd = std)) * (sqrt(pi)) / 2
  g_area =  g_area1 + g_area2

  # cat('Mass:', g_area1, g_area2, '\n')
  # cat('Acceptance rate:', f_area1 / g_area1, f_area2 / g_area2, '\n')
  return(f_area / g_area)
}

###################################
library(ggplot2)
library(dplyr)
library(purrr)

b <- seq(0, 1, length.out = 100)
as <- seq(1e-3, 1, length.out = 50)

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

