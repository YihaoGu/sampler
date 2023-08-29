get_acceptance_rate = function(a, b) {
  target = function(x) {
    x / (1 + x ^ 2) * exp(-a ^ 2 * (x - b) ^ 2)
  }
  
  sd = sqrt(0.5) / a
  m = b
  
  # breaking integrals into more parts helps numerical stability
  f_area1 = integrate(target, lower = 0, upper = m)$value
  f_area2 = integrate(target, lower = m, upper = Inf)$value
  f_area = f_area1 + f_area2
  stopifnot(f_area > 1e-10) # warning for underflow
  
  g_area =  pnorm(0, mean = b, sd = sd, lower.tail = FALSE) * sqrt(pi) / a / 2
  
  return(f_area / g_area)
}

###################################
library(ggplot2)
library(dplyr)
library(purrr)

b <- seq(0.5, 1, length.out = 100)
as <- seq(1, 5, length.out = 50)

# Create a dataframe
df <- map_df(as, ~{
  data.frame(
    a = .x,
    b = b,
    acceptance_rate = sapply(b, get_acceptance_rate, a = .x)
  )
})

# Plot
ggplot(df, aes(x=b, y=acceptance_rate, color=a, group=a)) +
  geom_line() +
  labs(x="c", y="Acceptance rate", color="a value") +
  theme_minimal() +
  scale_y_continuous(limits = c(0, 1)) +
  scale_color_gradient(low="#FCCF31", high="#F55555")