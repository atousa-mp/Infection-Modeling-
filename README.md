# SEIR Model for Covid-19 for the population of Sweden, Susceptible(S) = 10,499,990, Exposed (E) = 5, Infectious (I) = 5, Recovered (R) = 0 
# β (transmission rate) = 0.5, σ (incubation rate) = 1/5.2 days , γ (recovery rate) = 1/10 days
# Install required packages if not installed
if(!require(deSolve)) install.packages("deSolve")
if(!require(ggplot2)) install.packages("ggplot2")
if(!require(tidyr)) install.packages("tidyr")

library(deSolve)
library(ggplot2)
library(tidyr)

# Total population of Sweden
N <- 10500000

# Initial conditions
init <- c(
  S = N - 10,
  E = 5,
  I = 5,
  R = 0
)

# Parameters
parameters <- c(
  beta = 0.5,     # transmission rate
  sigma = 1/5.2,  # incubation rate
  gamma = 1/10    # recovery rate
)

# Time (days)
times <- seq(0, 200, by = 1)

# SEIR model function
seir_model <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {

    dS <- -beta * S * I / N
    dE <- beta * S * I / N - sigma * E
    dI <- sigma * E - gamma * I
    dR <- gamma * I

    return(list(c(dS, dE, dI, dR)))
  })
}

# Solve the system
out <- ode(y = init, times = times, func = seir_model, parms = parameters)
out <- as.data.frame(out)

# Convert to long format for plotting
out_long <- pivot_longer(out, cols = c("S","E","I","R"),
                         names_to = "Compartment",
                         values_to = "Population")

# Plot epidemic curves
ggplot(out_long, aes(x = time, y = Population, color = Compartment)) +
  geom_line(size = 1.2) +
  labs(
    title = "SEIR Model Simulation of COVID-19 in Sweden",
    x = "Days",
    y = "Population"
  ) +
  theme_minimal()
