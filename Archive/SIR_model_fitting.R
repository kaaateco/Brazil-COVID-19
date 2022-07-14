ls_sir <- function(start, end, data){

# Title:

# Authors:

# Packages ----
library(deSolve)
library(outbreaks)
library(gridExtra)
library(arm)
library(tidyverse)
library(dplyr)

# data_read()

# Set date ranges ----
date_initial = start
date_final = end
print("input")
# Import data ----
brazil =
  brazil %>%
  mutate(total_removed = deaths_total + recoveries_total,
         active_cases = cases_total - total_removed,
         daily_removed = total_removed - lag(total_removed),
         daily_infected = cases_total - lag(cases_total),
         daily_deaths = deaths_total - lag(deaths_total),
         daily_recoveries = recoveries_total - lag(recoveries_total),
         day = 1:n(),
         I = active_cases,
         R = total_removed)
print("end")
brazil =
  brazil %>%
  filter(I < 2e7, R > -3e5)
print(str(date_initial))
brazil =
  brazil %>%
  filter(date >= as.Date(date_initial), date <= as.Date(date_final))

# SIR solver ----
sir_1 = function(beta, gamma, I0, R0, times, N, lambda, mu) {
  # define SIR equations
  sir_equations = function(time, variables, parameters) {
    with(as.list(c(variables, parameters)), {
      dS = -beta * I * S/N + lambda * N - mu * S
      dI =  beta * I * S/N - gamma * I -mu * I
      dR =  gamma * I - mu*R
      return(list(c(dS, dI, dR)))
    })
  }
  # prepare input for ODE solver
  parameters_values = c(beta = beta, gamma = gamma)
  S0 = N - I0 - R0
  initial_values = c(S = S0, I = I0, R = R0)
  # solve system of ODEs
  out = ode(initial_values, times, sir_equations, parameters_values,method = "rk4")
  return(as.data.frame(out))
}

starting_param_val = log(c(1e-2,1e-5))
N = 212000000                                 # population size
lambda = mu = 0                            # birth/death rate
data = brazil                               # set the data set

beta = 4
gamma = 0.1
I0 = 5
R0 = 2

times = data$day - 70
print(str(times))
predictions = sir_1(beta = beta, gamma = gamma,                        # parameters
                    I0 = I0, R0 = R0,                                  # variables' intial values
                    times = times, N = N, lambda = lambda, mu = mu)

predictions =
  predictions %>%
  pivot_longer(cols = c(S, I, R))

ggplot(predictions) +
   aes(x = time, y = value, group = name, color = name) +
  geom_line()

# Calculate sums of squares ----
ss = function(beta, gamma, N, data, lambda, mu) {
  # starting cases and removals on day 1
  I0 = data$I[1]
  R0 = data$R[1]
  times = data$day
  # transform parameters so they are non-negative
  beta = exp(beta)
  gamma = exp(gamma)
  # generate predictions using parameters, starting values
  predictions = sir_1(beta = beta, gamma = gamma,                        # parameters
                      I0 = I0, R0 = R0,                                  # variables' intial values
                      times = times, N = N, lambda = lambda, mu = mu)    # time points
  # compute the sums of squares
  sum((predictions$I[-1] - data$I[-1])^2 + (predictions$R[-1] - data$R[-1])^2)
  #sum((predictions$I[-1] - data$I[-1])^2 )
}

# convenient wrapper to return sums of squares ----
ss2 = function(x, N, data, lambda, mu) {
  ss(beta = x[1], gamma = x[2], N = N, data = data, lambda = lambda, mu = mu)
}

# set starting values ----
starting_param_val = log(c(1e-2,1e-5))
N = 1366e6                                 # population size
lambda = mu = 0                            # birth/death rate
data = brazil                               # set the data set

# Optimization result ----
ss_optim = optim(starting_param_val, ss2, N = N, data = data, lambda = lambda,
                 mu = mu)

# Calculate R0 ----
pars = ss_optim$par                       # extract parameter estimates
R = exp(pars[1]) / exp(pars[2])           # compute R0 for SIR

# Predictions ----
predictions = sir_1(beta = exp(pars[1]), gamma = exp(pars[2]), I0 = data$I[1],
                    R0 = data$R[1], times = data$day, N = N, lambda = lambda,
                    mu = mu)              # generate predictions from the least
                                          # squares solution

# Collect predictions into data frame ----
pred_I = predictions$I
pred_R = predictions$R
date = seq(date_initial, date_final, by = 1)
pred = data.frame(date, pred_I, pred_R)

}

#
# # Plot results ----
# ci = c("#C79999")
# mn = c("#7C0000")
# date_breaks = "1 month"
#
# base = ggplot() +
#   xlab("") +
#   scale_x_date(
#     date_breaks = date_breaks,
#     labels = scales::date_format("%e %b")
#   ) +
#   theme_bw() +
#   theme(
#     axis.text.x = element_text(angle = 45, hjust = 1),
#     axis.text = element_text(size = 12),
#     axis.title = element_text(size = 12)
#   ) +
#   theme(legend.position = "right")
#
# p1 = base +
#   geom_line(mapping = aes(x = date, y = pred_I, color = colour),
#             data = pred, size = 1,color=mn) +
#   geom_bar(mapping = aes(x = date, y = I), stat = "identity",
#            data = data, width = 0.5, fill = 'steelblue', alpha = 0.7,
#   ) + xlim(date_initial, date_final)
#
# p1 = p1 + labs(y = "Active Cases")
# #ggsave("Cases_8months.pdf",p1,width=8, height=6)
#
# p2 = base +
#   geom_line(mapping = aes(x = date, y = pred_R, color = colour),
#     data = pred, size = 1,color = mn) +
#   geom_bar(mapping = aes(x = date, y=R), stat="identity",
#     data = data, width = 0.5, fill = 'steelblue', alpha = 0.7) +
#   xlim(date_initial, date_final)
#
# p2 = p2 + labs(y = "Removed")
# #ggsave("Removed_8months.pdf",p2,width=8, height=6)
#
# p = grid.arrange(p1, p2)
# #ggsave("Both_6.png", p, width = 8, height = 6)
#
# p
