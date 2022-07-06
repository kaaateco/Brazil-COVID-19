# Title:

# Authors:

# Packages ----
library(deSolve)
library(outbreaks)
library(gridExtra)
library(arm)
library(tidyverse)
library(bbmle)
library(dplyr)

# Set date ranges ----
date_initial = as.Date("2020-04-01")
date_final = as.Date("2020-06-01")

# Import data ----
# confirmed
case_url = "https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_confirmed_global.csv"
confirmed = read_csv(case_url)

# deaths
death_url = "https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_deaths_global.csv"
deaths = read_csv(death_url)

# recoveries
recovery_url = "https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_recovered_global.csv"
recoveries = read_csv(recovery_url)

# Process and merge data ----
brazil_confirmed = confirmed %>%
  filter(`Country/Region` == "Brazil") %>% #, is.na(`Province/State`)
  select(-`Province/State`, -Lat, -Long, -`Country/Region`)
brazil_confirmed = tibble(date = colnames(brazil_confirmed),
                           cases_total = as.numeric(brazil_confirmed[1, ]))

brazil_deaths = deaths %>%
  filter(`Country/Region` == "Brazil") %>% #, is.na(`Province/State`)
  select(-`Province/State`, -Lat, -Long, -`Country/Region`)
brazil_deaths = tibble(date = colnames(brazil_deaths),
                        deaths_total = as.numeric(brazil_deaths[1, ]))

brazil_recovered = recoveries %>%
  filter(`Country/Region` == "Brazil") %>% #, is.na(`Province/State`)
  select(-`Province/State`, -Lat, -Long, -`Country/Region`)
brazil_recovered = tibble(date = colnames(brazil_recovered),
                           recoveries_total = as.numeric(brazil_recovered[1, ]))

brazil = brazil_confirmed %>% left_join(brazil_recovered, by = "date")
brazil = brazil %>% left_join(brazil_deaths, by = "date")

brazil =
  brazil %>%
  mutate(date = as.Date(date, format = "%m/%d/%y"))

brazil =
  brazil %>%
  mutate(total_removed = deaths_total + recoveries_total,
         active_cases = cases_total - total_removed,
         day = 1:n(),
         I = active_cases,
         R = total_removed)

brazil =
  brazil %>%
  filter(I < 2e7, R > -3e5)

brazil =
  brazil %>%
  filter(date >= date_initial, date <= date_final)

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
N = 1366e6                                 # population size
lambda = mu = 0                            # birth/death rate
data = brazil                               # set the data set

beta = 2
gamma = 0.05
I0 = 5
R0 = 2

times = data$day - 70

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

##loglikelihood function for mle
logli = function(beta, gamma, N, dat, lambda, mu) {
  I0 = dat$I[1]
  R0 = dat$R[1]
  times = dat$day
  beta = exp(beta)
  gamma = exp(gamma)
  predictions = sir_1(beta = beta, gamma = gamma,   # parameters
                       I0 = I0, R0 = R0, # variables' intial values
                       times = times, N = N, lambda = lambda, mu = mu)
  ## negative of log likelihood
  -sum(dpois(x = dat$I, lambda = predictions$I, log = TRUE)) - sum(dpois(x = dat$R, lambda = predictions$R, log = TRUE))
}


method = "mle"

if(method == "ls"){
  # set starting values ----
  starting_param_val = log(c(1e-2,1e-5))
  N = 1366e6                                 # population size
  lambda = mu = 0                            # birth/death rate
  data = brazil                               # set the data set

  # Optimization result ----
  ss_optim = optim(starting_param_val, ss2, N = N, data = data, lambda = lambda,
                   mu = mu)

  # Obtain beta, gamma
  pars = ss_optim$par
}
if(method == "mle"){
  N=1366e6
  lambda=mu=0

  starting_param_val = list(beta = 1e-2, gamma = 1e-5)
  estimates_pois = mle2(minuslogl = logli,
                         start = lapply(starting_param_val, log), method = "BFGS",
                         data=list(dat = data, N = N, lambda = lambda, mu = mu))
  pars = as.numeric(coef(estimates_pois))
}




cbind(exp(pars), exp(confint(estimates_pois)))

R = as.numeric(exp(pars[1]) / exp(pars[2]))           # compute R0 for SIR

# Predictions ----
predictions = sir_1(beta = exp(pars[1]), gamma = exp(pars[2]), I0 = data$I[1],
                    R0 = data$R[1], times = data$day, N = N, lambda = lambda,
                    mu = mu)              # generate predictions from the least
                                          # squares solution

# Collect predictions into data frame ----
date = seq(date_initial, date_final, by = 1)
pred_I_med = round(predictions$I)
pred_R_med = round(predictions$R)

## 95% C.I using theoretical quantiles of the poisson distribution
cl = 0.95
cl = (1 - cl) / 2
lwrI = qpois(p = cl, lambda = pred_I_med)
uprI = qpois(p = 1 - cl, lambda = pred_I_med)
pred_I=data.frame(date,pred_I_med,lwrI,uprI)

lwrR = qpois(p = cl, lambda = pred_R_med)
uprR = qpois(p = 1 - cl, lambda = pred_R_med)
pred_R=data.frame(date,pred_R_med,lwrR,uprR)

# Plot results ----
ci = c("#C79999")
mn = c("#7C0000")
date_breaks = "1 month"

base = ggplot() +
  xlab("") +
  scale_x_date(
    date_breaks = date_breaks,
    labels = scales::date_format("%e %b")
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 12)
  ) +
  theme(legend.position = "right")

p1 = base +
  geom_line(mapping = aes(x = date, y = pred_I_med, color = colour),
            data = pred_I, size = 0.5, color = mn) +
  geom_ribbon(
    mapping = aes(x = date, ymin = lwrI, ymax = uprI),
    data = pred_I,
    size = 1, fill = ci, alpha = 0.8,
  ) +
  geom_bar(mapping = aes(x = date, y = I), stat = "identity",
           data = data, width = 0.5, fill = 'steelblue', alpha = 0.7,
  ) +
  xlim(date_initial, date_final)

p1 = p1 + labs(y = "Active Cases")
#ggsave("Cases_8months.pdf",p1,width=8, height=6)

p2 = base +
  geom_line(mapping = aes(x = date, y = pred_R_med, color = colour),
            data = pred_R, size = 1,color=mn) +
  ggplot2::geom_ribbon(
    mapping = ggplot2::aes(x = date, ymin = lwrR, ymax=uprR),
    data =pred_R,
    size = 1,fill=ci,alpha=0.8,
  )+
  geom_bar(mapping = aes(x = date, y = R), stat = "identity",
           data = data, width = 0.5, fill = 'steelblue', alpha = 0.7,
  ) +
  xlim(date_initial, date_final)
p2 = p2 + labs(y = "Removed")
#ggsave("Removed_8months.pdf",p2,width=8, height=6)

p = grid.arrange(p1, p2)
#ggsave("Both_6.png", p, width = 8, height = 6)

