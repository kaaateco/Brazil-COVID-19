# Title: SIR_divided_times
# Authors:
# Packages ----
library(deSolve)
library(outbreaks)
library(gridExtra)
library(arm)
library(tidyverse)
library(ggplot2)
periods <- function(date_initial, date_final,starting_par){
  # Set date ranges ----
  #date_initial = as.Date("2020-10-01")
  #date_final = as.Date("2021-04-01")
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
    dplyr::select(-`Province/State`, -Lat, -Long, -`Country/Region`)
  brazil_confirmed = tibble(date = colnames(brazil_confirmed),
                            cases_total = as.numeric(brazil_confirmed[1, ]))
  brazil_deaths = deaths %>%
    filter(`Country/Region` == "Brazil") %>% #, is.na(`Province/State`)
    dplyr::select(-`Province/State`, -Lat, -Long, -`Country/Region`)
  brazil_deaths = tibble(date = colnames(brazil_deaths),
                         deaths_total = as.numeric(brazil_deaths[1, ]))
  brazil_recovered = recoveries %>%
    filter(`Country/Region` == "Brazil") %>% #, is.na(`Province/State`)
    dplyr::select(-`Province/State`, -Lat, -Long, -`Country/Region`)
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
  sir_1 = function(beta, gamma, times, N, lambda, mu,I0,R0) {
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
                               # set the data set

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
                                                          # variables' intial values
                        times = times, N = N, lambda = lambda, mu = mu, I0=I0, R0=R0)    # time points
    # compute the sums of squares
    sum((predictions$I[-1] - data$I[-1])^2 + (predictions$R[-1] - data$R[-1])^2)
    #sum((predictions$I[-1] - data$I[-1])^2 )
  }
  # convenient wrapper to return sums of squares ----
  ss2 = function(x, N, data, lambda, mu) {
    ss(beta = x[1], gamma = x[2], N = N, data = data, lambda = lambda, mu = mu)
  }
  # set starting values ----
  starting_param_val = starting_par
  N = 212e6                                 # population size
  lambda = mu = 0                            # birth/death rate
  data = brazil   
  I0=data$I[1]
  R0=data$R[1]
  # set the data set
  # Optimization result ----
  ss_optim = optim(starting_param_val, ss2, N = N, data = data, lambda = lambda,
                   mu = mu)
  # Calculate R0 ----
  pars = ss_optim$par 
  # extract parameter estimates
  R = exp(pars[1]) / exp(pars[2])           # compute R0 for SIR
  # Predictions ----
  predictions = sir_1(beta = exp(pars[1]), gamma = exp(pars[2]), times = data$day, N = N, lambda = lambda,
                      mu = mu, I0 = I0, R0 = R0)              # generate predictions from the least
  # squares solution
  # Collect predictions into data frame ----
  pred_I = predictions$I
  pred_R = predictions$R
  date = seq(date_initial, date_final, by = 1)
  pred = data.frame(date, pred_I, pred_R)
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
    geom_line(mapping = aes(x = date, y = pred_I, color = colour),
              data = pred, size = 1,color=mn) +
    geom_bar(mapping = aes(x = date, y = I), stat = "identity",
             data = data, width = 0.5, fill = 'steelblue', alpha = 0.7,
    ) +
    xlim(date_initial, date_final)
  p1 = p1 + labs(y = "Active Cases")
  #ggsave("Cases_8months.pdf",p1,width=8, height=6)
  p2 = base +
    geom_line(mapping = aes(x = date, y = pred_R, color = colour),
              data = pred, size = 1,color = mn) +
    geom_bar(mapping = aes(x = date, y=R), stat="identity",
             data = data, width = 0.5, fill = 'steelblue', alpha = 0.7) +
    xlim(date_initial, date_final)
  p2 = p2 + labs(y = "Removed")
  #ggsave("Removed_8months.pdf",p2,width=8, height=6)
  p = grid.arrange(p1, p2)
  #ggsave("Both_6.png", p, width = 8, height = 6)
  return(list(pred,pars))
}


starting=log(c(1e-2,1e-5))
binder1 <- periods(as.Date("2020-11-01"), as.Date("2020-12-01"),starting)
binder2 <- periods(as.Date("2020-12-01"), as.Date("2021-01-01"),binder1[[2]])
binder3 <- periods(as.Date("2021-02-01"), as.Date("2021-03-01"),binder2[[2]])
binder4 <- periods(as.Date("2021-03-01"), as.Date("2021-04-01"),binder3[[2]])

bind <- rbind(binder1[[1]], binder2[[1]], binder3[[1]], binder4[[1]])
print(str(bind))
bind








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
  dplyr::select(-`Province/State`, -Lat, -Long, -`Country/Region`)
brazil_confirmed = tibble(date = colnames(brazil_confirmed),
                          cases_total = as.numeric(brazil_confirmed[1, ]))
brazil_deaths = deaths %>%
  filter(`Country/Region` == "Brazil") %>% #, is.na(`Province/State`)
  dplyr::select(-`Province/State`, -Lat, -Long, -`Country/Region`)
brazil_deaths = tibble(date = colnames(brazil_deaths),
                       deaths_total = as.numeric(brazil_deaths[1, ]))
brazil_recovered = recoveries %>%
  filter(`Country/Region` == "Brazil") %>% #, is.na(`Province/State`)
  dplyr::select(-`Province/State`, -Lat, -Long, -`Country/Region`)
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


brazil = brazil %>% 
  filter(date >= as.Date("2020-11-01"), date <= as.Date("2021-04-01"))



plot_periods = function(bind)
{ci = c("#C79999")
mn = c("#7C0000")
date_breaks = "1 month"
base = ggplot() +
  xlab("") +
  scale_x_date(
    date_breaks = date_breaks,
    labels = scales::date_format("%b %Y")
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 12)
  ) +
  theme(legend.position = "right")
p1 = base +
  geom_smooth(mapping = aes(x = date, y = pred_I, color = colour),
            data = bind, size = 1,color=mn, span = 0.2) +
  geom_bar(mapping = aes(x = date, y = I), stat = "identity",
           data = brazil, width = 0.5, fill = 'steelblue', alpha = 0.7,
  ) #+
#xlim(date_initial, date_final)
p1 = p1 + labs(y = "Active Cases")
#ggsave("Cases_8months.pdf",p1,width=8, height=6)
p2 = base +
  geom_smooth(mapping = aes(x = date, y = pred_R, color = colour),
            data = bind, size = 1,color = mn, span = 0.2) +
  geom_bar(mapping = aes(x = date, y= R), stat="identity",
           data = brazil, width = 0.5, fill = 'steelblue', alpha = 0.7) #+
#xlim(binder1, binder6)
p2 = p2 + labs(y = "Removed")
#("Removed_8months.pdf",p2,width=8, height=6)
p = grid.arrange(p1, p2)
#
grid.arrange(p1, p2)
}
plot_periods(bind)


#Need to add beta and gamma parameters