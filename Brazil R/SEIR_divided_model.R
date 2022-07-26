# Title: SEIR_Dividing_Times

# Authors:

# Packages ----
library(deSolve)
library(outbreaks)
library(gridExtra)
library(arm)
library(tidyverse)
library(dplyr)

periods_SEIR <- function(date_initial, date_final,starting_par){

# SIR solver ----
seir_1 = function(beta, gamma, I0, R0, E0, times, N, lambda, mu,sigma, K) {
  # define SIR equations
  seir_equations = function(time, variables, parameters) {
    with(as.list(c(variables, parameters)), {
      dS = -beta * I * S/N + lambda * N - mu * S
      dE = (beta * I * S)/N - (sigma * E) - (mu * E)
      dI =  sigma * E - gamma * I -mu * I
      dR =  gamma * I - mu*R
      return(list(c(dS, dE, dI, dR)))
    })
  }



  # prepare input for ODE solver
  parameters_values = c(beta = beta, gamma = gamma)
  E0 = K*I0
  S0 = N - I0 - R0- E0
  initial_values = c(S = S0,E = E0, I = I0, R = R0)

  # solve system of ODEs
  out = ode(initial_values, times, seir_equations, parameters_values,method = "rk4")
  return(as.data.frame(out))
}



# N = 212000000                         # population size
# lambda = mu = 0                       # birth/death rate
# data = brazil                         # set the data set
#
# beta = 4
# gamma = 0.1
# sigma = 1/5.6
# I0 = 5
# R0 = 2
#
# times = data$day - 70
#
# predictions = seir_1(beta = beta, gamma = gamma,                     # parameters
#                     I0 = I0, R0 = R0, E0 = E0,                       # variables' intial values
#                     times = times, N = N, lambda = lambda, mu = mu)
#
# predictions =
#   predictions %>%
#   pivot_longer(cols = c(S, E, I, R))
#
# ggplot(predictions) +
#    aes(x = time, y = value, group = name, color = name) +
#   geom_line()


# Calculate sums of squares ----
ss = function(beta, gamma, N, data, lambda, mu,sigma,K) {
  # starting cases and removals on day 1
  I0 = data$I[1]
  R0 = data$R[1]
  times = data$day

  # transform parameters so they are non-negative
  beta = exp(beta)
  gamma = exp(gamma)

  # generate predictions using parameters, starting values
  predictions = seir_1(beta = beta, gamma = gamma,              # parameters
                       I0 = I0, R0 = R0,                         # variables' intial values
                       times = times, N = N, lambda = lambda, mu = mu, sigma = sigma,K=K)    # time points
  # compute the sums of squares
  sum((predictions$I[-1] - data$I[-1])^2 + (predictions$R[-1] - data$R[-1])^2)
}


# convenient wrapper to return sums of squares ----
ss2 = function(x, N, data, lambda, mu, sigma,K) {
  ss(beta = x[1], gamma = x[2], N = N, data = data, lambda = lambda, mu = mu, sigma=sigma,K=K)
}


# set starting values ----
starting_param_val = starting_par
N = 212000000                             # population size
lambda = mu = 1/(75.88*365)               # birth/death rate
data = brazil                             # set the data set
sigma=1/5.6
I0=data$I[1]
R0=data$R[1]

# Optimization result ----


sse<-function(K){
  ss_optim = optim(starting_param_val, ss2, N = N, data = data, lambda = lambda,
                   mu = mu, sigma=sigma, K=K)
  return(ss_optim$value)
}
# Calculate R0 ----


k_seq=seq(0,3,length.out=30)
ssevec=rep(0,times=30)
for (i in 1:length(k_seq)){
  ssevec[i]=sse(k_seq[i])
}
K=k_seq[which.min(ssevec)]
ss_optim = optim(starting_param_val, ss2, N = N, data = data, lambda = lambda,
                 mu = mu, sigma=sigma,K=K)
pars = ss_optim$par                       # extract parameter estimates
R = exp(pars[1]) / exp(pars[2])           # compute R0 for SIR


# Predictions ----
predictions = seir_1(beta = exp(pars[1]), gamma = exp(pars[2]), I0 = I0,
                     R0 = R0, times = data$day, N = N, lambda = lambda,
                     mu = mu, sigma=sigma,K=K)
# generate predictions from the least squares solution


# Collect predictions into data frame ----
pred_E = predictions$E
pred_I = predictions$I
pred_R = predictions$R
date = seq(date_initial, date_final, by = 1)
pred = data.frame(date, pred_E, pred_I, pred_R)


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
  ) + xlim(date_initial, date_final)

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

p

return(list(pred,pars,K))
}


sse((starting_param_val, ss2, N = N, data = data, lambda = lambda,
     mu = mu, sigma=sigma,K=K))

starting_par=log(c(1e-2,1e-5))
binder1 <- periods_SEIR(as.Date("2020-11-01"), as.Date("2020-12-01"),starting_par)
binder2 <- periods_SEIR(as.Date("2020-12-01"), as.Date("2021-01-01"), c(binder1[[2]]: binder1[[3]]))
binder3 <- periods_SEIR(as.Date("2021-02-01"), as.Date("2021-03-01"),c(binder2[[2]]: binder2[[3]]))
binder4 <- periods_SEIR(as.Date("2021-03-01"), as.Date("2021-04-01"),c(binder3[[2]]: binder3[[3]]))

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



plot_periods = function(bind, brazil)
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
              data = bind, size = 1,color = "mediumblue", span = 0.2) +
  geom_bar(mapping = aes(x = date, y = I), stat = "identity",
           data = brazil, width = 0.5, fill = 'green4', alpha = 0.7,
  )

p1 = p1 + labs(y = "Active Cases")
#ggsave("Cases_8months.pdf",p1,width=8, height=6)

p2 = base +
  geom_smooth(mapping = aes(x = date, y = pred_R, color = colour),
              data = bind, size = 1, color = "mediumblue", span = 0.2) +
  geom_bar(mapping = aes(x = date, y= R), stat="identity",
           data = brazil, width = 0.5, fill = 'green4', alpha = 0.7) #+

p2 = p2 + labs(y = "Removed")

#p = grid.arrange(p1, p2)
#grid.arrange(p1, p2)

ecase_plot <- ggplotly(p1)
edeath_plot <- ggplotly(p2)

subplot(ecase_plot, edeath_plot, titleX = T, titleY = T, margin = .08,
        nrows = 2, shareX = F) %>%
  plotly::config(toImageButtonOptions = list(width = NULL))

}

plot_periods(bind, brazil)


#Need to add beta and gamma parameters
