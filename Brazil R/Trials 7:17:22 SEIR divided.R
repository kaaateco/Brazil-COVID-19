# Title: SEIR_Dividing_Times

# Authors:

# Packages ----
library(deSolve)
library(outbreaks)
library(gridExtra)
library(arm)
library(tidyverse)
library(dplyr)

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

periods_SEIR <- function(date_initial, date_final,starting_par){
  
  
  
  # Set date ranges ----
  #date_initial = as.Date("2020-11-01")
  #date_final = as.Date("2020-12-01")
  
  
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
    ss(beta = x[1], gamma = x[2], N = N, data = brazil, lambda = lambda, mu = mu, sigma=sigma,K=K)
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
  
  
  sse<-function(starting_param_val, ss2, N = N, data = brazil, lambda = lambda,
                mu = mu, sigma=sigma, K = K){
    ss_optim = optim(starting_param_val, ss2, N = N, data = brazil, lambda = lambda,
                     mu = mu, sigma=sigma,K=K)
    return(ss_optim$value)
  }
  # Calculate R0 ----
  
  
  k_seq=seq(0,3,length.out=30)
  ssevec=rep(0,times=30)
  for (i in 1:length(k_seq)){
    ssevec[i]=sse(starting_param_val, ss2, N = N, data = brazil, lambda = lambda,
                  mu = mu, sigma=sigma, K = k_seq[i])
  }
  
  K=k_seq[which.min(ssevec)]
  ss_optim = optim(starting_param_val, ss2, N = N, data = brazil, lambda = lambda,
                   mu = mu, sigma=sigma,K=K)
  pars = ss_optim$par                       # extract parameter estimates
  R = exp(pars[1]) / exp(pars[2])           # compute R0 for SIR
  
  
  # Predictions ----
  predictions = seir_1(beta = exp(pars[1]), gamma = exp(pars[2]), I0 = I0,
                       R0 = R0, times = data$day, N = N, lambda = lambda,
                       mu = mu, sigma=sigma,K=K)
  # generate predictions from the least squares solution
  pred=predictions[,c("I","R")]
  
  
  #PREDICTIONS!!!!!!!!!!
  
  ciband<-function(pred,sigma_l,sigma_u,sigma_m,beta,gamma,data,rep,K){
    pred_I_med=pred$I
    pred_R_med=pred$R
    sd=abs(1/sigma_l-1/sigma_u)/(2*1.96)
    pred_I=matrix(0,nrow=nrow(pred),ncol=rep)
    pred_R=matrix(0,nrow=nrow(pred),ncol=rep)
    date=data$date
    
    
    
    for(i in 1:rep){
      De=rnorm(1,mean=1/sigma_m, sd=sd)
      sigma=1/De
      predictions = seir_1(beta = beta, gamma = gamma, I0 = data$I[1],
                           R0 = data$R[1], times = data$day, N = N, lambda = lambda,
                           mu = mu, sigma=sigma, K = K)
      p_I=rpois(nrow(predictions),lambda=predictions$I)
      p_R=rpois(nrow(predictions),predictions$R)
      pred_I[,i]=p_I
      pred_R[,i]=p_R
    }
    lwrI = round(apply(pred_I,1,quantile, probs=0.025))
    uprI = round(apply(pred_I,1,quantile, probs=0.975))
    sd=(uprI-lwrI)/(1.96*2)
    lwrI=pred_I_med-1.96*sd
    uprI=pred_I_med+1.96*sd
    pred_I=data.frame(date,pred_I_med,lwrI,uprI)
    lwrR = round(apply(pred_R,1,quantile, probs=0.025))
    uprR = round(apply(pred_R,1,quantile, probs=0.975))
    sd=(uprR-lwrR)/(1.96*2)
    lwrR=pred_R_med-1.96*sd
    uprR=pred_R_med+1.96*sd
    pred_R=data.frame(date,pred_R_med,lwrR,uprR)
    return(list(pred_I,pred_R))
  }
  
  
  
  
  
  pred_ci<- ciband(pred = pred, sigma_l = 1/6.3, sigma_u = 1/5.0, sigma_m = 
                     1/5.6, beta = exp(pars[1]), gamma = exp(pars[2]), data = brazil, rep = 100, K= K)
  
  
  
  
  
  # Collect predictions into data frame ----
  
  
  
  pred_I = pred_ci[[1]]
  pred_R= pred_ci[[2]]
  
  
  
  
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
  
  p
  
  return(list(pred_I,pred_R,pars,K))
}

#print(str(data))
#sse(starting_param_val, ss2 = ss2, N = N, data = brazil, lambda = lambda,
    #mu = mu, sigma=sigma, K=K)


starting_par=log(c(1e-2,1e-5))
binder1 <- periods_SEIR(as.Date("2020-11-01"), as.Date("2020-11-30"),starting_par)
binder2 <- periods_SEIR(as.Date("2020-12-01"), as.Date("2020-12-31"), c(binder1[[3]][1], binder1[[3]][2]))
binder3 <- periods_SEIR(as.Date("2021-01-01"), as.Date("2021-01-31"), c(binder2[[3]][1], binder2[[3]][2]))
binder4 <- periods_SEIR(as.Date("2021-02-01"), as.Date("2021-02-28"),c(binder3[[3]][1], binder3[[3]][2]))
binder5 <- periods_SEIR(as.Date("2021-03-01"), as.Date("2021-03-31"),c(binder4[[3]][1], binder4[[3]][2]))

##Testing Period

I0 <- 1285159
R0 <- 11566006
N = 212000000
lambda = mu = 1/(75.88*365)
sigma=1/5.6
sigma_l = 1/6.3 
sigma_u = 1/5.0 
sigma_m = 1/5.6

testing_seir <- function(date_initial,date_final, beta, gamma, K) {
date=seq(date_initial,date_final,by=1)
test_length=length(seq(date_initial,date_final,by=1))
predictions_test = seir_1(beta = beta, gamma = gamma, I0 = I0,
                     R0 = R0, times = c(1:test_length), N = N, lambda = lambda,
                     mu = mu, sigma=sigma,K=K)


# generate predictions from the least squares solution
pred=predictions_test[,c("I","R")]


#PREDICTIONS!!!!!!!!!!
rep=100
ciband<-function(pred,sigma_l,sigma_u,sigma_m,beta,gamma,rep,K){
  pred_I_med=pred$I
  pred_R_med=pred$R
  sd=abs(1/sigma_l-1/sigma_u)/(2*1.96)
  pred_I=matrix(0,nrow=nrow(pred),ncol=rep)
  pred_R=matrix(0,nrow=nrow(pred),ncol=rep)
  
  
  
  for(i in 1:rep){
    De=rnorm(1,mean=1/sigma_m, sd=sd)
    sigma=1/De
    predictions_test = seir_1(beta = beta, gamma = gamma, I0 = I0,
                         R0 = R0, times = c(1:test_length), N = N, lambda = lambda,
                         mu = mu, sigma=sigma, K = K)
    p_I=rpois(nrow(predictions_test),lambda=predictions_test$I)
    p_R=rpois(nrow(predictions_test),predictions_test$R)
    pred_I[,i]=p_I
    pred_R[,i]=p_R
  }
  lwrI = round(apply(pred_I,1,quantile, probs=0.025))
  uprI = round(apply(pred_I,1,quantile, probs=0.975))
  sd=(uprI-lwrI)/(1.96*2)
  lwrI=pred_I_med-1.96*sd
  uprI=pred_I_med+1.96*sd
  pred_I=data.frame(date,pred_I_med,lwrI,uprI)
  lwrR = round(apply(pred_R,1,quantile, probs=0.025))
  uprR = round(apply(pred_R,1,quantile, probs=0.975))
  sd=(uprR-lwrR)/(1.96*2)
  lwrR=pred_R_med-1.96*sd
  uprR=pred_R_med+1.96*sd
  pred_R=data.frame(date,pred_R_med,lwrR,uprR)
  return(list(pred_I,pred_R))
}
testing_ci <- ciband(pred = pred, sigma_l = 1/6.3, sigma_u = 1/5.0, sigma_m = 
                        1/5.6, beta = beta, gamma = gamma, rep = rep, K= K)
return(testing_ci)
}



binder_test_period <- testing_seir(date_initial=as.Date("2021-04-01"), date_final=as.Date("2021-04-30"), 
                                   beta = exp(binder4[[3]][1]), gamma = exp(binder4[[3]][2]), K = binder4[[4]])
binder_test_period

pred_bind_I <- rbind(binder1[[1]], binder2[[1]], binder3[[1]], binder4[[1]], binder5[[1]], binder_test_period[[1]])
pred_bind_R <- rbind(binder1[[2]], binder2[[2]], binder3[[2]], binder4[[2]], binder5[[2]], binder_test_period[[2]])


bind_I <- rbind(binder1[[1]], binder2[[1]], binder3[[1]], binder4[[1]], binder5[[1]])
bind_R <- rbind(binder1[[2]], binder2[[2]], binder3[[2]], binder4[[2]], binder5[[2]])

bind_I
bind_R

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
         I = active_cases,
         R = total_removed)
brazil =
  brazil %>%
  filter(I < 2e7, R > -3e5)



brazil = brazil %>% 
  filter(date >= as.Date("2020-11-01"), date <= as.Date("2021-04-30"))
brazil$day=c(1:nrow(brazil))

brazil = brazil %>% 
  mutate(train_test = ifelse(as.Date(date) >= "2021-04-01", "test", "train"))

pred_bind_I = pred_bind_I %>% 
  mutate(train_test = ifelse(as.Date(date) >= "2021-04-01", "test", "train"))

pred_bind_R = pred_bind_R %>% 
  mutate(train_test = ifelse(as.Date(date) >= "2021-04-01", "test", "train"))


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
  theme(legend.position = "none") 
p1 = base +
  geom_smooth(mapping = aes(x = date, y = pred_I_med, color = train_test, fill = train_test),
              data = pred_bind_I, size = 1, span = 0.2) +
  ggplot2::geom_ribbon(
    mapping = ggplot2::aes(x = date, ymin = lwrI, ymax=uprI, color = train_test, fill = train_test),
    data = pred_bind_I,
    size = 1,fill=ci,alpha=0.8
  ) +
  geom_bar(mapping = aes(x = date, y = I, color = train_test, fill = train_test), stat = "identity",
           data = brazil, width = 0.1, fill = 'steelblue', alpha = 0.7, size = 0.1) +
  geom_vline(xintercept = as.Date("2021-04-01"))
#xlim(date_initial, date_final)
p1 = p1 + labs(y = "Active Cases")
#ggsave("Cases_8months.pdf",p1,width=8, height=6)
p2 = base +
  geom_smooth(mapping = aes(x = date, y = pred_R_med, color = train_test, fill = train_test),
              data = pred_bind_R, size = 1,color = mn, span = 0.2) +
  ggplot2::geom_ribbon(
    mapping = ggplot2::aes(x = date, ymin = lwrR, ymax=uprR, color = train_test, fill = train_test),
    data = pred_bind_R,
    size = 1,fill=ci,alpha=0.8,
  ) +
  geom_bar(mapping = aes(x = date, y= R, color = train_test, fill = train_test), stat="identity",
           data = brazil, width = 0.1, fill = 'steelblue', alpha = 0.7, size = 0.1) +
  geom_vline(xintercept = as.Date("2021-04-01")) 
#xlim(binder1, binder6)
p2 = p2 + labs(y = "Removed")
#("Removed_8months.pdf",p2,width=8, height=6)
p = grid.arrange(p1, p2)
#
grid.arrange(p1, p2)
}
plot_periods(pred_bind_I)

#SMAPE for testing period (I)

T_test = 30
SMAPE_SEIR_test <- function(model){
  brazil_smape = brazil %>% 
    filter(train_test == "test")
  pred_bind_I_smape = pred_bind_I %>% 
    filter(train_test == "test")
  abs_num <- abs(pred_bind_I_smape$pred_I_med - brazil_smape$I)
  abs_denom <- (abs(brazil_smape$I) + abs(pred_bind_I_smape$pred_I_med))/2
  summation <- sum(abs_num/abs_denom)
  return ((100/T_test)*(summation))
}

SMAPE_SEIR_test(bind_I)

T_train = 151
SMAPE_SEIR_train <- function(model){
  brazil_smape = brazil %>% 
    filter(train_test == "train")
  abs_num <- abs(bind_I$pred_I_med - brazil_smape$I)
  abs_denom <- (abs(brazil_smape$I) + abs(bind_I$pred_I_med))/2
  summation <- sum(abs_num/abs_denom)
  return ((100/T_train)*(summation))
}

SMAPE_SEIR_train(bind_I)





