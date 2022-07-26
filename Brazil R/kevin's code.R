source(“israel_functions.R”)
library(deSolve)
library(outbreaks)
library(gridExtra)
library(arm)
library(tidyverse)
library(bbmle)
library(zoo)
library(eSIR)
library(lubridate)
library(dplyr)
library(ggplot2)
get_R_data = function(){
  # confirmed
  case_url = “https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_confirmed_global.csv”
  confirmed = read_csv(case_url)
  # deaths
  death_url = “https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_deaths_global.csv”
  deaths = read_csv(death_url)
  # recoveries
  recovery_url = “https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_recovered_global.csv”
  recoveries = read_csv(recovery_url)
  # Next steps: Inspect the data with
  # Process the data with dplyr (or software of your choice)
  # Data format should be
  #      - each row is a day (with date)
  #      - the columns are daily cases, deaths, recoveries.
  # -----
  # Multiple regions under this country?
  confirmed %>% filter(`Country/Region` == “Israel”) %>%
    pull(`Province/State`) %>% unique()
  israel_confirmed = confirmed %>%
    filter(`Country/Region` == “Israel”) %>% #, is.na(`Province/State`)
    dplyr::select(-`Province/State`, -Lat, -Long, -`Country/Region`)
  israel_confirmed = tibble(date = colnames(israel_confirmed),
                            cases_total = as.numeric(israel_confirmed[1, ]))
  # process/merge deaths and recoveries the in the same fashion
  # -----
  israel_deaths = deaths %>%
    filter(`Country/Region` == “Israel”) %>% #, is.na(`Province/State`)
    dplyr::select(-`Province/State`, -Lat, -Long, -`Country/Region`)
  israel_deaths = tibble(date = colnames(israel_deaths),
                         deaths_total = as.numeric(israel_deaths[1, ]))
  israel_recovered = recoveries %>%
    filter(`Country/Region` == “Israel”) %>% #, is.na(`Province/State`)
    dplyr::select(-`Province/State`, -Lat, -Long, -`Country/Region`)
  israel_recovered = tibble(date = colnames(israel_recovered),
                            recoveries_total = as.numeric(israel_recovered[1, ]))
  israel = israel_confirmed %>% left_join(israel_recovered, by = “date”)
  israel = israel %>% left_join(israel_deaths, by = “date”)
  # Transform data for plotting
  as.Date(“1/24/20”, format = “%m/%d/%y”)
  israel =
    israel %>%
    mutate(date = as.Date(date, format = “%m/%d/%y”))
  # daily case, death, recoveries
  israel =
    israel %>%
    mutate(total_removed = deaths_total + recoveries_total,
           total_active_infected = cases_total - total_removed,
           day_index = 1:n(),
           daily_removed = total_removed - lag(total_removed),
           daily_infected = cases_total - lag(cases_total),
           I = daily_infected,
           R = daily_removed) %>% filter (R > -3e4) %>% filter (I >= 0)
  # Static plots
  ggplot(israel) +
    aes(x = date, y = I) +
    geom_col()
  israel =
    israel %>%
    filter(date > as.Date(“2020-12-27"),
           date < as.Date(“2021-04-18”))
  # set the mean serial interval
  serial_interval = 4.0
  serial_sd = (4.9 - 4.0) / 1.96
  # 95% interval = 3.1 to 4.9
  # for a five day window
  t_start = seq(2, nrow(israel) - 4)
  t_end = t_start + 4
  res <- EpiEstim::estimate_R(
    incid = israel$daily_infected,
    method = “parametric_si”,
    config = EpiEstim::make_config(list(
      mean_si             = serial_interval,
      std_si              = serial_sd,
      si_parametric_distr = “G”,
      t_start             = t_start,
      t_end               = t_end,
      seed                = 46342))
  )
  return (res)
}
#source: https://ourworldindata.org/covid-vaccinations?country=ISR
#x <- “11/14/2006 20:10"
get_eSIR_plots = function(){
  date0 <- as.Date(“2020-12-27”)
  date1 <- as.Date(“2021-01-07")
date2 <- as.Date(“2021-02-07”)
date3 <- as.Date(“2021-03-07")
  date4 <- as.Date(“2021-04-18”)
  date4_half <-  as.Date(“2021-04-10")
test_date = as.Date(“2021-05-24”)
date_initial = date0
date_final = test_date
israel <- israel_data_only(date_initial, date_final)
israel = israel %>%
  mutate(period = ifelse(date < date1, 1,
                         ifelse(date >= date1 & date < date2, 2,
                                ifelse(date >= date2 & date < date3, 3,
                                       ifelse(date >= date3 & date<date4, 4, ifelse(date>=date4, 5, 0))))))
#test_df = israel%>% filter(date >= as.Date(“2021-04-10"))
test_df = israel%>% filter(period==5)
train_df = israel %>% filter(period <5)
NI_complete = train_df$cases_total
RI_complete = train_df$total_removed
N=9449000
R <- RI_complete / N
Y <- NI_complete / N - R
change_time <- c(“01/07/2021”, “02/07/2021", “03/07/2021”)
pi0 = get_R_data()$R %>% mutate(death_07da = zoo::rollmeanr(`Mean(R)`, k = 7, fill = NA))
pi0 = get_R_data()$R
pi0 = pi0 %>% rename(Mean_R  = `Mean(R)`)
first_seven = mean(pi0$`Mean_R`[1:7])
pi0 = pi0 %>% mutate(pi = Mean_R/first_seven)
R = R[1:nrow(pi0)]
Y = Y[1:nrow(pi0)]
pi0$date = as.Date(date0)
for(i in 1:nrow(pi0)){
  pi0$date[i] = israel$date[i]
}
pi0 = pi0 %>% mutate(date = format(date, format = “%m/%d/%Y”))
pi0$Y = Y
pi0$R = R
change_time = pi0$date
change_time = change_time[2:length(change_time)-1]
#change_time = c(“01/07/2021”, “02/07/2021", “03/07/2021”)
#pi = c(1.25, 1, 1.5, 0.5)
vac = read_csv(“vaccination.csv”) %>% filter(date >= date_initial & date<=date4_half) %>% mutate(vac_inv = 1-people_vaccinated/N)
pi = vac$vac_inv
#model = tvt.eSIR(Y = pi0$Y, R = pi0$R, pi0 = pi0$pi, T_fin = 149, begin_str = pi0$date[1], save_plot_data = TRUE, change_time = change_time)
model = tvt.eSIR(Y = pi0$Y, R = pi0$R, pi0 = pi, T_fin = 149, begin_str = pi0$date[1], save_plot_data = TRUE, change_time = change_time)
I_proj = colMeans(model$theta_pp[, , 2], na.rm = TRUE)
R_proj = colMeans(model$theta_pp[, , 3], na.rm = T)
as.Date(seq(as.Date(“2021/04/10”), by = “day”, length.out = 5), format(“%B-%d-%y”))
a = seq(as.Date(“2021/04/10"), by = “day”, length.out = 45)
b = c()
for(i in 1:length(a)){
  b = c(b, format(a[i], format = “%m/%d/%y”))
}
b = as.Date(b, format = “%m/%d/%y”)
#eSIR_df = data.frame(date = b, preds = model$Y_mean, ymin = model$Y_band$lower, ymax = model$Y_band$upper, actual = test_df$I/N)
#eSIR_df = model$data_comp
#ggplot(data = eSIR_df, aes(x = date, y = preds)) + scale_x_date(date_breaks = “1 week”, labels = scales::date_format(“%e %b”))+ geom_line(colour = “black”) + geom_ribbon(mapping = aes(x = date, ymin = ymin, ymax = ymax), alpha = 0.6, fill = “grey70”) + geom_line(mapping = aes(x = date, y = actual), colour = “blue”)
eSIR_df_2 = data.frame(model$data_comp, date = israel$date, actual = israel$I/N)
# ggplot(data = eSIR_df_2) + scale_x_date(date_breaks = “2 week”, labels = scales::date_format(“%e %b”))  +
#  geom_line(mapping = aes(x = date, y = mean),color= “mean”) +
#  geom_line(mapping = aes(x = date, y = actual), color = “actual”) +
#  scale_colour_manual(“”,breaks = c(“mean”,“actual”),values = c(“mean”=“red”,“actual”=“black”)) +
#  geom_ribbon(mapping = aes(x = date, ymin = lower, ymax = upper), alpha = 0.6, fill = “grey70") +
          #  geom_vline(xintercept=as.numeric(israel$date[106]), linetype=“dotted”, color = “blue”, size=1.5)
          #
          dat = c(eSIR_df_2$mean, eSIR_df_2$actual)
        type = c(rep(c(“Prediction”),149), rep(c(“Actual”),149))
        rib_upper = c(eSIR_df_2$upper, eSIR_df_2$actual)
        rib_lower = c(eSIR_df_2$lower, eSIR_df_2$actual)
        date_long = rep(eSIR_df_2$date,2)
        eSIR_df_long = data.frame(dat,type,rib_upper,rib_lower,date_long)
        p_infected = ggplot(data = eSIR_df_long)+ scale_x_date(date_breaks = “1 month”,
                                                               labels = scales::date_format(“%e %b”)) +
          geom_line(mapping = aes(x = date_long, y = dat*N, color = type)) +
          geom_ribbon(mapping = aes(x = date_long, ymin = rib_lower*N, ymax = rib_upper*N,fill=type),alpha=0.4)+
          scale_color_manual(values = c(“red”,“black”)) + labs(x = “Date”, y = “Cases”) +
          geom_vline(xintercept=as.numeric(israel$date[106]), linetype=“dotted”, color = “blue”, size=1.5)
        eSIR_df_2 = data.frame(model$data_comp_R, date = israel$date, actual = israel$R/N)
        #ggplot(data = eSIR_df_2) + scale_x_date(date_breaks = “2 week”, labels = scales::date_format(“%e %b”))  + geom_ribbon(mapping = aes(x = date, ymin = lower, ymax = upper), alpha = 0.6, fill = “grey70”) + geom_line(mapping = aes(x = date, y = mean)) + geom_line(mapping = aes(x = date, y = actual), colour = “red”) + geom_vline(xintercept=as.numeric(israel$date[106]), linetype=“dotted”,
        #color = “blue”, size=1.5)
        dat = c(eSIR_df_2$mean, eSIR_df_2$actual)
        type = c(rep(c(“Prediction”),149), rep(c(“Actual”),149))
        rib_upper = c(eSIR_df_2$upper, eSIR_df_2$actual)
        rib_lower = c(eSIR_df_2$lower, eSIR_df_2$actual)
        date_long = rep(eSIR_df_2$date,2)
        eSIR_df_long = data.frame(dat,type,rib_upper,rib_lower,date_long)
        p_removed = ggplot(data = eSIR_df_long)+ scale_x_date(date_breaks = “1 month”,
                                                              labels = scales::date_format(“%e %b”)) +
          geom_line(mapping = aes(x = date_long, y = dat*N, color = type)) +
          geom_ribbon(mapping = aes(x = date_long, ymin = rib_lower*N, ymax = rib_upper*N,fill=type),alpha=0.4)+
          scale_color_manual(values = c(“red”,“black”)) + labs(x = “Date”, y = “Removed”) +
          geom_vline(xintercept=as.numeric(israel$date[106]), linetype=“dotted”, color = “blue”, size=1.5)
        return(list(p_infected = p_infected, p_removed = p_removed))
}
