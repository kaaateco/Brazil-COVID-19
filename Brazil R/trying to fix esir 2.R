# Data of COVID-19 can be found in the following R packages:
# removed = 6,282,184
# library(devtools)
# install_github("GuangchuangYu/nCov2019")
# library(nCov2019)
# install_github(qingyuanzhao/2019-nCov-Data)
# library("2019-nCov-Data")
set.seed(20192020)
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
         I = active_cases,
         R = total_removed)
brazil =
  brazil %>%
  filter(I < 2e7, R > -3e5)



brazil = brazil %>%
  filter(date >= as.Date("2020-11-01"), date <= as.Date("2021-04-01"))
brazil$day=c(1:nrow(brazil))

# Brazil province data Nov 01 -> April 01
# Cumulative number of infected cases
NI_complete <- brazil$cases_total
# Cumulative number of removed cases
RI_complete <- brazil$total_removed
N <- 215633189
R <- RI_complete / N
Y <- NI_complete / N - R #Nov 01 -> April 01

### Step function of pi(t)
change_time <- c("11/01/2020", "01/14/2021", "01/17/2021", "04/01/2021")
pi0 <- c(1.0, 0.9, 0.5, 0.5,0.5)
res.step <- tvt.eSIR(Y,
                     R,
                     begin_str = "11/01/2020",
                     T_fin = 181,
                     pi0 = pi0,
                     change_time = change_time,
                     dic = FALSE,
                     casename = "Brazil_step",
                     save_files = FALSE,
                     save_mcmc = FALSE,
                     save_plot_data = TRUE, # can plot them
                     M = 5e3,
                     nburnin = 2e3 # additional to M
)
#> The follow-up is from 01/13/20 to 07/30/20 and the last observed date is 02/11/20.
#> Running for step-function pi(t)
#> Compiling model graph
#>    Resolving undeclared variables
#>    Allocating nodes
#> Graph information:
#>    Observed stochastic nodes: 60
#>    Unobserved stochastic nodes: 36
#>    Total graph size: 1875
#>
#> Initializing model

res.step$plot_infection

# If we change the time unit to 5 days
# Cumulative number of infected cases
NI_complete2 <- NI_complete[seq(1, length(NI_complete), 5)]
# Cumulative number of removed cases
RI_complete2 <- RI_complete[seq(1, length(RI_complete), 5)]
N2 <- 215633189
R2 <- RI_complete2 / N2
Y2 <- NI_complete2 / N2 - R2 #Nov 01 -> April 01
change_time2 <- c("11/01/2020", "01/14/2021", "01/17/2021", "04/01/2021")
pi02 <- c(1.0, 0.9, 0.5, 0.5,0.5)
res.step2 <- tvt.eSIR2(Y2,
                       R2,
                       begin_str = "11/01/2020",
                       T_fin = 181,
                       pi0 = pi02,
                       change_time = change_time,
                       dic = FALSE,
                       casename = "Brazil_step2",
                       save_files = FALSE,
                       save_mcmc = FALSE,
                       save_plot_data = TRUE,
                       M = 5e3,
                       nburnin = 2e3,
                       time_unit = 5 #new feature!
)
#> The follow-up is from 01/13/20 to 07/26/20 and the last observed date is 02/07/20.
#> Running for step-function pi(t)
#> Compiling model graph
#>    Resolving undeclared variables
#>    Allocating nodes
#> Graph information:
#>    Observed stochastic nodes: 12
#>    Unobserved stochastic nodes: 12
#>    Total graph size: 395
#>
#> Initializing model

res.step2$plot_infection

res.step2$plot_removed

res.step2$spaghetti_plot

res.step2$dic_val
#> NULL
#res.step$gelman_diag_list
### continuous exponential function of pi(t)
res.exp <- tvt.eSIR(Y,
                    R,
                    begin_str = "11/01/2020",
                    death_in_R = brazil$deaths_total/(brazil$recoveries_total + brazil$deaths_total), #0.4
                    T_fin = 181,
                    exponential = TRUE,
                    dic = FALSE,
                    lambda0 = 0.05,
                    casename = "Brazil_exp",
                    save_files = FALSE,
                    save_mcmc = FALSE,
                    save_plot_data = TRUE,
                    M = 5e3,
                    nburnin = 2e3)
#> The follow-up is from 01/13/20 to 07/30/20 and the last observed date is 02/11/20.
#> Running for exponential-function pi(t)
#> Compiling model graph
#>    Resolving undeclared variables
#>    Allocating nodes
#> Graph information:
#>    Observed stochastic nodes: 60
#>    Unobserved stochastic nodes: 36
#>    Total graph size: 1875
#>
#> Initializing model

res.exp$plot_infection

res.exp$spaghetti_plot
