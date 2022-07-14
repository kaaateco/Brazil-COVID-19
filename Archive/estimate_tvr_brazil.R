# Read raw COVID-19 data from JHU
# Gist: https://gist.github.com/mkleinsa/9b4fb698d348d5d4d2e475253334bb2b

# Github folder:
# https://github.com/CSSEGISandData/COVID-19/tree/master/csse_covid_19_data/csse_covid_19_time_series

# load packages
# Install if needed: install.packages("tidyverse")
library(tidyverse)
library(ggplot2)
library(plotly)

# confirmed
case_url = "https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_confirmed_global.csv"
confirmed = read_csv(case_url)

# deaths
death_url = "https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_deaths_global.csv"
deaths = read_csv(death_url)

# recoveries
recovery_url = "https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_recovered_global.csv"
recoveries = read_csv(recovery_url)

# Next steps: Inspect the data with
View(confirmed)

# Process the data with dplyr (or software of your choice)
# Data format should be
#      - each row is a day (with date)
#      - the columns are daily cases, deaths, recoveries.

# -----

# Which countries are present in the data?
confirmed %>% pull(`Country/Region`) %>% unique()

# Multiple regions under this country?
confirmed %>% filter(`Country/Region` == "Brazil") %>%
  pull(`Province/State`) %>% unique()

brazil_confirmed = confirmed %>%
  filter(`Country/Region` == "Brazil") %>% #, is.na(`Province/State`)
  dplyr::select(-`Province/State`, -Lat, -Long, -`Country/Region`)
brazil_confirmed = tibble(date = colnames(brazil_confirmed),
                           cases_total = as.numeric(brazil_confirmed[1, ]))

# process/merge deaths and recoveries the in the same fasion

# -----

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

# Transform data for plotting
as.Date("1/01/20", format = "%m/%d/%y")
brazil =
  brazil %>%
  mutate(date = as.Date(date, format = "%m/%d/%y"))

# daily case, death, recoveries
brazil =
  brazil %>%
  mutate(total_removed = deaths_total + recoveries_total,
         total_infected = cases_total - total_removed,
         day_index = 1:n(),
         daily_removed = total_removed - lag(total_removed),
         daily_infected = total_infected - lag(total_infected),
         I = daily_infected,
         R = daily_removed)

brazil =
  brazil %>%
  filter(I < 3e5)

#brazil =
 # brazil %>%
  #filter(date > as.Date("2021-08-05"))

# Static plots
ggplot(brazil) +
  aes(x = date, y = I) +
  geom_col()

brazil =
  brazil %>%
  filter(date > as.Date("2020-04-01"),
         date < as.Date("2020-06-01"))


# fixing negative values

impu <- function(index){
  braz_pos <- brazil %>% filter(daily_infected > 0)
  brazil$daily_infected(index) <- mean(braz$daily_infected(index-1), brazil$daily_infected(index+1))
}

a = brazil %>%
  mutate(daily_infected = ifelse(daily_infected < 0, (lag(daily_infected) + lead(daily_infected))/2, daily_infected))


brazil$daily_infected[13] <- mean(brazil$daily_infected[12], brazil$daily_infected[15])

# set the mean serial interval
serial_interval = 4

# for a five day window
t_start = seq(2, nrow(brazil) - 4)
t_end = t_start + 4

res <- EpiEstim::estimate_R(
  incid = brazil$daily_infected,
  method = "parametric_si",
  config = EpiEstim::make_config(list(
    mean_si             = serial_interval,
    std_si              = 2,
    si_parametric_distr = "G",
    t_start             = t_start,
    t_end               = t_end,
    seed                = 46342))
)

plot(res)

start_date = brazil$date[1]

# fancy plot
plt_data <- tibble(
  date_num = res$dates
) %>% left_join(
  res$R, by = c("date_num" = "t_end")
) %>%
  dplyr::select(
    date_num, t_start, r = `Mean(R)`, lower = `Quantile.0.025(R)`, upper = `Quantile.0.975(R)`
  ) %>%
  add_column(date = brazil$date) %>%
  mutate(
    text = paste0("Date: ", format(date, format = '%b %d'), "<br>R: ",
                  format(round(r, 2), nsmall = 2), "<br>CI: ",
                  paste0("[", format(round(lower, 2), nsmall = 2), ", ",
                         format(round(upper, 2), nsmall = 2), "]"))
  ) %>%
  filter(!is.na(r))

cap <- paste0("\uA9 COV-IND-19 Study Group. Last updated: ",
              format(Sys.Date(), format = "%b %e"), sep = ' ')
axis_title_font <- list(size = 16)
tickfont        <- list(size = 16)

p <- plot_ly(plt_data, x = ~date, y = ~r, type = "scatter", mode = "lines",
             line = list(color = "rgb(54, 163, 11)", width = 5),
             hoverinfo = "text",
             text   = ~text) %>%
  add_markers(data = plt_data, x = ~date, y = ~r, mode = "marker",
              marker = list(color = "rgb(38, 38, 38)", symbol = 3)) %>%
  add_ribbons(ymin = ~lower,
              ymax = ~upper,
              line = list(color = 'rgba(54, 163, 11, 0.05)'),
              fillcolor = 'rgba(54, 163, 11, 0.2)',
              hoverinfo = "none") %>%
  layout(
    title = list(text = cap, xanchor = "left", x = 0),
    xaxis = list(title = "Date", titlefont = axis_title_font,
                 tickfont = tickfont, zeroline = T),
    yaxis = list(title = "R(t)", titlefont = axis_title_font,
                 tickfont = tickfont, zeroline = T),
    shapes = list(
      type = "line", xref = "paper", yref = "data",
      x0 = 0, x1 = 1, y0 = 1, y1 = 1,
      line = list(color = "rgba(255, 153, 51, 0.5)")
    ),
    showlegend = FALSE
  ) %>%
  plotly::config(toImageButtonOptions = list(width = NULL, height = NULL))

p

