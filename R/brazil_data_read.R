data_read <- function(start, end){

# Read raw COVID-19 data from JHU

# Github folder:
# https://github.com/CSSEGISandData/COVID-19/tree/master/csse_covid_19_data/csse_covid_19_time_series

# load packages
library(tidyverse)
library(ggplot2)
library(plotly)
library(htmlwidgets)

# confirmed
case_url = "https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_confirmed_global.csv"
confirmed = read_csv(case_url)

# recovered
# Exercise: repeat this code for recoveries
recovered_url = "https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_recovered_global.csv"
recovered = read_csv(recovered_url)

# died
# Exercise: repeat this code for deaths
death_url = "https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_deaths_global.csv"
deaths = read_csv(death_url)

# Next steps: Inspect the data with
## View(confirmed)

# Process the data with dplyr (or software of your choice)
# Data format should be
#      - each row is a day (with date)
#      - the columns are daily cases, deaths, recoveries.

confirmed %>% pull(`Country/Region`) %>% unique()

confirmed %>% pull(`Country/Region`) %>% unique()


confirmed %>% filter(`Country/Region` == "Brazil") %>%
  pull(`Province/State`) %>% unique()

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

brazil_recovered = recovered %>%
  filter(`Country/Region` == "Brazil") %>% #, is.na(`Province/State`)
  dplyr::select(-`Province/State`, -Lat, -Long, -`Country/Region`)
brazil_recovered = tibble(date = colnames(brazil_recovered),
                           recoveries_total = as.numeric(brazil_recovered[1, ]))

brazil = brazil_confirmed %>% left_join(brazil_recovered, by = "date")
brazil = brazil %>% left_join(brazil_deaths, by = "date")

# Transform data for plotting
#as.Date("1/24/20", format = "%m/%d/%y")

brazil =
  brazil %>%
  mutate(date = as.Date(date, format = "%m/%d/%y"))

brazil =
  brazil %>%
  filter(date > as.Date(start)-1,
         date < as.Date(end))


# daily case, death, recoveries
brazil =
  brazil %>%
  mutate(total_removed = deaths_total + recoveries_total,
         daily_removed = total_removed - lag(total_removed),
         R = total_removed - lag(total_removed),
         I = cases_total - total_removed,
         daily_new_cases = cases_total- lag(cases_total),
         daily_new_recoveries = recoveries_total - lag(recoveries_total),
         daily_new_deaths = deaths_total - lag(deaths_total),
         day_index = 1:n()
         )

#brazil =
 # brazil %>% mutate(text_cases = paste0("Cases: ", daily_new_cases),
  #                  text_recoveries = paste0("Recoveries: ", daily_new_recoveries),
   #                 text_deaths = paste0("Deaths: ", daily_new_deaths))

 brazil =
   brazil %>% filter(I > 0, R > 0)



# Static plots
ggplot(brazil) +
  aes(x = date, y = I) +
  geom_col()

ggplot(brazil) +
  aes(x = date, y = R) +
  geom_col()


# stacked geom col:
stack = bind_rows(
  brazil %>% dplyr::select(date, val = daily_new_cases, ) %>%
    mutate(name = "New Cases"),
  brazil %>% dplyr::select(date, val = daily_new_deaths) %>%
    mutate(name = "Fatalities"),
  brazil %>% dplyr::select(date, val = daily_new_recoveries) %>%
    mutate(name = "Recoveries")
  %>% filter("New Cases" > 0)
)

stack <- stack %>% mutate(text = paste0(name, ": ", val))


ggplot(stack) +
  aes(x = date, y = val, fill = name) +
  geom_col()

return(list(brazil, stack))

}

brazil <- data_read("2020-11-01", "2021-05-01")[[1]]
stack <- data_read("2020-11-01", "2021-05-01")[[2]]

p3 <- plot_ly(stack, x = ~date, y =~val, color = ~name, text = ~text, colors = c("gold", "mediumblue", "green4"),
              type = "bar", hoverinfo = "text") %>%
  layout(barmode = "stack", title = list(xanchor = "left", x = 0), legend =
           list(orientation = "v", font = list(size = 16)), hovermode = "text",
         autosize = FALSE, width = 1125, height = 500,
         yaxis = list(title = "Daily Counts"), xaxis = list(title = "Date")) %>%
  plotly::config(toImageButtonOptions = list(width = NULL, height = NULL))

p3
saveWidget(p3, "daily.html")

#
# plot_ly(brazil, x = ~date, y =~I, type = "bar", hoverinfo = "text") %>%
#   layout(barmode = "stack", title = list(xanchor = "left", x = 0), legend =
#            list(orientation = "h", font = list(size = 16))) %>%
#   plotly::config(toImageButtonOptions = list(width = NULL, height = NULL))
#
#
# plot_ly(brazil, x = ~date, y =~R, type = "bar", hoverinfo = "text") %>%
#   layout(barmode = "stack", title = list(xanchor = "left", x = 0), legend =
#            list(orientation = "h", font = list(size = 16))) %>%
#   plotly::config(toImageButtonOptions = list(width = NULL, height = NULL))
#
# plot_ly(stack, x = ~date, y =~val, color = ~name,
#         type = "bar", hoverinfo = "text") %>%
#   layout(barmode = "stack", title = list(xanchor = "left", x = 0), legend =
#            list(orientation = "h", font = list(size = 16))) %>%
#   plotly::config(toImageButtonOptions = list(width = NULL, height = NULL))
