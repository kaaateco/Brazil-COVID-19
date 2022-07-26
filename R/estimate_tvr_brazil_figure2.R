
library(tidyverse)
library(ggplot2)
library(plotly)
library(RcppRoll)
library(htmlwidgets)

#source("Shiny/brazil_data_read.R")

# Transform data for plotting

as.Date("1/24/20", format = "%m/%d/%y")

start = "2020-11-01"
end = "2021-05-01"
brazil <- data_read(start, end)

brazil <- brazil[[1]]

brazil =
  brazil %>% mutate(date = as.Date(date, format = "%m/%d/%y"))

brazil = brazil[-1,]

# Static plots
# ggplot(brazil) +
#   aes(x = date, y = I) +
#   geom_col()



# Removing negative values

# brazil <- brazil %>%
#   mutate(daily_new_cases = ifelse(daily_new_cases < 0, NA, daily_new_cases)) %>%
#   mutate(daily_new_cases2 = stats::filter(daily_new_cases, rep(1/7, 7))) %>%
#   mutate(I = ifelse(I < 0, NA, I)) %>%
#   mutate(I2 = stats::filter(I, rep(1/7, 7)))
#
# n1 <- which(is.na(brazil$daily_new_cases))
# n2 <- which(!is.na(brazil$daily_new_cases2))
# n3 <- length(n1)
#
# for(i in 1:n3) {
#   xx <- max(which(n2 < n1[i]))
#   brazil$daily_new_cases[n1[i]] <- brazil$daily_new_cases2[n2[xx]]
# }
#
# m1 <- which(is.na(brazil$I))
# m2 <- which(!is.na(brazil$I2))
# m3 <- length(m1)

# for(i in 1:m3) {
#   yy <- max(which(m2 < m1[i]))
#   brazil$I[m1[i]] <- brazil$I2[m2[yy]]
# }


# set the mean serial interval
serial_interval = 4

# for a five day window
t_start = seq(2, nrow(brazil) - 4)
t_end = t_start + 4

res <- EpiEstim::estimate_R(
  incid = brazil$daily_new_cases,
  method = "parametric_si",
  config = EpiEstim::make_config(list(
    mean_si             = serial_interval,
    std_si              = 3,
    si_parametric_distr = "G",
    t_start             = t_start,
    t_end               = t_end,
    seed                = 1234))
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
  ) %>% filter(!is.na(r))

# cap <- paste0("\uA9 COV-IND-19 Study Group. Last updated: ",
            #  format(Sys.Date(), format = "%b %e"), sep = ' ')
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
  layout(width = 1125, height = 500,
    title = list(xanchor = "left", x = 0),
    xaxis = list(title = "Date", titlefont = axis_title_font,
                 tickfont = tickfont, zeroline = T),
    yaxis = list(title = "R(t)", titlefont = axis_title_font,
                 tickfont = tickfont, zeroline = T),
    shapes = list(
      type = "line", xref = "paper", yref = "data",
      x0 = 0, x1 = 1, y0 = 1, y1 = 1, width = 1000, height = 550,
      line = list(color = "mediumblue")
    ),
    showlegend = FALSE
  ) %>%
  plotly::config(toImageButtonOptions = list(width = NULL, height = NULL))

tvr_plot <- p

tvr_plot
saveWidget(tvr_plot, "tvr.html")

