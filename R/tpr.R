# library(data.table)
# library(plotly)
# library(deSolve)
# library(outbreaks)
# library(gridExtra)
# library(arm)
# library(tidyverse)
# library(ggplot2)
#
# world_covid_data <- fread("C:\\Users\\katec\\Documents\\UMich BDSI-SIBS\\Infectious Disease\\owid-covid-data.csv")
# world_covid_data <-  world_covid_data %>%
#   filter(location == "Brazil") %>%
#   mutate(tpr = total_cases/total_tests) %>%
#   filter(date >= as.Date("2020-11-01"), date <= as.Date("2021-04-01"))
#
# ci = c("#C79999")
# mn = c("#7C0000")
# date_breaks = "1 month"
# base = ggplot() +
#   xlab("") +
#   scale_x_date(
#     date_breaks = date_breaks,
#     labels = scales::date_format("%b %Y")
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
#   geom_smooth(mapping = aes(x = date, y = tpr, color = colour),
#               data = world_covid_data, size = 1,color="green4", span = 0.2,
#               fill = "gold") +
#   labs(y="Test Positivity Rate")
#
# ggplotly(p1)
