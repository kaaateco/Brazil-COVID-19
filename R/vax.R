library(data.table)
library(plotly)
library(deSolve)
library(outbreaks)
library(gridExtra)
library(arm)
library(tidyverse)
library(ggplot2)
library(htmlwidgets)

world_covid_data <- fread("Data/owid-covid-data.csv")

world_covid_data <-  world_covid_data %>%
  filter(location == "Brazil") %>%
  mutate(vac = people_vaccinated, Population = 215633189, full = people_fully_vaccinated,
         first_dose = vac/population, full_dose = full/population) %>%
  filter(date >= as.Date("2020-01-01"), date <= as.Date("2022-07-01"))

#world_covid_data <- world_covid_data %>% mutate(text = paste0(name, ": ", val))

# world_covid_data = bind_rows(
#   world_covid_data %>% dplyr::select(date, val = first_dose) %>%
#     mutate(name = "First Dose"),
#   world_covid_data %>% dplyr::select(date, val = full_dose) %>%
#     mutate(name = "Second Dose"),
# )

vax_p1 <- plot_ly(world_covid_data, x = ~date, y = ~first_dose, name = "First Dose", type = 'scatter', mode = 'lines',
                line = list(color = 'mediumblue', width = 4, shape = "spline")) %>%
  add_trace(y = ~full_dose, name = 'Second Dose', line = list(color = "green", width = 4)) %>%
  layout(xaxis = list(title = "Date"),
          yaxis = list(title = "Proportion of Population", range = c(0, 1)),
         width = 1150, height = 475)

vax_p1
saveWidget(vax_p1, "vax.html")


# ci = c("#C79999")
# mn = c("#7C0000")
# date_breaks = "1 month"
#
# base = ggplot() +
#   xlab("") +
#   scale_x_date(
#     date_breaks = date_breaks,
#     labels = scales::date_format("%b %Y")
#   ) +
#   theme_bw() +
#   theme(
#     axis.text.x = element_text(hjust = 1),
#     axis.text = element_text(size = 12),
#     axis.title = element_text(size = 12)
#   ) +
#   theme(legend.position = "right")
#
# vax_p1 = base +
#   # geom_smooth(se = FALSE, mapping = aes(x = date, y = vac, color = colour),
#   #             data = world_covid_data, size = 1,color=mn, span = 0.2) +
#   geom_smooth(se = FALSE, mapping = aes(x = date, y = (vac / population)*100),
#               data = world_covid_data, size = 2,color="green4", span = 0.2) +
#   geom_smooth(se = FALSE, mapping = aes(x = date, y = (full / population)*100),
#               data = world_covid_data, size = 2,color="mediumblue", span = 0.2) +
#   labs(y = "Percent of Population Vaccinated")
#
# vax_p1 <- ggplotly(vax_p1, width = 1125, height = 525) %>%
#   add_trace(y = ~(vac / population), name = '1st Dose')
#
# #vax_p1.update_yaxes(range=[0, 20])
#
# vax_p1
