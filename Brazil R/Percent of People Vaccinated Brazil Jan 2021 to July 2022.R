library(data.table)
library(plotly)
library(deSolve)
library(outbreaks)
library(gridExtra)
library(arm)
library(tidyverse)
library(ggplot2)
world_covid_data <- fread("C:\\Users\\shado\\Downloads\\owid-covid-data.csv")
world_covid_data <-  world_covid_data %>% 
  filter(location == "Brazil") %>% 
  mutate(vac = people_vaccinated, Population = 215633189, full = people_fully_vaccinated) %>% 
  filter(date >= as.Date("2021-01-17"), date <= as.Date("2022-07-01"))

ci = c("#C79999")
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
  # geom_smooth(se = FALSE, mapping = aes(x = date, y = vac, color = colour),
  #             data = world_covid_data, size = 1,color=mn, span = 0.2) +
  geom_smooth(se = FALSE, mapping = aes(x = date, y = (vac / population)*100, color = colour),
             data = world_covid_data, size = 2,color="Green", span = 0.2) +
  geom_smooth(se = FALSE, mapping = aes(x = date, y = (full / population)*100, color = colour),
            data = world_covid_data, size = 2,color="Blue", span = 0.2) +
  labs(y = "Percent of Population Vaccinated")

p1 
