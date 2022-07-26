library(shiny)
library(dplyr)
library(plotly)
library(bslib)
library(tidyverse)

#source("Shiny/brazil_data_read.R")
#source("Shiny/estimate_tvr_brazil_figure2.R")
#source("Shiny/trying_divide_times.R")
#source("Shiny/SEIR_divided_model.R")

ui <- navbarPage("COVID-19 in Brazil (11/1/20 - 4/1/21)",
        tabPanel("Daily Cases",
          fluidPage(theme = bslib::bs_theme(bootswatch = "simplex"),

            titlePanel("Daily COVID-19 Cases, Fatalities, and Recoveries"),
              mainPanel(
                plotlyOutput("casePlot")
            #    plotlyOutput("tprPlot")
              ))),
        tabPanel("Instantaneous R",
              mainPanel(
                plotlyOutput("rtPlot")
              )),
        tabPanel("SIR",
              mainPanel(
                plotlyOutput("lssirPlot")
              ))
)

server <- function(input, output) {

    brazil <- data_read("2020-11-01", "2021-04-01")[[1]]
    stack <- data_read("2020-11-01", "2021-04-01")[[2]]

  # Full Case plotting

  output$casePlot <- renderPlotly({

    p3 <- plot_ly(stack, x = ~date, y =~val, color = ~name, text = ~text, colors = c("mediumblue", "gold", "green4"),
            type = "bar", hoverinfo = "text") %>%
      layout(barmode = "stack", title = list(xanchor = "left", x = 0), legend =
               list(orientation = "v", font = list(size = 16)), hovermode = "text",
                autosize = FALSE, width = 1000, height = 550,
             xaxis = list(title = "Date")), yaxis = list(title = "Daily Counts")) %>%
      plotly::config(toImageButtonOptions = list(width = NULL, height = NULL))
    p3
  })


  # Rt Estimation
  output$rtPlot <- renderPlotly({p})


  # LS SIR Model
  brazil$R <- brazil$total_removed
  output$lssirPlot <- renderPlotly({plot_periods(bind, brazil)
  })


  # LS SIR Model
  output$lsseirPlot <- renderPlotly({plot_periods(bind, brazil)
  })
}

# Run the application
shinyApp(ui = ui, server = server)
