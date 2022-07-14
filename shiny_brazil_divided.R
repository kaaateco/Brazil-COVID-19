library(shiny)
library(dplyr)
library(plotly)

source("Shiny/brazil_data_read.R")
source("estimate_tvr_brazil_figure2.R")
source("trying_divide_times.R")
source("tpr.R")
# source("SIRmle.R")

ui <- navbarPage("COVID-19 in Brazil",
        tabPanel("Plots",
          fluidPage(theme = bslib::bs_theme(bootswatch = "united"),

            titlePanel("COVID-19 in Brazil 11/1/20-4/1/21"),
              mainPanel(
                plotlyOutput("casePlot"),
                plotOutput("rtPlot"),
                plotlyOutput("lssirPlot"),
                plotlyOutput("tprPlot")
                            ))),
        tabPanel("Metrics"),
        tabPanel("References")
)

server <- function(input, output) {

    brazil <- data_read("2020-11-01", "2021-04-01")

    stack <- bind_rows(
      brazil %>% dplyr::select(date, val = daily_new_cases) %>%
        mutate(name = "New Cases"),
      brazil %>% dplyr::select(date, val = daily_new_deaths) %>%
        mutate(name = "Fatalities"),
      brazil %>% dplyr::select(date, val = daily_new_recoveries) %>%
        mutate(name = "Recoveries")
    )

  # Full Case plotting

    # output$casePlot <- renderPlotly({
    #   plot_ly(dataplot(), x = ~date, y =~val, color = ~name,
    #           type = "bar") %>%
    #     layout(barmode = "stack", title = list(xanchor = "left", x = 0),
    #            xaxis = list(title = "Date", titlefont = axis_title_font),
    #            yaxis = list(title = "Daily Counts", titlefont = axis_title_font),
    #            legend = list(orientation = "v", font = list(size = 16)), hovermode = "x unified") %>%
    #     plotly::config(toImageButtonOptions = list(width = NULL, height = NULL))

  output$casePlot <- renderPlotly({
    p1 <- plot_ly(brazil, x = ~date, y =~I, text = ~text_cases, type = "bar", hoverinfo = "text") %>%
      layout(barmode = "stack", title = list(xanchor = "left", x = 0), legend =
               list(orientation = "h", font = list(size = 16))) %>%
      plotly::config(toImageButtonOptions = list(width = NULL, height = NULL))

    p2 <- plot_ly(brazil, x = ~date, y =~daily_new_recoveries, text = ~text_recoveries, type = "bar", hoverinfo = "text") %>%
      layout(barmode = "stack", title = list(xanchor = "left", x = 0), legend =
               list(orientation = "h", font = list(size = 16))) %>%
      plotly::config(toImageButtonOptions = list(width = NULL, height = NULL))

    p3 <- plot_ly(stack, x = ~date, y =~val, color = ~name, colors = c("mediumblue", "gold", "green4"),
            type = "bar", hoverinfo = "text") %>%
      layout(barmode = "stack", title = list(xanchor = "left", x = 0), legend =
               list(orientation = "h", font = list(size = 16)), hovermode = "x unified") %>%
      plotly::config(toImageButtonOptions = list(width = NULL, height = NULL))

    p3
  })


  # Rt Estimation
  output$rtPlot <- renderPlotly({p})


  # LS SIR Model
  brazil$R <- brazil$total_removed
  output$lssirPlot <- renderPlotly({plot_periods(bind, brazil)

  })

  # TPR Plot
  output$tprPlot <- renderPlotly({p1})

  #
  #
  # # MLE SIR Model
  # output$mlesirPlot <- renderPlotly({p1})
}

# Run the application
shinyApp(ui = ui, server = server)
