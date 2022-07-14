library(shiny)
library(dplyr)
library(plotly)

 source("brazil_data_read.R")
 source("estimate_tvr_brazil_figure2.R")
 source("SIR_model_fitting_func.R")
# source("SIRmle.R")

# Define UI for application
ui <- fluidPage(

  # Application title
  titlePanel("COVID-19 in Brazil 11/1/20-4/1/21"),

     mainPanel(
       plotlyOutput("casePlot"),
       plotOutput("rtPlot"),
       # plotlyOutput("lssirPlot"),
       # plotlyOutput("mlesirPlot")
       )
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

    ls_sir_pred <- ls_sir("2020-11-01", "2021-04-01", brazil)


  # Full Case plotting
  output$casePlot <- renderPlotly({
    p1 <- plot_ly(brazil, x = ~date, y =~I, type = "bar", hoverinfo = "text") %>%
      layout(barmode = "stack", title = list(xanchor = "left", x = 0), legend =
               list(orientation = "h", font = list(size = 16))) %>%
      plotly::config(toImageButtonOptions = list(width = NULL, height = NULL))

    p2 <- plot_ly(brazil, x = ~date, y =~R, type = "bar", hoverinfo = "text") %>%
      layout(barmode = "stack", title = list(xanchor = "left", x = 0), legend =
               list(orientation = "h", font = list(size = 16))) %>%
      plotly::config(toImageButtonOptions = list(width = NULL, height = NULL))

    p3 <- plot_ly(stack, x = ~date, y =~val, color = ~name,
            type = "bar", hoverinfo = "text") %>%
      layout(barmode = "stack", title = list(xanchor = "left", x = 0), legend =
               list(orientation = "h", font = list(size = 16))) %>%
      plotly::config(toImageButtonOptions = list(width = NULL, height = NULL))
    p3
  })

  # LS SIR Model
   output$lssirPlot <- renderPlot({
     # Plot results ----
     ci = c("#C79999")
     mn = c("#7C0000")
     date_breaks = "1 month"

     base = ggplot() +
       xlab("") +
       scale_x_date(
         date_breaks = date_breaks,
         labels = scales::date_format("%e %b")
       ) +
       theme_bw() +
       theme(
         axis.text.x = element_text(angle = 45, hjust = 1),
         axis.text = element_text(size = 12),
         axis.title = element_text(size = 12)
       ) +
       theme(legend.position = "right")

     print("p1")
     p1 = base +
       geom_line(mapping = aes(x = date, y = pred_I, color = colour),
                 data = ls_sir_pred, size = 1,color=mn) +
       geom_bar(mapping = aes(x = date, y = I), stat = "identity",
                data = data, width = 0.5, fill = 'steelblue', alpha = 0.7,
       ) + xlim(date_initial, date_final)

     p1 = p1 + labs(y = "Active Cases")
     #ggsave("Cases_8months.pdf",p1,width=8, height=6)

     p2 = base +
       geom_line(mapping = aes(x = date, y = pred_R, color = colour),
                 data = pred, size = 1,color = mn) +
       geom_bar(mapping = aes(x = date, y=R), stat="identity",
                data = data, width = 0.5, fill = 'steelblue', alpha = 0.7) +
       xlim(date_initial, date_final)

     p2 = p2 + labs(y = "Removed")
     #ggsave("Removed_8months.pdf",p2,width=8, height=6)

     p = grid.arrange(p1, p2)
     #ggsave("Both_6.png", p, width = 8, height = 6)

     p
   })


  # # Rt Estimation
  # output$lssirPlot <- renderPlotly({})
  #
  #
  # # MLE SIR Model
  # output$mlesirPlot <- renderPlotly({p1})
}

# Run the application
shinyApp(ui = ui, server = server)
