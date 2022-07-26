library(shiny)
library(dplyr)
library(plotly)
library(bslib)
library(tidyverse)

ui <- navbarPage("COVID-19 in Brazil",
        tabPanel("About",
          fluidPage(theme = bslib::bs_theme(bootswatch = "flatly"),
            titlePanel("Exploring COVID-19 in Brazil"),
              mainPanel("Since the beginning of 2020, there have been 570 million cases of COVID-19 worldwide and over 6.38 million deaths. For countries that mishandled the pandemic early on, the impacts have been devastating. The leadership in Brazil had arguably one of the worst responses to COVID-19 in the world, with the country having the third highest amount of cases and the second highest amount of deaths. With public  data from the Johns Hopkins Coronavirus Resource Center, we are able to dive deeper into the case dynamics in Brazil and better understand how policy neglect and a failure of public leadership enabled one of the worst outbreaks of COVID-19 in the world."
              ))),
        tabPanel("Daily Cases",
            titlePanel("Daily COVID-19 Cases, Fatalities, and Recoveries"),
              mainPanel(
                includeHTML("daily.html")
               # plotlyOutput("casePlot")
              )),
        tabPanel("R(t)",
            titlePanel("R(t) Over Time"),
              mainPanel(
                includeHTML("tvr.html")
                #plotlyOutput("rtPlot")
              )),
        tabPanel("SIR",
              titlePanel("Observed vs. SIR Estimated Cases"),
              mainPanel(
                includeHTML("sir.html")
                #plotlyOutput("lssirPlot")
              )),
        tabPanel("SEIR",
                 titlePanel("Observed vs. SEIR Estimated Cases"),
                 mainPanel(
                   includeHTML("seir.html")
                  # plotlyOutput("lsseirPlot")
                 )),
        tabPanel("eSIR",
                 titlePanel("Observed vs. eSIR Estimated Cases"),
                # mainPanel(
                 #  plotlyOutput("lssirPlot")
                 ),
        tabPanel("Vaccination",
                 titlePanel("Vaccination Rate Over Time"),
                 mainPanel(
                   includeHTML("vax.html")
                   #plotlyOutput("vaxPlot")
              ))
)

server <- function(input, output) {

    # brazil <- data_read("2020-11-01", "2021-05-01")[[1]]
    # stack <- data_read("2020-11-01", "2021-05-01")[[2]]

  # Full Case plotting
  # output$casePlot <- renderPlotly({
  #   p3 <- plot_ly(stack, x = ~date, y =~val, color = ~name, text = ~text, colors = c("gold", "mediumblue", "green4"),
  #           type = "bar", hoverinfo = "text") %>%
  #     layout(barmode = "stack", title = list(xanchor = "left", x = 0), legend =
  #              list(orientation = "v", font = list(size = 16)), hovermode = "text",
  #               autosize = FALSE, width = 1125, height = 500,
  #            yaxis = list(title = "Daily Counts"), xaxis = list(title = "Date")) %>%
  #     plotly::config(toImageButtonOptions = list(width = NULL, height = NULL))
  #   p3
  #})


  # Rt Estimation
  #output$rtPlot <- renderPlotly({p})


  # Vaccination Plot
  #output$vaxPlot <- renderPlotly({vax_p1})


  # LS SIR Model
 # brazil$R <- brazil$total_removed
  #output$lssirPlot <- renderPlotly({plot_periods(bind, brazil)
  #})


  # LS SEIR Model
  #output$lsseirPlot <- renderUI({HTML("seir.html")})

    #renderPlotly({eplot_periods(pred_bind_I)


  # LS eSIR Model
  #output$lsesirPlot <- renderPlotly({eplot_periods(pred_bind_I)
 # })
}

# Run the application
shinyApp(ui = ui, server = server)
