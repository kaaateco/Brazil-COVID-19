library(shiny)
library(dplyr)

source("brazil_data_read.R")

server <- plot_ly(brazil, x = ~date, y =~I, type = "bar", hoverinfo = "text") %>%
  layout(barmode = "stack", title = list(xanchor = "left", x = 0), legend =
           list(orientation = "h", font = list(size = 16))) %>%
  plotly::config(toImageButtonOptions = list(width = NULL, height = NULL))

# Define UI for application that draws a histogram
ui <- fluidPage(

  # Application title
  titlePanel("COVID-19 in Brazil 4/1/20-6/1/20"),

  imageOutput(
    outputId,
    width = "100%")

#   # Sidebar with a slider input for number of bins
#   sidebarLayout(
#     sidebarPanel(
#       sliderInput("bins",
#                   "Number of bins:",
#                   min = 1,
#                   max = 50,
#                   value = 30)
#     ),
#
#     # Show a plot of the generated distribution
     mainPanel(
       plotOutput("distPlot")
     )
   )
# )
#
# # Define server logic required to draw a histogram
server <- function(input, output) {

  output$distPlot <- renderPlot({
    # generate bins based on input$bins from ui.R
    x    <- faithful[, 2]
    bins <- seq(min(x), max(x), length.out = input$bins + 1)

    # draw the histogram with the specified number of bins
    hist(x, breaks = bins, col = 'darkblue', border = 'white')
  })
 }

# Run the application
shinyApp(ui = ui, server = server)
