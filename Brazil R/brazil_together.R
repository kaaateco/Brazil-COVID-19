
# New Cases, Recoveries, Fatalities plot
source("brazil_data_read.R")

plot_ly(brazil, x = ~date, y =~I, type = "bar", hoverinfo = "text") %>%
  layout(barmode = "stack", title = list(xanchor = "left", x = 0), legend =
           list(orientation = "h", font = list(size = 16))) %>%
  plotly::config(toImageButtonOptions = list(width = NULL, height = NULL))


plot_ly(brazil, x = ~date, y =~R, type = "bar", hoverinfo = "text") %>%
  layout(barmode = "stack", title = list(xanchor = "left", x = 0), legend =
           list(orientation = "h", font = list(size = 16))) %>%
  plotly::config(toImageButtonOptions = list(width = NULL, height = NULL))

plot_ly(stack, x = ~date, y =~val, color = ~name,
        type = "bar", hoverinfo = "text") %>%
  layout(barmode = "stack", title = list(xanchor = "left", x = 0), legend =
           list(orientation = "h", font = list(size = 16))) %>%
  plotly::config(toImageButtonOptions = list(width = NULL, height = NULL))



# R estimation
source("estimate_tvr_brazil.R")




# LS SIR model fitting
source("SIR_model_fitting.R")




# MLE SIR model fitting
source("SIRmle.R")
