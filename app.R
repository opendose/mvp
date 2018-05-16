library(shiny)
library(stringr)
library(lattice)
library(ggplot2)
library(dplyr)
source("model.R")

# Define UI for app that draws a histogram ----
ui <- fluidPage(
  titlePanel("Bayesian vancomycin"),
  sidebarLayout(
    sidebarPanel(
      radioButtons("modsel", "Model source", choices=c("mod.vanc_roberts2011", "mod.vanc_thomson2009")),
      hr(),
      numericInput("crcl", "Creatinine clearance (mL/min; urinary for Roberts, C-G for Thomson)", 90, 0, 200),
      numericInput("tbw", "Total body weight (kg)", 70, 0, 200),
      hr(),
      h2("Prior (blue)"),
      textAreaInput("doses", "Doses ('hours mg', one per line)", value="0 1500", rows=3),
      h2("Posterior (red)"),
      textAreaInput("tdm", "Drug levels ('hours mg/L', one per line)", value="12 15", rows=3)
    ),
    mainPanel(
      h2("output"),
      textOutput("report"),
      plotOutput("plot")
    )
  )
)


server <- function(input, output) {
  Rbasemod <- reactive(get(input$modsel))
  
  Rdoses <- reactive(
    read.delim(text=input$doses, sep=" ", header=FALSE, col.names=c("t", "amt")) %>% mutate(rate=600)
  )

  Rtdm <- reactive(
    read.delim(text=input$tdm, sep=" ", header=FALSE, col.names=c("t", "y"))
  )
  
  Rprior <- reactive(
    Rbasemod() %>% param(TBW=input$tbw, CRCL=input$crcl) %>% ev(do.call("ev", Rdoses()))
  )
  
  Rfit <- reactive(
    Rprior() %>% tdm(y=Rtdm()$y, t=Rtdm()$t)
  )
  
  Rposterior <- reactive(
    Rprior() %>% hack.mod.for.fit(Rfit())
  )
  
  output$plot <- renderPlot(
    ggplot() +
      geom_line(data=Rprior() %>% remove.mod.uncertainty %>% mrgsim %>% as.data.frame, aes(x=time, y=DV), color='blue', label='Prior') +
      geom_line(data=Rposterior() %>% remove.mod.uncertainty %>% mrgsim %>% as.data.frame, aes(x=time, y=DV), color='red', label='Posterior') +
      geom_point(data=Rtdm(), aes(x=t, y=y), color='red', label='Drug levels')
  )
}


shinyApp(ui = ui, server = server)
