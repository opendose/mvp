library(shiny)
library(stringr)
library(lattice)
library(ggplot2)
library(dplyr)
source("model.R")

# Define UI for app that draws a histogram ----
ui <- fluidPage(
  titlePanel("Bayesian BPG"),
  sidebarLayout(
    sidebarPanel(
      numericInput("bmi", "Body mass index (BMI; threshold=25)", 20, 0, 60),
      numericInput("ffm", "Fat-free mass (kg)", 40, 0, 200),
      hr(),
      h2("Prior (blue)"),
      numericInput("dose", "Dose (mg)", 900, 0),
      h2("Posterior (red)"),
      textAreaInput("tdm", "Drug levels ('days ug/L', one per line)", value="30 20", rows=3)
    ),
    mainPanel(
      h2("output"),
      plotOutput("plot"),
      plotOutput("plot2"),
      textOutput("report")
    )
  )
)


# Constants
simduration <- 30


server <- function(input, output) {
  Rbasemod <- reactive(mod.bpg_perscomm_hand2018)
  
  Rtdm <- reactive(
    read.delim(text=input$tdm, sep=" ", header=FALSE, col.names=c("t", "y"))
  )
  
  Rprior <- reactive(
    Rbasemod() %>% param(FFM=input$ffm, BMI=input$bmi) %>% ev(t=0, amt=input$dose * 1000)
  )
  
  Rfit <- reactive(
    Rprior() %>% tdm(y=Rtdm()$y, t=Rtdm()$t)
  )
  
  Rposterior <- reactive(
    Rprior() %>% hack.mod.for.fit(Rfit())
  )
  
  output$plot <- renderPlot(
    ggplot() +
      geom_line(data=Rprior() %>% drop.re %>% update(end=simduration) %>% mrgsim %>% as.data.frame, aes(x=time, y=DV), color='blue', label='Prior') +
      geom_line(data=Rposterior() %>% drop.re %>% update(end=simduration) %>% mrgsim %>% as.data.frame, aes(x=time, y=DV), color='red', label='Posterior') +
      geom_point(data=Rtdm(), aes(x=t, y=y), color='red', label='Drug levels')
  )
  
  output$plot2 <- renderPlot(
    
    Rposterior() %>% drop.re %>% update(end=simduration) %>% mrgsim %>% plot
  )
  
  output$report <- renderPrint({
    mypost <- Rposterior() %>% drop.re %>% mrgsim %>% as.data.frame %>% head(1)
    myprior <- Rprior() %>% drop.re %>% mrgsim %>% as.data.frame %>% head(1)
    
    allcols <- colnames(mypost)
    goodcols <- grep("^(V\\d?|CL)$", allcols, value=TRUE, perl=TRUE)
    
    mypost <- mypost[goodcols]
    myprior <- myprior[goodcols]
    
    for(col in goodcols) {
      print(paste(col, "=", mypost[col], "population mean", myprior[col]))
    }
  })
}


shinyApp(ui = ui, server = server)
