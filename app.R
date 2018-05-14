library(shiny)
library(stringr)
source("model.R")

# Define UI for app that draws a histogram ----
ui <- fluidPage(
  titlePanel("Bayesian vancomycin"),
  sidebarLayout(
    sidebarPanel(
      sliderInput("crcl", "Cockcroft-Gault est CrCL (mL/min)", 0, 200, 90),
      sliderInput("tbw", "Total body weight (kg)", 0, 200, 70),
      hr(),
      h2("Prior"),
      textAreaInput("doses", "Doses ('hours milligrams', one per line)", value="0 1500", rows=3),
      h2("Posterior"),
      textAreaInput("tdm", "Drug levels ('hours ml/L', one per line)", value="12 15", rows=3)
    ),
    mainPanel(
      h2("output"),
      textOutput("report"),
      plotOutput("plot")
    )
  )
)


server <- function(input, output) {
  Rdoses <- reactive(
    read.delim(text=input$doses, sep=" ", header=FALSE, col.names=c("t", "amt"))
  )

  Rtdm <- reactive(
    read.delim(text=input$tdm, sep=" ", header=FALSE, col.names=c("t", "y"))
  )
  
  Rprior <- reactive(
    mod.vanc_roberts2011 %>% param(TBW=input$tbw, CRCL=input$crcl) %>% ev(do.call("ev", Rdoses()))
  )
  
  Rfit <- reactive(
    Rprior() %>% tdm(y=Rtdm()$y, t=Rtdm()$t)
  )
  
  Rposterior <- reactive(
    Rprior() %>% hack.mod.for.fit(Rfit())
  )
  
  output$report <- renderPrint({
    print(Rposterior())
  })
  output$plot <- renderPlot(plot(Rposterior() %>% mrgsim(nid=20), DV~.))
}


shinyApp(ui = ui, server = server)
