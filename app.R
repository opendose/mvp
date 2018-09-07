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
      textAreaInput("tdm", "Drug levels (format each line as 'days ug/L')", value="", rows=3)
    ),
    mainPanel(
      h2("output"),
      plotOutput("plot"),
      fluidRow(
        column(6, plotOutput("eta1")),
        column(6, plotOutput("eta2"))
      ),
      textOutput("report")
    )
  )
)


# Constants
simduration <- 30


extractp <- function(mod, prm) {
  (mod %>% drop.re %>% update(end=0) %>% mrgsim %>% as.data.frame)[[prm]][1]
}


mkbell <- function(val, covar.name, mean, logstdev, covar.name2=NULL, mean2=NULL, logstdev2=NULL) {
  bellfunc <- function(x){dnorm(x=log(x/mean)/logstdev)}
  bellcurve <- stat_function(fun=bellfunc, geom="ribbon", mapping=aes(ymin=0, ymax=..y..), fill="grey")
  theplot <- ggplot() + bellcurve

  lside <- mean/exp(3.5*logstdev)
  rside <- mean*exp(3.5*logstdev)
  
  val <- val %>% pmin(rside) %>% pmax(lside) # keep the eta within the visible region
  theplot <- theplot + geom_vline(xintercept=val, color="black")
  
  if(is.null(covar.name2)) {
    theplot <- theplot + scale_x_log10()
  } else {
    theplot <- theplot + scale_x_log10(sec.axis=sec_axis(~.*5, name=covar.name2))
  }
  
  theplot <- theplot + aes(c(lside, rside)) + xlab(covar.name) + ylab("")
  
  theplot
}


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
  
  
  
  Rmainplot <- reactive({
    theplot <- ggplot() +
      geom_line(data=Rprior() %>% drop.re %>% update(end=simduration) %>% mrgsim %>% as.data.frame, aes(x=time, y=DV), color='blue', label='Prior') +
      geom_point(data=Rtdm(), aes(x=t, y=y), color='red', label='Drug levels')
    
    if(nrow(Rtdm()) > 0) {
      theplot <- theplot +
        geom_line(data=Rposterior() %>% drop.re %>% update(end=simduration) %>% mrgsim %>% as.data.frame, aes(x=time, y=DV), color='red', label='Posterior')
    }
    
    theplot
  })
  
  output$plot <- renderPlot(Rmainplot())
  
  Retasource <- reactive({if(nrow(Rtdm()) > 0) {Rposterior()} else {Rprior()}})
  
  output$eta1 <- renderPlot({
    mean <- extractp(Rprior(), "INDV")
    logstdev <- log(extractp(Rprior() %>% param(FORCEETA2=1), "INDV") / mean)
    val <- extractp(Retasource(), "INDV")
    mkbell(val, "Volume of distribution (L)", mean, logstdev)
  })
  
  output$eta2 <- renderPlot({
    mean <- extractp(Rprior(), "INDT12DIS")
    logstdev <- log(extractp(Rprior() %>% param(FORCEETA1=1), "INDT12DIS") / mean)
    val <- extractp(Retasource(), "INDT12DIS")
    mkbell(val, "Distribution half-life (days)", mean, logstdev, "Other covariate")
  })
  
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
