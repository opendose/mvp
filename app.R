library(shiny)
library(stringr)
library(lattice)
library(ggplot2)
library(dplyr)
library(data.table)
library(mrgsolve)
library(magrittr)
source("opendose.R")

mod <- mread("bpg_perscomm_hand2018")
mod <- opendose_modrewriteETA(mod) # recompiles the model with ETAn params
mod <- mod %>% update(delta=1/24/6)

# Define UI for app that draws a histogram ----
ui <- fluidPage(
  titlePanel("Bayesian BPG"),
  sidebarLayout(
    sidebarPanel(
      numericInput("bmi", "Body mass index (BMI; threshold=25)", 20, 0, 60),
      numericInput("ffm", "Fat-free mass (kg)", 40, 0, 200),
      hr(),
      numericInput("dose", "Dose (mg)", 900, 0),
      hr(),
      textAreaInput("tdm", "Drug levels (format each line as 'days ug/L')", value="", rows=3),
      hr(),
      checkboxInput("showensemble", "Ensemble of estimates from TDM")
    ),
    mainPanel(
      h2("output"),
      plotOutput("plot"),
      plotOutput("popPlot")
    )
  )
)



server <- function(input, output) {
  Rbasemod <- reactive(mod)
  
  Rtdm <- reactive({read.delim(text=input$tdm, sep=" ", header=FALSE, col.names=c("time", "Y"))})
  Rtdmactive <- reactive({nrow(Rtdm()) > 0})
  Rpriormod <- reactive({Rbasemod() %>% param(FFM=input$ffm, BMI=input$bmi) %>% ev(amt=input$dose * 1000) %>% update(end=30)})
  Rpriorsim <- reactive({Rpriormod() %>% mrgsim %>% as.data.frame})
  Rebe <- reactive({opendose_gamble(Rpriormod(), Rtdm()$time, Rtdm()$Y, methods="ebe")$ebe})
  Rensemble <- reactive({opendose_gamble(Rpriormod(), Rtdm()$time, Rtdm()$Y, methods=c("montecarlo_prior", "montecarlo_posterior"))})

  Rpharmacodynamic <- reactive({
    df <- Rpriorsim() %>% subset(F >= 20)

    if(nrow(df) >= 2) {
      alltime <- df$time
      l <- head(alltime, 1)
      r <- tail(alltime, 1)
      list(
        geom_rect(aes(xmin=l, xmax=r, ymin=0, ymax=20), fill='darkolivegreen2'),
        annotate(geom="text", x=(l+r)/2, y=10, label=sprintf('Days over 20ug/L = %.1f', r-l), color="black")
      )
    } else {
      list()
    }
  })
  
  Rplot <- reactive({
    theplot <- ggplot()
    
    for(el in Rpharmacodynamic()) {
      theplot <- theplot + el
    }
    
    if(Rtdmactive()) {
      if(input$showensemble) {
        ensemble <- Rensemble()$montecarlo_posterior %>% head(50)
        ensemble$ID <- seq(nrow(ensemble))
        modresults <- apply(ensemble, 1, function(etalist) {
          etalist %<>% as.list
          names(etalist) <- colnames(ensemble)
          theid <- etalist$ID
          etalist$ID <- NULL
          df <- Rpriormod() %>% param(etalist) %>% mrgsim %>% as.data.frame
          df$ID <- theid
          df
        })
        modresults <- rbindlist(modresults)
        modresults$ID <- as.factor(modresults$ID)
        theplot <- theplot + geom_line(data=modresults, aes(x=time, y=F, group=ID), color="gray", size=0.5)
      }
      theplot <- theplot + geom_line(data=Rpriormod() %>% param(Rebe()) %>% mrgsim %>% as.data.frame, aes(x=time, y=F), color='red')
      theplot <- theplot + geom_point(data=Rtdm(), aes(x=time, y=Y), color='green', size=3)
    }
    
    theplot <- theplot + geom_line(data=Rpriorsim(), aes(x=time, y=F), color='blue')
    
    theplot <- theplot +
      scale_x_continuous(name="Time (days)", expand=c(0, 0.1)) +
      scale_y_continuous(name="Level (ug/L)", limits=c(0,50), expand=c(0, 0.1))

    theplot
  })
  
  Rpopplot <- reactive({
    omega <- Rpriormod() %>% omat %>% as.matrix
    stdevs <- omega %>% diag %>% sqrt
    scal <- 3
    gd <- expand.grid(ETA1=seq(-scal*stdevs[1], scal*stdevs[1], scal*stdevs[1]/100),
                      ETA2=seq(-scal*stdevs[2], scal*stdevs[2], scal*stdevs[2]/100))
    gd$PRIORPDF <- dmvnorm(as.matrix(gd), sigma=omega)/dmvnorm(c(0,0), sigma=omega)
    theplot <- ggplot() + geom_raster(data=gd, aes(x=ETA1, y=ETA2, fill=PRIORPDF)) +
      scale_fill_gradientn(colours = topo.colors(100))
    
    if(Rtdmactive()) {
      if(input$showensemble) {
        theplot <- theplot + geom_point(data=Rensemble()$montecarlo_posterior, aes(x=ETA1, y=ETA2), color='gray', size=0.2)
      }
      theplot <- theplot + geom_point(data=as.data.frame(Rebe()), aes(x=ETA1, y=ETA2), color='red', size=5)
    }
    
    theplot <- theplot +
      scale_x_continuous(name="ETA1 (half-lives)", limits=c(-scal*stdevs[1], scal*stdevs[1]), expand=c(0, 0)) +
      scale_y_continuous(name="ETA2 (Vdist)", limits=c(-scal*stdevs[2], scal*stdevs[2]), expand=c(0, 0))

    theplot
  })
  
  output$plot <- renderPlot(Rplot())
  output$popPlot <- renderPlot(Rpopplot())
}


shinyApp(ui = ui, server = server)
