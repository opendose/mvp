library(shiny)
library(stringr)
library(lattice)
library(ggplot2)
library(dplyr)
library(data.table)
library(mrgsolve)
library(magrittr)
source("opendose.R")

ffm.model <- function(age, weight, height, male) {
  bmi <- weight / height^2
  if (male) {
    (0.88+(1-0.88)/(1+(age/13.4)^-12.7))*(9270*weight)/(6680+(216*bmi))
  } else {
    (1.11+(1-1.11)/(1+(age/7.1)^-1.1))*(9270*weight)/(8780+(244*bmi))
  }
}

mod <- mread("bpg_perscomm_hand2018")
mod <- opendose_modrewriteETA(mod) # recompiles the model with ETAn params
mod <- mod %>% update(delta=1/24/6)

ui <- fluidPage(
  title = "BPG",
  h1("Calculator: benzathine benzylpenicillin (BPG) in children"),
  hr(),
  fluidRow(
    column(4,
           wellPanel(
             h4("Patient"),
             conditionalPanel("!input.manualinfo",
                              numericInput("optionalage", "Age (years)", 14, 0, 20),
                              numericInput("optionalheight", "Height (cm)", 161, 0, 200),
                              numericInput("optionalweight", "Weight (kg)", 63, 0, 200),
                              radioButtons("optionalsex", "Sex", c("Female"="f", "Male"="m"), inline=T)
             ),
             conditionalPanel("input.manualinfo",
                              numericInput("optionalffm", "Fat-free mass (kg)", 40, 0, 200),
                              checkboxInput("optionaloverweight", "Overweight (BMI > 25 ∴ CL ↑86%)")
             )
           )
    ),
    column(4,
           wellPanel(
             h4("Dose"),
             numericInput("dose", "Intramuscular BPG (mg)", 450, 0, 3600, 450)
             # radioButtons("mxhistory", NULL, c("First dose"="first", "Steady state"="steady")),
             # conditionalPanel("input.mxhistory == 'steady'",
             #                  numericInput("interval", "Time between last two doses (days)", 21, 7, 49, 7)
             # )
           )
    ),
    column(4,
           wellPanel(
             h4(checkboxInput("tdmactive", "Level")),
             conditionalPanel("input.tdmactive",
               numericInput("tdmtime", "Time since last dose (days)", 0, 0, 28, 0.25),
               numericInput("tdmlevel", "Measured plasma BPG (µg/L)", 0, 0, 200, 0.1)
             )
           )
    )
  ),
  hr(),
  tabsetPanel(
    tabPanel("Results", plotOutput("plot")),
    tabPanel("Technical settings",
             conditionalPanel("input.password != 'opendose'",
                              passwordInput("password", "Password ('opendose')")
             ),
             conditionalPanel("input.password == 'opendose'",
                              checkboxInput("manualinfo", "Enter BMI and FFM manually"),
                              checkboxInput("showensemble", "Show Monte Carlo posterior estimates (slow)"),
                              htmlOutput("report"),
                              plotOutput("popplot")
             )
    )
  )
)

server <- function(input, output) {
  # individual covariates
  Rffm <- reactive({
    if (input$manualinfo) {
      input$optionalffm
    } else {
      ffm.model(
        age=input$optionalage,
        weight=input$optionalweight,
        height=input$optionalheight/100,
        male=(input$optionalsex == "m")
      )
    }
  })
  
  Roverweight <- reactive({
    as.numeric(
      if (input$manualinfo) {
        input$optionaloverweight
      } else {
        input$optionalweight/(input$optionalheight/100)^2 > 25
      }
    )
  })
  
  Rbasemod <- reactive(mod)

  Rtdm <- reactive(data.frame(time=input$tdmtime, Y=input$tdmlevel))
  Rtdmactive <- reactive(input$tdmactive)
  Rpriormod <- reactive({Rbasemod() %>% param(FFM=Rffm(), OVWT=Roverweight()) %>% ev(amt=input$dose * 1000) %>% update(end=30)})
  Rebe <- reactive({opendose_gamble(Rpriormod(), Rtdm()$time, Rtdm()$Y, methods="ebe")$ebe})
  Rensemble <- reactive({opendose_gamble(Rpriormod(), Rtdm()$time, Rtdm()$Y, methods=c("montecarlo_prior", "montecarlo_posterior"))})

  Rplot <- reactive({
    theplot <- ggplot()

    if(Rtdmactive() && input$showensemble) {
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

    if(Rtdmactive()) {
      theplot <- theplot + geom_point(data=Rtdm(), aes(x=time, y=Y), size=3)
    }

    if (Rtdmactive()) {
      best <- Rpriormod() %>% param(Rebe()) %>% mrgsim %>% as.data.frame
      best$method <- "With drug level"
      secondbest <- Rpriormod() %>% mrgsim %>% as.data.frame
      secondbest$method <- "Without drug level"
      bigtable <- rbind(best, secondbest)
      theplot <- theplot + theme(legend.position = c(1, 1), legend.justification = c(1, 1), legend.title=element_blank())
    } else {
      best <- Rpriormod() %>% mrgsim %>% as.data.frame
      best$method <- "No drug level"
      bigtable <- best
      theplot <- theplot + theme(legend.position="none")
    }
    bigtable$method %<>% factor

    
    df <- best %>% subset(F >= 20)
    
    if(nrow(df) >= 2) {
      alltime <- df$time
      l <- head(alltime, 1)
      r <- tail(alltime, 1)
      theplot <- theplot +
        geom_rect(aes(xmin=l, xmax=r, ymin=0, ymax=20), fill='darkolivegreen2') +
        annotate(geom="text", x=(l+r)/2, y=10, label=sprintf('Time over 20 µg/L = %.1f days', r-l), color="black")
    }
    
    theplot <- theplot + geom_line(data=bigtable, aes(x=time, y=F, linetype=method), size=1)
    

    theplot <- theplot +
      scale_x_continuous(name="Time since last dose (days)", expand=c(0, 0.1)) +
      scale_y_continuous(name="Plasma BPG concentration (µg/L)", limits=c(0,50), expand=c(0, 0.1))

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

  output$popplot <- renderPlot(Rpopplot())
  output$report <- reactive({paste(c(
    sprintf("Calculated FFM = %.2f kg", Rffm()),
    sprintf("Calculated BMI = %.1f", if (input$manualinfo) -1 else input$optionalweight/(input$optionalheight/100)^2),
    sprintf("Categorical overweight (>25) = %s", if (Roverweight()) "yes" else "no"),
    sprintf("ETA: %s", if (Rtdmactive()) Rebe() else "tdm not active"),
    ""
  ), collapse="<br />")})
}

shinyApp(ui = ui, server = server)
