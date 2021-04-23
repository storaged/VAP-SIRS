library(shiny)
library(shinyWidgets)
library(EpiDynamics)
library(ggplot2)
library(stringr)
library(deSolve)
library(MASS)
library(rvest) 
library(gridExtra)
library(reshape2)
library(zoo)
library(grid)
library(RColorBrewer)
library(dplyr)
library(DT)
library(shinydashboard)
library(promises)
library(future)
library(Cairo)
library(ggpubr)
library(tidyverse)
library(scales)
library(MullerPlot)
library(ggmuller)
source("html_helpers.R")


plan(multisession)


CONSTANT <- list(beta=1/6,gamma=1/6)

GAMMA <- 1/6
R_MAX <- 4
F_MAX <- 1; F_MIN <- 0
BETA_MAX <- R_MAX * GAMMA
TRANSMISSION <- BETA_MAX * F_MAX
NU_2 <- 500
# Define UI for application that draws a histogram

ui <- dashboardBody(
  tabBox(width = 12,
         tabPanel("VIP-SIRS",
                  #fluidPage(
                    fluidPage(
                      get_html_head() , shinyjs::useShinyjs(),
                      titlePanel("SIRS model with vaccination"),
                      get_short_description(),
                      fluidRow(
                        column(4,
                               mainPanel(h4("Vaccination:")),
                               checkboxInput("fast_vacc", "fast (100% in 200 days)", value = F),
                               checkboxInput("variant", "Fast loss of immunity (~200 days)", value = F)),
                        column(4,
                               mainPanel(h4("Only two betas:")),
                               checkboxInput("two_betas", withMathJax("Use only $\\beta_{S_V,I_V}^{I_V}$ others are equal") , value = T),
                               checkboxInput("fix_gamma", withMathJax("Fix: $\\gamma = \\frac{1}{6}$") , value = T),
                               checkboxInput("fix_reduction", withMathJax("Fix: $r = 0.8$") , value = F),
                               checkboxInput("fix_nu1star", withMathJax("Fix: $\\nu_1^{*} = \\nu_1$") , value = F)),
                        column(4,
                               mainPanel(h4("Contacts:")),
                               checkboxInput("check1", "Use means", value = F),
                               checkboxInput("go_crazy", withMathJax("Go crazy: $\\beta_{S_V,I_V}^{I_V} = 4 \\cdot \\beta_{S,I}^{I}$"), value = F))
                        ),
                    
                      fluidRow(
                        mainPanel(h4("Selected parameters' values: "))
                        ),
                      uiOutput("current_parameters"),
                      fluidRow( 
                        column(4,
                               checkboxInput("log_scale", "Use log-scale on Y axis", value = F)),
                        column(4,
                               selectInput("plot_type", "Select plot type:",
                                           c("Daily incidents" = "daily",
                                             "Infections per 1M" = "per1M"))),
                        column(4,
                               checkboxInput("proportional_mixing", "Use proportional mixing", value = F))
                      ),
                      fluidRow(
                        mainPanel(h4("Representative signature settings: "))
                      ),
                      fluidRow( 
                        column(2,
                               checkboxInput("signatureI", "Signature I (- +)", value = F)),
                        column(2,
                               checkboxInput("signatureII", "Signature II (+ - +)", value = F)),
                        column(2,
                               checkboxInput("signatureIII", "Signature III (+)", value = F)),
                        column(2,
                               checkboxInput("signatureIV", "Signature IV (-)", value = F)),
                        column(2,
                               checkboxInput("signatureV", "Signature V (+ -)", value = F))
                      ),
                      fluidRow(
                        plotOutput("extDistPlotv2", height = "800px")
                        ),
                      fluidRow(
                        column(2, 
                               downloadButton("downloadData", "Download .csv")),
                        column(2, downloadButton("downloadPlot", "Download plot"))
                      ),
                      
                      
                      fluidRow(
                        column(2,
                               sliderInput("infectivity",withMathJax("$t$") , min = 0,max = 1,value = TRANSMISSION, step = 0.001),
                               
                               sliderInput("beta11",withMathJax("$\\beta_{S_N,I_N}^{I_N}$"), min = 0,max = 1,value = 1/7, step = 0.001),
                               sliderInput("beta12",withMathJax("$\\beta_{S_N,I_N}^{I_V}$") , min = 0,max = 1,value = 2.5/7, step = 0.001),
                               sliderInput("beta13",withMathJax("$\\beta_{S_N,I_N}^{I_D}$") , min = 0,max = 1,value = 2.5/7, step = 0.001)
                               
                               ),
                        column(2,
                               sliderInput("freedom",withMathJax("$f$") , min = 0,max = 1,value = 0.37, step = 0.01),
                               
                               sliderInput("beta21",withMathJax("$\\beta_{S_V,I_V}^{I_N}$"), min = 0,max = 1,value = 2.5/7, step = 0.001),
                               sliderInput("beta22",withMathJax("$\\beta_{S_V,I_V}^{I_V}$"), min = 0,max = 1,value = 4/7, step = 0.001),
                               sliderInput("beta23",withMathJax("$\\beta_{S_V,I_V}^{I_D}$"), min = 0,max = 1,value = 4/7, step = 0.001)
                               
                               ),
                        column(2,
                               sliderInput("freedom_V",withMathJax("$f_V$") , min = 0,max = 1,value = 0.95, step = 0.01),
                              
                               sliderInput("beta31",withMathJax("$\\beta_{S_D,I_D}^{I_N}$"), min = 0,max = 1,value = 2.5/7, step = 0.001),
                               sliderInput("beta32",withMathJax("$\\beta_{S_D,I_D}^{I_V}$"), min = 0,max = 1,value = 4/7, step = 0.001),
                               sliderInput("beta33",withMathJax("$\\beta_{S_D,I_D}^{I_D}$"), min = 0,max = 1,value = 4/7, step = 0.001)
                               
                               ),
                        column(2,
                               sliderInput("reduction",withMathJax("$r$") , min = 0,max = 1,value = 1, step = 0.001), 
                               sliderInput("nu1",withMathJax("$\\nu_1$"), min = 0,max = 0.05, value = 0.004, step = 0.0001),
                               sliderInput("nu1star",withMathJax("$\\nu_1^{*}$"), min = 0,max = 0.05,value = 0.004, step = 0.0001),
                               sliderInput("nu2",withMathJax("$\\nu_2$"), min = 0,max = 0.05,value = 1/NU_2, step = 0.0001),
                               sliderInput("kappa", withMathJax("$\\kappa$"),  min = 0,max = 0.1,value = 1/180, step = 0.0001)
                               
                               ),
                        column(2,
                               sliderInput("sceptics", withMathJax("$a\\, (\\xi)$"),  min = 0,max = 1,value = 0.1, step = 0.01),
                               sliderInput("gamma", withMathJax("$\\gamma$"),  min = 0,max = 1,value = GAMMA, step = 0.001),
                               sliderInput("alpha1", withMathJax("$\\alpha_1$"),  min = 0,max = 1,value = 0.9, step = 0.01),
                               sliderInput("alpha2", withMathJax("$\\alpha_2$"),  min = 0,max = 1,value = 0.9, step = 0.01)
                               ),
                        column(2,
                               sliderInput("S", withMathJax("$S = (1-\\xi ) \\cdot S_N + \\xi \\cdot S_D$"), min = 0,  max = 100000, value = 100000, step=1000),
                               sliderInput("Sv1", withMathJax("$S_V$"), min = 0,  max = 100000, value = 0, step=1000),
                               sliderInput("I", withMathJax("$I = (1-\\xi ) \\cdot I_N + \\xi \\cdot I_D$"),   min = 0, max = 20000,  value = 1, step=1),
                               sliderInput("R", withMathJax("$R = (1-\\xi ) \\cdot R_N + \\xi \\cdot R_D$"),  min = 0,  max = 100000, value = 0, step=100 ),
                               sliderInput("V", withMathJax("$V$"),  min = 0,  max = 100000, value = 0 , step=100 ),
                               
                               sliderInput("lookup", withMathJax("$t_i$"),  min = 0,  max = 2*365, value = 400, step=1 )
                               )
                        )
                      )
                    #)
                  ),
                  
         tabPanel("VIP-SIRS: Description",
                  fluidPage(
                    get_html_head()
                    , shinyjs::useShinyjs(),
                    titlePanel("VIP-SIRS: Formulation of the model"),
                    
                    fluidRow(
                      mainPanel(h4("System of ODEs: "))
                    ),
                    fluidRow(
                      column(6,
                             uiOutput("ext_current_ode_p1")),
                      column(6,
                             uiOutput("ext_current_ode_p2"))
                    ), 
                    fluidRow(
                      mainPanel(h4("Parameters description: "))
                    ),
                    uiOutput("description")
                    )
         )
    )
)

go_crazy <- function(input, output, session){
  updateSliderInput(session, "beta22", value = min(1, 4 * input$beta11) ) 
  average_betas(input, output, session)
}

signatures_df <- data.frame(f = 1 - c(0.8, 0.63, 0.3, 0.8, 0.5), 
                            fv = 1 - c(0.05, 0.05, 0.05, 0.4, 0.4))
signatures_f <- lapply(1:5, function(i) 
  function(input, output, session){ 
    updateSliderInput(session, "freedom", value = signatures_df[i,"f"]); 
    updateSliderInput(session, "freedom_V", value = signatures_df[i,"fv"])  
})
names(signatures_f) <- paste("signature", as.roman(1:5), sep = "")

average_betas <- function(input, output, session){
  updateSliderInput(session, "beta12", value = mean(c(input$beta11, input$beta22))) 
  updateSliderInput(session, "beta21", value = mean(c(input$beta11, input$beta22))) 
}

# Define server logic required to draw a histogram
server <- function(input, output, session) {
  
  observeEvent(input$signatureI,{ if (input$signatureI) signatures_f$signatureI(input, output, session) })
  observeEvent(input$signatureII,{ if (input$signatureII) signatures_f$signatureII(input, output, session) })
  observeEvent(input$signatureIII,{ if (input$signatureIII) signatures_f$signatureIII(input, output, session) })
  observeEvent(input$signatureIV,{ if (input$signatureIV) signatures_f$signatureIV(input, output, session) })
  observeEvent(input$signatureV,{ if (input$signatureV) signatures_f$signatureV(input, output, session) })
  

  observeEvent(input$check1,{
    if (input$check1 & !input$go_crazy) { 
      toggleState("beta12") 
      toggleState("beta21")  
      average_betas(input, output, session)
    } else { 
      toggleState("beta12") 
      toggleState("beta21") 
    }
  })
  
  observeEvent(input$SIR.stage,{
    if (input$SIR.stage) { 
      shinyjs::show(id="SIR.beta2") 
      shinyjs::show(id="SIR.gamma2") 
      shinyjs::show(id="SIR.kappa2") 
      shinyjs::show(id="SIR.restriction") 
      #average_betas(input, output, session)
    } else { 
      shinyjs::hide(id="SIR.beta2") 
      shinyjs::hide(id="SIR.gamma2") 
      shinyjs::hide(id="SIR.kappa2") 
      shinyjs::hide(id="SIR.restriction") 
    }
  })
  
  observeEvent(input$SR.solution,{
    if (input$SR.solution) { 
      shinyjs::disable(id="SR.kappa") 
      updateSliderInput(session, "SR.kappa", value =  0 ) 
      shinyjs::show(id="SR.S.equation") 
      output$SR.S.equation <- renderUI({ 
        withMathJax(paste0("Rozwiązanie równania: ", 
                          '$$\\frac{dS}{dt}   = - \\gamma \\cdot S:', 
                          "\\quad S = e^{-\\gamma \\cdot \\ t + C}$$")) 
        })
    } else { 
      shinyjs::enable(id="SR.kappa") 
      shinyjs::hide(id="SR.S.equation") 
    }
  })
  
  
  
  
  
  observeEvent(input$go_crazy,{
    if (input$go_crazy) { 
      toggleState(id="beta22") 
      toggleState(id="beta12") 
      toggleState(id="beta21") 
      go_crazy(input, output, session)
    } else {
      toggleState(id="beta22") 
      toggleState(id="beta12") 
      toggleState(id="beta21") 
    }
    
  })
  
  update_betas <- function(input, output, session){
    updateSliderInput(session, "beta12", value =  input$beta11 ) 
    updateSliderInput(session, "beta13", value =  input$beta11 ) 
    updateSliderInput(session, "beta21", value =  input$beta11 ) 
    updateSliderInput(session, "beta23", value =  input$beta11 ) 
    updateSliderInput(session, "beta31", value =  input$beta11 ) 
    updateSliderInput(session, "beta32", value =  input$beta11 ) 
    updateSliderInput(session, "beta33", value =  input$beta11 ) 
  }
  
  observeEvent(input$two_betas,{
    toggleState(id="beta12") 
    toggleState(id="beta13")
    toggleState(id="beta21")
    toggleState(id="beta23")
    toggleState(id="beta31")
    toggleState(id="beta32")
    toggleState(id="beta33")
  })
  
  update_on_freedom <- function(input, output, session){
    updateSliderInput(session, "beta11", value =  input$freedom * input$infectivity ) 
    updateSliderInput(session, "beta22", value =  input$reduction * input$freedom_V * input$infectivity ) 
  }
  
  observeEvent(input$fix_gamma,{
    if (input$two_betas) { 
      if(input$fix_gamma){
        shinyjs::disable(id="gamma") 
        updateSliderInput( session, "gamma", value = 1/6 ) 
      } else {
        shinyjs::enable(id="gamma")
      }
    }
  })
  
  observeEvent(input$fix_reduction,{
    if (input$two_betas) { 
      if(input$fix_reduction) {
        shinyjs::disable(id="reduction") 
        updateSliderInput( session, "reduction", value = 0.8 ) 
      } else {
        shinyjs::enable(id="reduction")
      }
    }
  })
  
  observeEvent(input$fix_nu1star,{
      if(input$fix_nu1star) {
        shinyjs::disable(id="nu1star") 
        updateSliderInput( session, "nu1star", value = input$nu1 ) 
      } else {
        shinyjs::enable(id="nu1star")
      }
  })
  
  observeEvent(input$beta11,{
    if (input$check1 & !input$go_crazy) { 
      average_betas(input, output, session)
    } 
    if (input$go_crazy) { 
      go_crazy(input, output, session)
    } 
    if (input$two_betas) { 
      update_betas(input, output, session)
    } 
  })
  
  observeEvent(input$nu1,{
    if (input$fix_nu1star) { 
      updateSliderInput(session, "nu1star", value = input$nu1 ) 
    } 
  })
  
  observeEvent(input$freedom | input$freedom_V | input$infectivity | input$reduction ,{
    if (input$two_betas) { 
      update_on_freedom(input, output, session)
    } 
  })
  
  observeEvent(input$beta22,{
    if (input$check1 & !input$go_crazy) { 
      average_betas(input, output, session)
    } 
  })
  
  observeEvent(input$fast_vacc,{
    if (input$fast_vacc) { 
      updateSliderInput( session, "nu1", value = 1/200 ) 
    } #else {
      #updateSliderInput( session, "nu1", value = 1/800 ) 
    #}
  })
  
  observeEvent(input$variant,{
    if (input$variant) { 
      updateSliderInput( session, "nu2", value = 1/200 ) 
    } #else {
      #updateSliderInput( session, "nu2", value = 1/400 ) 
    #}
  })
  
  
  
  output$description <- renderUI({
    fluidRow(
      column(6, 
             withMathJax(
               sprintf('$$\\beta_{X,Y}^{Z} \\hbox{ - contact rate, influence of } Z \\hbox{ onto the flow from } X \\hbox{ to } Y $$\n
                       $$t \\hbox{ - the infectivity of the dominant virus variant} $$ \n
                       $$f_V, f \\hbox{ - restrictions-related freedom of groups: vaccinated, others} $$ \n
                       $$\\gamma \\hbox{ - recovery rate } $$ \n
                       $$\\kappa \\hbox{ - becoming susceptible after recovery rate } $$ ')
               )
               ),
      column(6, 
             withMathJax(
               sprintf('
                       $$\\nu_1 \\hbox{ - vaccination rate }$$ \n
                       $$\\nu_2 \\hbox{ - loss of immunity rate }$$ \n
                       $$\\xi \\hbox{ - fraction of population that deny vaccination } $$ \n
                       $$\\alpha_i \\hbox{ - vaccination efficacy}$$' )
               )
               )
               )
  })
  
  output$current_parameters <- renderUI({
    fluidRow(
      column(6, 
             withMathJax(
               sprintf('$$\\beta_{S_N,I_N}^{I_N} = %.03f,  \\beta_{S_N,I_N}^{I_V} = %.03f, \\beta_{S_N,I_N}^{I_D} = %.03f$$ \n 
                       $$ \\beta_{S_V,I_V}^{I_N} = %.03f,  \\beta_{S_V,I_V}^{I_V} = %.03f, \\beta_{S_V,I_V}^{I_D} = %.03f $$ \n
                       $$ \\beta_{S_D,I_D}^{I_N} = %.03f,  \\beta_{S_D,I_D}^{I_V} = %.03f, \\beta_{S_D,I_D}^{I_D} = %.03f $$', 
                       input$beta11, input$beta12, input$beta13, 
                       input$beta21, input$beta22, input$beta23, 
                       input$beta31, input$beta32, input$beta33 
               )
             )
      ),
      column(6, 
             withMathJax(
               sprintf('$$ \\nu_1 = %.04f, \\nu_2 = %.04f; \\quad $$ \n
                       $$ \\gamma = %.04f, \\kappa = %.04f; \\quad $$ \n
                       $$ \\alpha_1 =  %.04f,  \\alpha_2 =  %.04f $$', 
                       input$nu1, input$nu2, 
                       input$gamma, input$kappa, 
                       input$alpha1, input$alpha2 )
             )
      )
    )
  })
  
  
  output$SIR.ode <- renderUI({
    withMathJax(
      sprintf('\n$$\\frac{dS}{dt}   = -  \\beta  \\cdot I \\cdot S + \\kappa \\cdot  R$$ \n 
              $$\\frac{dI}{dt}   =  \\beta \\cdot I \\cdot  S - \\gamma \\cdot I $$ \n 
              $$\\frac{dR}{dt}   =  \\gamma_1 \\cdot I - \\kappa \\cdot R$$
              ')
      
      )
  })
  
  output$SR.ode <- renderUI({
    withMathJax(
      sprintf('\n$$\\frac{dS}{dt}   = - \\gamma \\cdot S + \\kappa \\cdot  R$$ \n 
              $$\\frac{dR}{dt}   =  \\gamma \\cdot S - \\kappa \\cdot R$$
              ')
      
      )
  })
  
  ## combine results into a single vector
  output$ext_current_ode_p1 <- renderUI({
    withMathJax(
      sprintf('$$\\frac{dS_N}{dt}   = - ( \\beta_{S_N,I_N}^{I_N}  \\cdot I_N + \\beta_{S_N,I_N}^{I_V} \\cdot I_V + \\beta_{S_N,I_N}^{I_D} \\cdot I_D ) \\cdot S_N - $$ \n 
              $$\\quad \\quad - \\nu_1 \\cdot S_N + \\kappa \\cdot  R_N$$ \n 
              $$\\frac{dS_V}{dt} = - ( \\beta_{S_V,I_V}^{I_N}  \\cdot I_N + \\beta_{S_V,I_V}^{I_V} \\cdot I_V + \\beta_{S_V,I_V}^{I_D} \\cdot I_D ) \\cdot S_V + $$ \n 
              $$\\quad \\quad + \\nu_1 \\cdot (1 - \\alpha) \\cdot S_N + \\nu_2 \\cdot V  - \\nu_1 \\cdot \\alpha \\cdot S_V  $$ \n 
              $$\\frac{dS_D}{dt} = - ( \\beta_{S_D,I_D}^{I_N}  \\cdot I_N + \\beta_{S_D,I_D}^{I_V} \\cdot I_V + \\beta_{S_D,I_D}^{I_D} \\cdot I_D) \\cdot S_D + $$ \n 
              $$\\quad \\quad + \\kappa \\cdot R_D $$ \n 
              $$\\quad$$
              $$\\frac{dV}{dt}   =  \\nu_1 \\cdot \\alpha \\cdot S_N + \\nu_1 \\cdot \\alpha \\cdot S_V - \\nu_2 \\cdot V $$  
              ')
      )
  })
  
  output$ext_current_ode_p2 <- renderUI({
    withMathJax(
      sprintf('$$\\frac{dI_N}{dt}   =  ( \\beta_{S_N,I_N}^{I_N}  \\cdot I_N + \\beta_{S_N,I_N}^{I_V} \\cdot I_V + \\beta_{S_N,I_N}^{I_D} \\cdot I_D ) \\cdot  S_N - \\gamma \\cdot I_N $$ \n 
              $$\\frac{dI_V}{dt} =  ( \\beta_{S_V,I_V}^{I_N}   \\cdot I_N + \\beta_{S_V,I_V}^{I_V} \\cdot I_V + \\beta_{S_V,I_V}^{I_D} \\cdot I_D ) \\cdot S_V - \\gamma \\cdot I_V$$ \n 
              $$\\frac{dI_D}{dt} =  ( \\beta_{S_D,I_D}^{I_N}   \\cdot I_N + \\beta_{S_D,I_D}^{I_V} \\cdot I_V + \\beta_{S_D,I_D}^{I_D} \\cdot I_D ) \\cdot S_D - \\gamma \\cdot I_D$$ \n 
              $$\\quad$$
              $$\\frac{dR_N}{dt}   =  \\gamma \\cdot I_N - \\kappa \\cdot R_N$$
              $$\\frac{dR_V}{dt}   =  \\gamma \\cdot I_V - \\kappa \\cdot R_V$$
              $$\\frac{dR_D}{dt}   =  \\gamma \\cdot I_D - \\kappa \\cdot R_D$$
              ')
      )
  })
  
  output$SIR.current.parameters <- renderUI({

    if(input$SIR.stage){
      withMathJax( sprintf('$$R^1_{0}= \\beta \\cdot \\gamma^{-1} = %.04f$$ \n $$R^2_{0}= \\beta_2 \\cdot \\gamma^{-1}_2 = %.04f$$', 
                         input$SIR.beta/input$SIR.gamma, input$SIR.beta2/input$SIR.gamma2))
    } else {
      withMathJax( sprintf('$$R^1_{0}= \\beta \\cdot \\gamma^{-1} = %.04f$$', 
                           input$SIR.beta/input$SIR.gamma))
    }
    
  })
  
  
  
  ode_curves_1 <- reactive({
    x_init_1  <- get_starting_values(input$dateRange[1]); myN <- 36e6
    x_init_1[1]  <- ( (x_init_1[1] + x_init_1[2] + x_init_1[3]) * myN - input$E - input$I)/myN
    x_init_1[2]  <- input$E/myN
    x_init_1[3]  <- input$I/myN
    my.params.final <- c(0,0.175,1/3.5,1/7, 
                         #ifelse(input$IsBeta, input$beta, 0.175), 
                         #ifelse(input$IsSigma, 1/input$sigma, 1/3.5),
                         #ifelse(input$IsGamma, 1/input$gamma, 1/7),
                         0.00267717, 0.01075589)
    
    estimate.curves.d.rec(my.params.final = my.params.final, my.est_params = NULL, my.constants = NULL, 
                          my.to_estimate = NULL, my.x_init = x_init_1, my.times = 1:365, 
                          my.start_date = input$dateRange[1], my.SEIR_model = SEIR.model.D.Rec)
  })
  
  output$sse_value <- renderUI({
    if(exists("ode_curves_1")){
      time_span <- seq(from = input$dateRange[1], to = input$dateRange[2], by=1)
      
      R_daily <- c(ode_curves_1()$R[1], diff(ode_curves_1()$R))[ ode_curves_1()$dates %in% time_span ]
      cases_daily <- c(wiki.tab_t0$cases[1], diff(wiki.tab_t0$cases))[wiki.tab_t0$dates %in% time_span]

      
      sse.value <- my.sse(R_daily[1:length(cases_daily)][-1], cases_daily[-1])
    } else {
      sse.value <- -1
    }
    withMathJax(
      sprintf('$$SSE = %s$$', format(round(as.numeric(sse.value), 1), nsmall=1, big.mark="\\,"))
    )
  })
  
  
  output$estimation_plot <- renderPlot({
    opt_calculation(E = input$Eest, I = input$Iest, 
                    start_date = input$dateRangeest[1], end_date = input$dateRangeest[2]) %...>% {
                      estimated_df <- .
                      getCurrentPlot(estimated_df, data_start = input$dateRangeest[1], data_end = input$dateRangeest[2])
                    }
  })
  
  
  output$simpleSIRplot <- renderPlot({
    
    times <- seq( from=0, to=input$SIR.t, by=1 )

    
    xstart <- c( S=input$SIR.S/sum(input$SIR.S + input$SIR.I + input$SIR.R),
                 I=input$SIR.I/sum(input$SIR.S + input$SIR.I + input$SIR.R),
                 R=input$SIR.R/sum(input$SIR.S + input$SIR.I + input$SIR.R) )

    
    parms1 <- c( SIR.beta  = input$SIR.beta,
                 SIR.gamma = input$SIR.gamma, 
                 SIR.kappa = input$SIR.kappa )

    
    #my.caption <- paste0("Scenario: ", 
    #                     "Go crazy: ", ifelse(input$go_crazy, "T", "F"), 
    #                     "; Vaccination: ", ifelse(input$fast_vacc, "fast", "slow"), 
    #                     "; Immunity: ", ifelse(input$variant, "200 days", "400 days"))
    #x    <- faithful[, 2] 
    #bins <- seq(min(x), max(x), length.out = input$bins + 1)
    
    #getCurrentPlot(ode_curves_1(), data_start = input$dateRange[1], data_end = input$dateRange[2])
    ode(
      func=simple.sir.model,
      y=xstart,
      times=times,
      parms=parms1
    ) %>% as.data.frame() -> out.pre
    

    if(input$SIR.stage){
      parms1 <- c( SIR.beta  = input$SIR.beta2, SIR.gamma = input$SIR.gamma2, SIR.kappa = input$SIR.kappa2 )
      xstart <- c( S=out.pre[out.pre$time == input$SIR.restriction, "S"],
                   I=out.pre[out.pre$time == input$SIR.restriction, "I"],
                   R=out.pre[out.pre$time == input$SIR.restriction, "R"] )

      times <- seq( from=input$SIR.restriction, to=input$SIR.t , by=1 )
      ode(
        func=simple.sir.model,
        y=xstart,
        times=times,
        parms=parms1
      ) %>% as.data.frame() -> out.pre.2
      out.pre.2 %>% melt(., id = "time") -> dat.2
      
    }
      
    out.pre %>% melt(., id = "time") -> out
    
    dat <- out
    
    dat %>% #  %>%
      ggplot(aes(x=time,y=value, color=variable)) +
      scale_color_manual("Populacje: ",
                         values = c("darkblue", "red2", "orange") #, "maroon", "orange", "yellow2", "darkgreen"),
      ) + 
  #    scale_linetype_manual("Efficiacy: ", values = c("solid", "dotted")) +
     scale_y_continuous(labels = scales::percent_format(accuracy = 0.5)) +
      theme_classic() +
      labs(x = "czas [dni]", y = "Procent społeczeństwa" ) +
      theme(text = element_text(size=18), 
            legend.text=element_text(size=18), 
            legend.position="top", 
            legend.box="vertical") -> g #+ 
    if(input$SIR.stage){
      g <- g + 
        geom_line(linetype = "dashed") +
        geom_line(data = dat.2, aes(x=time,y=value, color=variable)) + 
        geom_vline(xintercept = input$SIR.restriction, color = "orange2", linetype = "dotted", size = 1.1)
      
    } else {
      g <- g + geom_line(size=1, linetype = "solid") 
    }
    g
      #labs(x='time (days)', y='fraction of population', captiono = my.caption) 
  })
  
  output$simpleSRplot <- renderPlot({
    
    times <- seq( from=0, to=input$SR.t, by=1 )
    
    xstart <- c( S=input$SR.S/sum(input$SR.S + input$SR.R),
                 R=input$SR.R/sum(input$SR.S + input$SR.R) )
    
    parms1 <- c( SR.gamma  = input$SR.gamma,
                 SR.kappa = input$SR.kappa )
    
    ode(
      func=simple.sr.model,
      y=xstart,
      times=times,
      parms=parms1
    ) %>% as.data.frame() -> out.pre
    
    out.pre %>% melt(., id = "time") -> out
    
    dat <- out
    
    dat %>% #  %>%
      ggplot(aes(x=time,y=round(value*(input$SR.S + input$SR.R)), color=variable)) +
         scale_color_manual("Populacje: ",
                             values = c("darkblue", "red2") #, "maroon", "orange", "yellow2", "darkgreen"),
         ) + 
      #    scale_linetype_manual("Efficiacy: ", values = c("solid", "dotted")) +
      #scale_y_continuous(labels = scales::percent_format(accuracy = 0.5)) +
      theme_classic() +
     
      labs(x = "czas [dni]", y = "Liczba osobników" ) +
      theme(text = element_text(size=18), 
            legend.position="top", 
            legend.text=element_text(size=18), 
            legend.box="vertical") -> g #+ 
      g <- g + geom_line(size=1, linetype = "solid") 
    
      if(input$SR.solution){
        g <- g +  geom_line(data = data.frame(x = 0:(input$SR.t), 
                                              y = (input$SR.S + input$SR.R) * exp( - input$SR.gamma * (0:(input$SR.t)) + input$SR.S/(input$SR.S + input$SR.R) - 1)), 
                            aes(x, y), 
                            color = "red3", linetype="dotted", size = 2) 
      }
    g
    #labs(x='time (days)', y='fraction of population', captiono = my.caption) 
  })
  
  output$distPlot <- renderPlot({
    # generate bins based on input$bins from ui.R
    
    times <- seq(from=0,to=2*365,by=1/2)
    
    xstart <- c(S=input$S/sum(input$S + input$I + input$R),Sv=0,
                V=0,
                I=input$I/sum(input$S + input$I + input$R),Iv=0,
                R=input$R/sum(input$S + input$I + input$R))
    gamma <- 1/6
    R0    <- 2.5
    beta  <- 0.25
    parms1 <- c(beta11=input$beta11,beta12=input$beta12,beta21=input$beta21,beta22=input$beta22,
                gamma_1=input$gamma, gamma_2=input$gamma,
                nu_1=input$nu1,nu_2=input$nu2,
                kappa=input$kappa,
                alpha=input$alpha1)
    
    parms2<- c(beta11=input$beta11,beta12=input$beta12,beta21=input$beta21,beta22=input$beta22,
               gamma_1=input$gamma, gamma_2=input$gamma,
               nu_1=input$nu1,nu_2=input$nu2,
               kappa=input$kappa,
               alpha=input$alpha2)
    
    my.caption <- paste0("Scenario: ", 
                         "Go crazy: ", ifelse(input$go_crazy, "T", "F"), 
                         "; Vaccination: ", ifelse(input$fast_vacc, "fast", "slow"), 
                         "; Immunity: ", ifelse(input$variant, "200 days", "400 days"))
    #x    <- faithful[, 2] 
    #bins <- seq(min(x), max(x), length.out = input$bins + 1)
    
    #getCurrentPlot(ode_curves_1(), data_start = input$dateRange[1], data_end = input$dateRange[2])
    ode(
      func=closed.sir.model,
      y=xstart,
      times=times,
      parms=parms1
    ) %>%
      as.data.frame() -> out
    
    ode(
      func=closed.sir.model,
      y=xstart,
      times=times,
      parms=parms2
    ) %>%
      as.data.frame() -> out2
    # draw the histogram with the specified number of bins
    #hist(x, breaks = bins, col = 'darkgray', border = 'white')
    
    
    ## usuwamy duze kompartmenty
    ## Zostaje I oraz Iv
    head(out2)
    dat <- rbind(
      cbind(out[, c("time", "I", "Iv")], alpha = input$alpha1), 
      cbind(out2[, c("time", "I", "Iv")], alpha = input$alpha2)) %>%
      #.[,c(("time", "I", "Iv")] %>%
      gather(variable,value,-c(time,alpha) ) 
    
    
    dat$alpha <- factor(dat$alpha)
    levels(dat$alpha) <-  paste0(round(as.numeric(levels(dat$alpha))*100), "%")
    
    
    dat$variable <- factor(dat$variable, labels = c("Not vaccinated (I)", "After vaccination (Iv)"))
    
    dat %>% #  %>%
      ggplot(aes(x=time,y=value, color=variable, group = interaction(alpha, variable), linetype = alpha)) +
      scale_color_manual("Infections: ",
                         values = c("darkblue", "red2", "orange") #, "maroon", "orange", "yellow2", "darkgreen"),
      ) + 
      scale_linetype_manual("Efficacy: ", values = c("solid", "dotted")) +
      scale_y_continuous(labels = scales::percent_format(accuracy = 0.5)) +
      geom_line( size=1) +
      theme_classic()+
      theme(text = element_text(size=18), 
            legend.text=element_text(size=18), 
            legend.position="top", 
            legend.box="vertical") + 
      labs(x='time (days)', y='fraction of population')
  })
  
  output$extDistPlot <- renderPlot({
    # generate bins based on input$bins from ui.R
    
    times <- seq(from=0,to=2*365,by=1)
    
    xstart <- c(S  = (1 - input$sceptics) * input$S/sum(input$S + input$I + input$R), 
                Sv = 0, 
                Sd = (input$sceptics) * input$S/sum(input$S + input$I + input$R),
                V  = 0,
                I  = (1 - input$sceptics) * input$I/sum(input$S + input$I + input$R), 
                Iv = 0, 
                Id = (input$sceptics) * input$I/sum(input$S + input$I + input$R),
                R  = (1 - input$sceptics) * input$R/sum(input$S + input$I + input$R), 
                Rv = 0,
                Rd = (input$sceptics) * input$R/sum(input$S + input$I + input$R) )
    #gamma <- 1/6
    #R0    <- 2.5
    #  beta  <- 0.25
    parms1 <- c(beta11=input$beta11,beta12=input$beta12,beta13=input$beta13,
                beta21=input$beta21,beta22=input$beta22,beta23=input$beta23,
                beta31=input$beta31,beta32=input$beta32,beta33=input$beta33,
                gamma_1=input$gamma, gamma_2=input$gamma,
                nu_1=input$nu1,nu_2=input$nu2,
                kappa=input$kappa, proportional_mixing = input$proportional_mixing,
                alpha=input$alpha1)
    
    parms2<- c(beta11=input$beta11,beta12=input$beta12,beta13=input$beta13,
               beta21=input$beta21,beta22=input$beta22,beta23=input$beta23,
               beta31=input$beta31,beta32=input$beta32,beta33=input$beta33,
               gamma_1=input$gamma, gamma_2=input$gamma,
               nu_1=input$nu1,nu_2=input$nu2,
               kappa=input$kappa, proportional_mixing = input$proportional_mixing,
               alpha=input$alpha2)
    
    my.caption <- paste0("Scenario: ", 
                         "Go crazy: ", ifelse(input$go_crazy, "T", "F"), 
                         "; Vaccination: ", ifelse(input$fast_vacc, "fast", "slow"), 
                         "; Immunity: ", ifelse(input$variant, "200 days", "400 days"))
    #x    <- faithful[, 2] 
    #bins <- seq(min(x), max(x), length.out = input$bins + 1)
    
    #getCurrentPlot(ode_curves_1(), data_start = input$dateRange[1], data_end = input$dateRange[2])
    ode(
      func=ext.closed.sir.model,
      y=xstart,
      times=times,
      parms=parms1
    ) %>%
      as.data.frame() -> out
    
    ode(
      func=ext.closed.sir.model,
      y=xstart,
      times=times,
      parms=parms2
    ) %>%
      as.data.frame() -> out2
    # draw the histogram with the specified number of bins
    #hist(x, breaks = bins, col = 'darkgray', border = 'white')

    
    ## usuwamy duze kompartmenty
    ## Zostaje I oraz Iv
    
    dat <- rbind(
      cbind(out[, c("time", "I", "Iv", "Id")], alpha = input$alpha1), 
      cbind(out2[, c("time", "I", "Iv", "Id")], alpha = input$alpha2)) %>%
      #.[,c(("time", "I", "Iv")] %>%
      gather(variable,value,-c(time,alpha) ) 
    
    dat$alpha <- factor(dat$alpha)
    my.labels <- paste0(as.double(levels(dat$alpha)) * 100, "%")
    levels(dat$alpha) <- my.labels
    
    dat$variable <- factor(dat$variable)
    levels(dat$variable) <- c("Not vaccinated (In)", "Anti-vaccination (Id)", "After vaccination(Iv)")
    In.max <- 1e6*out[out$time==input$lookup,]$I ; In.max.2 <- 1e6*out2[out2$time==input$lookup,]$I 
    In.max <- ifelse (In.max < 1, round(In.max, 3), round(In.max))
    In.max.2 <- ifelse (In.max.2 < 1, round(In.max.2, 3), round(In.max.2))
    
    Iv.max <- 1e6*out[out$time==input$lookup,]$Iv; Iv.max.2 <- 1e6*out2[out2$time==input$lookup,]$Iv
    Iv.max <- ifelse (Iv.max < 1, round(Iv.max, 3), round(Iv.max))
    Iv.max.2 <- ifelse (Iv.max.2 < 1, round(Iv.max.2, 3), round(Iv.max.2))
    
    Id.max <- 1e6*out[out$time==input$lookup,]$Id; Id.max.2 <- 1e6*out2[out2$time==input$lookup,]$Id 
    Id.max <- ifelse (Id.max < 1, round(Id.max, 3), round(Id.max))
    Id.max.2 <- ifelse (Id.max.2 < 1, round(Id.max.2, 3), round(Id.max.2))
    #dat$value[dat$value == 0] <- NA
    if(length(levels(dat$alpha)) > 1){
      g <- dat %>% #  %>%
      ggplot(aes(x=time,y=value * 1e6, color=variable, group = interaction(alpha, variable), linetype = alpha)) +
      scale_linetype_manual("Efficacy: ", values = c("solid", "dotted")) 
    } else {
      g <- dat %>% #  %>%
        ggplot(aes(x=time,y=value * 1e6, color=variable, group = interaction(alpha, variable))) 
    }
      g <- g + scale_color_manual("Infections: ",
                         values = c("darkblue", "red2", "orange") #, "maroon", "orange", "yellow2", "darkgreen"),
      ) + 
      #scale_y_continuous(labels = scales::percent_format(accuracy = 0.5)) +
      geom_line(size = 1) + 
      geom_vline(xintercept = input$lookup) + 
      annotate("label", 
               alpha = 0.5, 
               x=input$lookup, 
               y=max(In.max, Iv.max, Id.max), 
               label= paste0("In = ", In.max , " | ", In.max.2, 
                             "\nIv = ", Iv.max, " | ", Iv.max.2, 
                             "\nId = ",Id.max, " | ", Id.max.2), 
               size = 5, hjust = "inward", vjust = "inward") + 
      
      
      #geom_text(data = data.frame(x=400, y = Inf, lab = paste0("I=",out[out$time==400*2,]$I,
      #                  "I_V=",out[out$time==400*2,]$Iv,
      #                  "I_D=",out[out$time==400*2,]$Id )), mapping = aes(x, y, label=lab)) + 
      theme_classic()+
      theme(text = element_text(size=18), 
            legend.text=element_text(size=18), 
            legend.position="top", 
            legend.box="vertical") + 
      labs(x='time (days)', 
           y = "Infected per milion"
           #y='fraction of population', 
           #caption = my.caption
           )
    
     if(input$log_scale){
       g <- g + scale_y_log10( breaks = trans_breaks("log10", function(x) 10^x),
                              labels = trans_format("log10", math_format(10^.x))) 
     }
    
     g
  })
  
  output$extDistPlotv2 <- renderPlot({
    
    plots <- plotInput()
    
    ggarrange(plots[["g.cum"]], plots[["g.daily"]], plots[["g.prop.mix"]],
              plots[["g.muller"]],  nrow = 2, ncol = 2) 
  })
  
  plotInputDownload <- reactive({
    plots <- plotInput()
    ggarrange(plots[["g.daily"]], plots[["g.muller"]], nrow = 2, ncol = 1) 
  })
  
  plotInput <- reactive({
    # generate bins based on input$bins from ui.R
    
    times <- seq(from=0,to=2*365,by=1)
    
    Npop <- sum(input$S + input$I + input$R + input$Sv1 + input$V)
    xstart <- c(S  = (1 - input$sceptics) * input$S/Npop, 
                Sv1 = input$Sv1/Npop, Sv2 = 0,
                Sd = (input$sceptics) * input$S/Npop,
                V  = input$V/Npop,
                I  = (1 - input$sceptics) * input$I/Npop, 
                Iv1 = 0, Iv2 = 0, 
                Id = (input$sceptics) * input$I/Npop,
                R  = (1 - input$sceptics) * input$R/Npop, 
                Rv = 0,
                Rd = (input$sceptics) * input$R/Npop )
    #gamma <- 1/6
    #R0    <- 2.5
    #  beta  <- 0.25
    parms1 <- c(beta11=input$beta11,beta12=input$beta12,beta13=input$beta13,
                beta21=input$beta21,beta22=input$beta22,beta23=input$beta23,
                beta31=input$beta31,beta32=input$beta32,beta33=input$beta33,
                gamma_1=input$gamma, gamma_2=input$gamma,
                nu_1=input$nu1,nu_2=input$nu2,nu_1star=input$nu1star,
                kappa=input$kappa,a=input$sceptics, proportional_mixing = input$proportional_mixing,
                alpha=input$alpha1)
    
    parms2<- parms1
    parms2["alpha"] = input$alpha2
               #c(beta11=input$beta11,beta12=input$beta12,beta13=input$beta13,
               #beta21=input$beta21,beta22=input$beta22,beta23=input$beta23,
               #beta31=input$beta31,beta32=input$beta32,beta33=input$beta33,
               #gamma_1=input$gamma, gamma_2=input$gamma,
               #nu_1=input$nu1,nu_2=input$nu2,nu_1star=input$nu1star,
               #kappa=input$kappa,a=input$sceptics, proportional_mixing = input$proportional_mixing,
               #alpha=input$alpha2)
    
    parms3<- parms1
    parms3["proportional_mixing"] = T
    
    my.caption <- paste0("Scenario: ", 
                         "Go crazy: ", ifelse(input$go_crazy, "T", "F"), 
                         "; Vaccination: ", ifelse(input$fast_vacc, "fast", "slow"), 
                         "; Immunity: ", ifelse(input$variant, "200 days", "400 days"))
    #x    <- faithful[, 2] 
    #bins <- seq(min(x), max(x), length.out = input$bins + 1)
    
    #getCurrentPlot(ode_curves_1(), data_start = input$dateRange[1], data_end = input$dateRange[2])
    ode(
      func=ext.closed.sir.model.v2,
      y=xstart,
      times=times,
      parms=parms1
    ) %>%
      as.data.frame() -> out
    
    muller.out <- out
    muller.out$S <- muller.out$S + muller.out$Sd
    muller.out$I <- muller.out$I + muller.out$Id + muller.out$Iv1 + muller.out$Iv2
    muller.out$R <- muller.out$R + muller.out$Rd + muller.out$Rv
    muller.out <- muller.out[, c("time", "S", "Sv1", "Sv2", "V", "I", "R")]
    
    muller.out %>% gather(variable,value,-c(time) ) -> melted.out
    colnames(melted.out) <- c("Generation", "Group_id", "Frequency")
    melted.out$Unique_id <- paste0(melted.out$Group_id, "_", melted.out$Generation)
    melted.out$Identity <- factor(melted.out$Group_id, levels = c("S", "Sv1", "Sv2", "V", "I", "R"))
    melted.out$Population <- melted.out$Frequency * 1e6
    
    
    ode(
      func=ext.closed.sir.model.v2,
      y=xstart,
      times=times,
      parms=parms2
    ) %>%
      as.data.frame() -> out2
    
    ode(
      func=ext.closed.sir.model.v2,
      y=xstart,
      times=times,
      parms=parms3
    ) %>%
      as.data.frame() -> out3
    # draw the histogram with the specified number of bins
    #hist(x, breaks = bins, col = 'darkgray', border = 'white')
    
    ## usuwamy duze kompartmenty
    ## Zostaje I oraz Iv
    
    dat <- rbind(
      cbind(out[, c("time", "Iv1", "Iv2")], I = rowSums(out[, c("I", "Id")]), Isum = rowSums(out[, c("I", "Iv1", "Iv2", "Id")]), alpha = input$alpha1), 
      cbind(out2[, c("time", "Iv1", "Iv2")],I = rowSums(out2[, c("I", "Id")]), Isum = rowSums(out2[, c("I", "Iv1", "Iv2", "Id")]), alpha = input$alpha2)) %>%
      
      #cbind(out[, c("time", "I", "Iv1", "Iv2", "Id")], Isum = rowSums(out[, c("I", "Iv1", "Iv2", "Id")]), alpha = input$alpha1), 
      #cbind(out2[, c("time", "I", "Iv1", "Iv2","Id")], Isum = rowSums(out2[, c("I", "Iv1", "Iv2", "Id")]), alpha = input$alpha2)) %>%
      #.[,c(("time", "I", "Iv")] %>%
      gather(variable,value,-c(time,alpha) ) 
    
    dat3 <- 
      cbind(out3[, c("time", "Iv1", "Iv2")],
            I = rowSums(out3[, c("I", "Id")]), 
            Isum = rowSums(out3[, c("I", "Iv1", "Iv2", "Id")]), 
            alpha = parms3["alpha"]) %>%
      gather(variable,value,-c(time,alpha) ) 
    
    dat3$variable <- factor(dat3$variable, c("I", "Iv1", "Iv2", "Isum"))
    print(levels(dat3$variable))
    levels(dat3$variable) <- c("Not vaccinated (I)", "After vaccination (Iv1)", "After vaccination (Iv2)", "All")
    
    
    tmp.dat <- cbind(out[, c("Iv1", "Iv2")], 
                     I = rowSums(out[, c("I", "Id")]), 
                     Isum = rowSums(out[, c("I", "Iv1", "Iv2", "Id")]))
    
    tmp.dat2 <- cbind(out2[, c("Iv1", "Iv2")], 
                     I = rowSums(out2[, c("I", "Id")]), 
                     Isum = rowSums(out2[, c("I", "Iv1", "Iv2", "Id")]))
    
    dat.diff <- rbind(
      cbind(time = out[-1, c("time")], 
        alpha = input$alpha1, 
        data.frame(diff(as.matrix(
          tmp.dat 
          ))) + (1/6) * tmp.dat[-1,] 
        ),
      cbind(time = out2[-1, c("time")], 
            alpha = input$alpha2, 
            data.frame(diff(as.matrix(
              tmp.dat2 
            ))) + (1/6) * tmp.dat2[-1,] 
      )) %>%
      gather(variable,value,-c(time,alpha) )
    
    print(head(dat.diff))
    
    dat.diff$alpha <- factor(dat.diff$alpha)
    my.labels <- paste0(as.double(levels(dat.diff$alpha)) * 100, "%")
    levels(dat.diff$alpha) <- my.labels
    dat.diff$variable <- factor(dat.diff$variable, levels = c("I", "Iv1", "Iv2", "Isum"))
    print(levels(dat.diff$variable))
    
    
    dat$alpha <- factor(dat$alpha)
    my.labels <- paste0(as.double(levels(dat$alpha)) * 100, "%")
    levels(dat$alpha) <- my.labels
    dat$variable <- factor(dat$variable, levels = c("I", "Iv1", "Iv2", "Isum"))
    print(levels(dat$variable))
    #levels(dat$variable) <- c("Not vaccinated (In)", "Anti-vaccination (Id)", "Total infected", "After vaccination(Iv1)", "After vaccination(Iv2)")
    levels(dat$variable) <- c("Not vaccinated (I)", "After vaccination (Iv1)", "After vaccination (Iv2)", "All")
    
    In.max <- 1e6*out[out$time==input$lookup,]$I ; In.max.2 <- 1e6*out2[out2$time==input$lookup,]$I 
    In.max <- ifelse (In.max < 1, round(In.max, 3), round(In.max))
    In.max.2 <- ifelse (In.max.2 < 1, round(In.max.2, 3), round(In.max.2))
    
    Iv1.max <- 1e6*out[out$time==input$lookup,]$Iv1; Iv1.max.2 <- 1e6*out2[out2$time==input$lookup,]$Iv1
    Iv1.max <- ifelse (Iv1.max < 1, round(Iv1.max, 3), round(Iv1.max))
    Iv1.max.2 <- ifelse (Iv1.max.2 < 1, round(Iv1.max.2, 3), round(Iv1.max.2))
    Iv2.max <- 1e6*out[out$time==input$lookup,]$Iv2; Iv2.max.2 <- 1e6*out2[out2$time==input$lookup,]$Iv2
    Iv2.max <- ifelse (Iv2.max < 1, round(Iv2.max, 3), round(Iv2.max))
    Iv2.max.2 <- ifelse (Iv2.max.2 < 1, round(Iv2.max.2, 3), round(Iv2.max.2))
    
    #Id.max <- 1e6*out[out$time==input$lookup,]$Id; Id.max.2 <- 1e6*out2[out2$time==input$lookup,]$Id 
    #Id.max <- ifelse (Id.max < 1, round(Id.max, 3), round(Id.max))
    #Id.max.2 <- ifelse (Id.max.2 < 1, round(Id.max.2, 3), round(Id.max.2))
    
    Isum.max <- 1e6*rowSums(out[, c("I", "Iv1", "Iv2")])[input$lookup]; 
    Isum.max.2 <- 1e6*rowSums(out2[, c("I", "Iv1", "Iv2")])[input$lookup] 
    Isum.max <- ifelse (Isum.max < 1, round(Isum.max, 3), round(Isum.max))
    Isum.max.2 <- ifelse (Isum.max.2 < 1, round(Isum.max.2, 3), round(Isum.max.2))
    

    mytable <- cbind(sites=c("Isum","I","Iv1","Iv2"),
                     data.frame(a=c(Isum.max, In.max, Iv1.max, Iv2.max), 
                                b=c(Isum.max.2, In.max.2, Iv1.max.2, Iv2.max.2)))

    names(mytable) <- c("group", my.labels)
    rownames(mytable) <- NULL
    
    my.label.coords <- get.coords.for.label(input$lookup, 
                                  xmin = min(dat$time), xmax = max(dat$time),
                                  ymin = 0, ymax = max(dat$value) * 1e6)
    
    #dat$value[dat$value == 0] <- NA
    
    if(length(levels(dat$alpha)) > 1){
      g <- dat %>% #  %>%
        ggplot(aes(x=time,y=value * 1e6, color=variable, group = interaction(alpha, variable), linetype = alpha)) +
        scale_linetype_manual("Efficacy: ", values = c("solid", "dotted")) 
    } else {
      g <- dat %>% #  %>%
        ggplot(aes(x=time,y=value * 1e6, color=variable, group = interaction(alpha, variable))) 
    }
    g <- g + 
      scale_color_manual("Group: ",
                         labels = expression(I, I[1], I[2], "All"),
                         values = c("#e05f38", "#e08138", "#e0ae38", "#EE6D56", "darkblue") #, "maroon", "orange", "yellow2", "darkgreen"),
      ) + 
      #scale_linetype_manual("Efficacy: ", values = c("solid", "dotted")) +
      #scale_y_continuous(labels = scales::percent_format(accuracy = 0.5)) +
      geom_line(size = 1) + 
      #### Annotation table with values
      #geom_vline(xintercept = input$lookup) + 
      #annotation_custom(tableGrob(mytable, rows = NULL, 
      #                            theme = ttheme_default(
      #                              core=list(bg_params = list(fill = "lightgrey", col=NA, alpha = 0.8),
      #                                        fg_params=list(fontface=3)))), 
      #                  xmin = my.label.coords$xmin, 
      #                  xmax = my.label.coords$xmax,
      #                  ymax = my.label.coords$ymax, 
      #                  ymin = my.label.coords$ymin) +
      theme_minimal() +
      guides(color = guide_legend(ncol = 1, override.aes = list(size = 2))) + 
      theme(text = element_text(size=18),
            legend.text=element_text(size=18),
            legend.key.width = unit(1,"cm"), 
            legend.position="right", legend.box="vertical") + 
      labs(x='Time (days)', 
           y = "Infected in 1M pop.",
           caption = "Proportional mixing"
           #y='fraction of population', 
           #caption = my.caption
           ) 
    
    if(length(levels(dat.diff$alpha)) > 1){
      g.diff <- dat.diff  %>% #  %>%
        ggplot(aes(x=time, y=value * 1e6, 
                   color=variable, 
                   group = interaction(alpha, variable), 
                   linetype = alpha)) +
        scale_linetype_manual("Efficacy: ", values = c("solid", "dotted")) 
    } else {
      g.diff <- dat.diff  %>% #  %>%
        ggplot(aes(x=time,y=value * 1e6, color=variable, group = interaction(alpha, variable))) 
    }
    g.diff <- g.diff + #dat.diff %>% #  %>%
      #ggplot(aes(x=time,y=value * 1e6, color=variable, group = interaction(alpha, variable), linetype = alpha)) +
      scale_color_manual("Group:",
                         labels = expression(I, I[1], I[2], "All"),
                         values = c("#dee500", "#ff458f", "#5e3c99", "#ff3747", "darkblue") #, "maroon", "orange", "yellow2", "darkgreen"),
      ) + 
      #scale_linetype_manual("Efficacy: ", values = c("solid", "dotted")) +
      #scale_y_continuous(labels = scales::percent_format(accuracy = 0.5)) +
      geom_line(size = 1, alpha = 0.8, linetype = "dashed") +
      geom_line(size = 1, data = subset(dat.diff, variable == "Isum"), linetype = "solid", color = "#ff3747") +
      #geom_bar(stat = "identity", position=position_dodge()) +
      #geom_vline(xintercept = input$lookup) + 
      #annotation_custom(tableGrob(mytable, rows = NULL, 
      #                            theme = ttheme_default(
      #                              core=list(bg_params = list(fill = "lightgrey", col=NA, alpha = 0.8),
      #                                        fg_params=list(fontface=3)))), 
      #                  xmin = my.label.coords$xmin, 
      #                  xmax = my.label.coords$xmax,
      #                  ymax = my.label.coords$ymax, 
      #                  ymin = my.label.coords$ymin) +
      theme_minimal() +
      guides(color = guide_legend(ncol = 1, override.aes = list(size = 2))) + 
      theme(text = element_text(size=18), 
            legend.text=element_text(size=18),
            legend.key.width = unit(1,"cm"), 
            legend.position="right", 
            legend.box="vertical") + 
      labs(x='Time (days)', 
           y = "Daily incidents per 1M pop."
           #caption = "Proportional mixing"
           #y='fraction of population', 
           #caption = my.caption
      )
    if(input$log_scale){
      g.diff <- g.diff + 
        scale_y_log10( breaks = trans_breaks("log10", function(x) 10^x),
                       labels = trans_format("log10", math_format(10^.x)), 
                       limits = c(10^-10, 10^5)) + 
        annotation_logticks(sides="l") +
        labs( y = "Daily incidents per 1M pop. [log]")
    } 
   
    
    g.prop.mix <- dat3 %>% 
      ggplot(aes(x=time,y=value * 1e6, color=variable, group = variable)) +
      scale_color_manual("Group:",
                         labels = expression(I, I[1], I[2], "All"),
                         values = c("#e05f38", "#e08138", "#e0ae38", "#EE6D56", "darkblue") #, "maroon", "orange", "yellow2", "darkgreen"),
      ) + 
      #scale_linetype_manual("Efficiacy: ", values = c("solid", "dotted")) +
      #scale_y_continuous(labels = scales::percent_format(accuracy = 0.5)) +
      geom_line(size = 1) + 
      #### Annotation table with values
      #geom_vline(xintercept = input$lookup) + 
      #annotation_custom(tableGrob(mytable, rows = NULL, 
      #                            theme = ttheme_default(
      #                              core=list(bg_params = list(fill = "lightgrey", col=NA, alpha = 0.8),
      #                                        fg_params=list(fontface=3)))), 
      #                  xmin = my.label.coords$xmin, 
      #                  xmax = my.label.coords$xmax,
      #                  ymax = my.label.coords$ymax, 
      #                  ymin = my.label.coords$ymin) +
      theme_minimal() +
      guides(color = guide_legend(ncol = 1, override.aes = list(size = 2))) + 
      theme(text = element_text(size=18), 
            legend.text=element_text(size=18),
            legend.key.width = unit(1,"cm"), 
            legend.position="right", 
            legend.box="vertical") + 
      labs(x='Time (days)', 
           y = "Infected in 1M pop.",
           #y='fraction of population', 
           caption = "Preferential mixing"
      )
    
    if(input$log_scale){
      g <- g + 
        scale_y_log10( breaks = trans_breaks("log10", function(x) 10^x),
                       labels = trans_format("log10", math_format(10^.x)), 
                       limits = c(10^-10, 10^5)) + 
        annotation_logticks(sides="l") +
        labs( y = "Infected in 1M pop. [log]")
    } 
    #else if (input$per100k) {
    #  g <- g 
    #} 
    
   g.muller <- Muller_plot(melted.out, 
                             add_legend = TRUE) + 
                theme_minimal() + 
                guides(fill = guide_legend(ncol = 1)) + 
                scale_y_continuous(labels = scales::percent_format(accuracy = 1)) + 
                scale_fill_manual("Group:", values = c("#FFF97E", "#EFBB40", "#F29837", "#81D452", "#EE6D56", "#79BFF9"), 
                                  labels = expression(S, S[1], S[2], V, I, R)) +
                theme(text = element_text(size=18),
                      axis.title.x =  element_blank(),
                      legend.text=element_text(size=18),
                      legend.position="right", 
                      legend.box="vertical") + 
                labs(x = "Time (days)", 
                     y = "Fraction of population")  
   list(g.cum = g, g.muller = g.muller, g.daily = g.diff, g.prop.mix = g.prop.mix)
   
                
  })
  
  output$tbl = renderDT({
    data_selected <- ode_curves_1()[ode_curves_1()$dates %in% seq(from = input$dateRange[1], to = input$dateRange[2], by=1), ]
    data_selected[,c(2:7)] <- round(data_selected[,c(2:7)])
    tbl <- datatable(
      data_selected, options = list(lengthChange = FALSE)
    )
    tbl
  })
  
  output$downloadData <- downloadHandler(
    filename = function() {
      paste(paste0("V-SIR_data_", Sys.Date()), ".csv", sep = "")
    },
    content = function(file) {
      write.csv(datasetInput(), file, row.names = FALSE)
    }
  )
  
  output$downloadPlot <- downloadHandler(
    filename = paste(paste0("V-SIR_plot_", Sys.Date()), ".png", sep = ""),
    content = function(file) {
      device <- function(..., width, height) {
        grDevices::png(..., width = width, height = height, pointsize = 9, 
                       res = 300, units = "in")
      }
      ggsave(file, plot = plotInputDownload(), device = device, width = 8, height = 8, units = "in")
  })
  
  datasetInput <- reactive({
    # generate bins based on input$bins from ui.R
    
    times <- seq(from=0,to=2*365,by=1)
    
    xstart <- c(S  = (1 - input$sceptics) * input$S/sum(input$S + input$I + input$R), 
                Sv1 = 0, Sv2 = 0,
                Sd = (input$sceptics) * input$S/sum(input$S + input$I + input$R),
                V  = 0,
                I  = (1 - input$sceptics) * input$I/sum(input$S + input$I + input$R), 
                Iv1 = 0, Iv2 = 0, 
                Id = (input$sceptics) * input$I/sum(input$S + input$I + input$R),
                R  = (1 - input$sceptics) * input$R/sum(input$S + input$I + input$R), 
                Rv = 0,
                Rd = (input$sceptics) * input$R/sum(input$S + input$I + input$R) )
    
    parms1 <- c(beta11=input$beta11,beta12=input$beta12,beta13=input$beta13,
                beta21=input$beta21,beta22=input$beta22,beta23=input$beta23,
                beta31=input$beta31,beta32=input$beta32,beta33=input$beta33,
                gamma_1=input$gamma, gamma_2=input$gamma,
                nu_1=input$nu1,nu_2=input$nu2,nu_1star=input$nu1star,
                kappa=input$kappa,a=input$sceptics,
                alpha=input$alpha1)
    
    ode(
      func=ext.closed.sir.model.v2,
      y=xstart,
      times=times,
      parms=parms1
    ) %>%
      as.data.frame() -> out
    
    out[,c("time", "S", "Sv1", "Sv2", "V")]
    
  })
  
  
  }

simple.sir.model <- function (t, x, params) {
  ## first extract the state variables
  N <- sum(x[1:3])
  S <- x[1]/N
  I <- x[2]/N
  R <- x[3]/N
  
  SIR.beta  <- params["SIR.beta"];  
  SIR.gamma <- params["SIR.gamma"];
  SIR.kappa <- params["SIR.kappa"]
  
  ## now code the model equations
  dSdt  <- - SIR.beta * S * I  + SIR.kappa * R
  dIdt  <-   SIR.beta * I * S  - SIR.gamma * I
  dRdt  <-   SIR.gamma * I     - SIR.kappa * R
  
  ## combine results into a single vector
  dxdt  <- c(dSdt, dIdt, dRdt)
  
  ## return result as a list!
  list(dxdt)
}

simple.sr.model <- function (t, x, params) {
  ## first extract the state variables
  N <- sum(x[1:2])
  S <- x[1]/N
  R <- x[2]/N
  
  SR.gamma  <- params["SR.gamma"];  
  SR.kappa <- params["SR.kappa"];
  
  ## now code the model equations
  dSdt  <- - SR.gamma * S  + SR.kappa * R
  dRdt  <-   SR.gamma * S  - SR.kappa * R
  
  ## combine results into a single vector
  dxdt  <- c(dSdt, dRdt)
  
  ## return result as a list!
  list(dxdt)
}


closed.sir.model <- function (t, x, params) {
  ## first extract the state variables
  N <- sum(x[1:6])
  S <- x[1]/N
  Sv <- x[2]/N
  V <- x[3]/N
  I <- x[4]/N
  Iv <- x[5]/N
  R <- x[6]/N
  
  ## now extract the parameters
  beta11 <- params["beta11"]; beta12 <- params["beta12"]; beta21 <- params["beta21"]; beta22 <- params["beta22"]; 
  gamma_1 <- params["gamma_1"]; gamma_2 <- params["gamma_2"];
  nu_1 <- params["nu_1"]; nu_2 <- params["nu_2"]
  kappa <- params["kappa"]
  alpha <- params["alpha"]
  
  ## now code the model equations
  dSdt  <- - beta11 * S  * I - beta12 * S  * Iv - nu_1 * S + kappa * R
  dSvdt <- - beta21 * Sv * I - beta22 * Sv * Iv + nu_2 * V + nu_1 * (1-alpha) * S
  
  dVdt  <- nu_1 * alpha * S - nu_2 * V
  
  dIdt  <- (beta11 * I + beta12 * Iv) * S  - gamma_1 * I
  dIvdt <- (beta21 * I + beta22 * Iv) * Sv - gamma_2 * Iv
  
  dRdt <- gamma_1 * I + gamma_2 * Iv - kappa * R
  ## combine results into a single vector
  dxdt <- c(dSdt,dSvdt,dVdt,dIdt,dIvdt,dRdt)
  ## return result as a list!
  list(dxdt)
}

#EXTENDED V-SIR MODEL
ext.closed.sir.model <- function (t, x, params) {
  ## first extract the state variables
  N   <- sum(x[1:10])
  S   <- x['S']/N
  Sv  <- x['Sv']/N
  Sd  <- x['Sd']/N
  V   <- x['V']/N
  I   <- x['I']/N
  Iv  <- x['Iv']/N
  Id  <- x['Id']/N
  R   <- x['R']/N
  Rv  <- x['Rv']/N
  Rd  <- x['Rd']/N
  
  ## now extract the parameters
  beta11 <- params["beta11"]; beta12 <- params["beta12"]; beta13 <- params["beta13"]; 
  beta21 <- params["beta21"]; beta22 <- params["beta22"]; beta23 <- params["beta23"]; 
  beta31 <- params["beta31"]; beta32 <- params["beta32"]; beta33 <- params["beta33"]; 
  gamma_1 <- params["gamma_1"]; gamma_2 <- params["gamma_2"];
  nu_1 <- params["nu_1"]; nu_2 <- params["nu_2"]
  kappa <- params["kappa"]
  alpha <- params["alpha"]
  
  ## now code the model equations
  dSdt  <-  - beta11 * S  * I - beta12 * S  * Iv - beta13 * S  * Id - nu_1 * S + kappa * R
  dSvdt <-  - beta21 * Sv * I - beta22 * Sv * Iv - beta23 * Sv * Id + ## all to Infected after vaccination (Iv)
    nu_2 * V + ## lost immunity after vaccination
    nu_1 * (1-alpha) * S - ## immunity not gained after first vaccination
    nu_1 * Sv + ## re-vaccination
    nu_1 * (1-alpha) * Sv ## simmunity not gained after re-vaccination
  dSddt <-  - beta31 * Sd * I - beta32 * Sd * Iv - beta33 * Sd * Id + kappa * Rd
  
  dVdt  <-  nu_1 * alpha * S - ## immune after first vaccination
    nu_2 * V + ## lost immunity 
    nu_1 * alpha * Sv ## immune after re-vaccination
  
  dIdt  <- (beta11 * I + beta12 * Iv + beta13 * Id) * S  - gamma_1 * I
  dIvdt <- (beta21 * I + beta22 * Iv + beta23 * Id) * Sv - gamma_1 * Iv
  dIddt <- (beta31 * I + beta32 * Iv + beta33 * Id) * Sd - gamma_1 * Id
  
  dRdt  <- gamma_1 * I   - kappa * R
  dRvdt <- gamma_1 * Iv  - kappa * Rv
  dRddt <- gamma_1 * Id  - kappa * Rd
  ## combine results into a single vector
  dxdt <- c(dSdt,dSvdt,dSddt,dVdt,dIdt,dIvdt,dIddt,dRdt,dRvdt,dRddt)
  ## return result as a list!
  list(dxdt)
}

ext.closed.sir.model.v2 <- function (t, x, params) {
  ## first extract the state variables
  N   <- sum(x[1:12])
  Sn   <- x['S']/N
  Sv1  <- x['Sv1']/N
  Sv2  <- x['Sv2']/N
  Sd  <- x['Sd']/N
  V   <- x['V']/N
  In   <- x['I']/N
  Iv1  <- x['Iv1']/N
  Iv2  <- x['Iv2']/N
  Id  <- x['Id']/N
  Rn   <- x['R']/N
  Rv  <- x['Rv']/N
  Rd  <- x['Rd']/N
  
  beta11  <- params["beta11"];  beta12   <- params["beta12"]; beta13 <- params["beta13"]; 
  beta21  <- params["beta21"];  beta22   <- params["beta22"]; beta23 <- params["beta23"]; 
  beta31  <- params["beta31"];  beta32   <- params["beta32"]; beta33 <- params["beta33"]; 
  gamma_1 <- params["gamma_1"]; gamma_2  <- params["gamma_2"];
  nu_1    <- params["nu_1"];    nu_1star <- params["nu_1star"]; nu_2 <- params["nu_2"]
  kappa   <- params["kappa"]; proportional_mixing <- params["proportional_mixing"]
  alpha   <- params["alpha"]; a <- params["a"]
  

  ## now code the model equations
  Sd <- a - Id - Rd
  I  <- Id + In
  Iv <- Iv1 + Iv2
  Sv <- Sv1 + Sv2
  S  <- Sd + Sn
  prop_mix_V <- 1
  prop_mix <- 1
  if(proportional_mixing){
    prop_mix_V <- beta22 / (beta11 * S + beta22 * (1 -S))
    prop_mix <- beta11 / (beta11 * S + beta22 * (1 -S))
  }
  dSddt  <-  - beta11 * prop_mix * Sd  *  I - 
    beta12 * prop_mix *  Sd  * Iv + 
    kappa * Rd
  dSndt  <-  - beta11 * prop_mix * Sn  *  I - 
    beta12 * prop_mix * Sn  * Iv - 
    nu_1 * Sn + 
    kappa * Rn
  dSv1dt <-  nu_1star * (1 - alpha) * Sv2 + 
    nu_1 * (1 - alpha) * Sn -
    nu_2 * Sv1 - 
    beta21 * prop_mix_V * Sv1 *  I - 
    beta22 * prop_mix_V * Sv1 * Iv 
  dSv2dt <-  - nu_1star * Sv2 +
    nu_2 * V + 
    nu_2 * Sv1 -
    beta21 * prop_mix_V * Sv2 *  I - beta22 * prop_mix_V * Sv2 * Iv + 
    kappa * Rv
  
  dVdt  <-  nu_1 * alpha * Sn +
    nu_1star * alpha * Sv2 -
    nu_2 * V +
    nu_1star * Rv +
    nu_1 * Rn
    
  dIddt  <- (beta11 * I + beta12 * Iv) * prop_mix * Sd  - gamma_1 * Id
  dIndt  <- (beta11 * I + beta12 * Iv) * prop_mix * Sn  - gamma_1 * In
  dIv1dt <- (beta21 * I + beta22 * Iv) * prop_mix_V * Sv1 - gamma_1 * Iv1
  dIv2dt <- (beta21 * I + beta22 * Iv) * prop_mix_V * Sv2 - gamma_1 * Iv2
  
  dRndt <- gamma_1 * In  - kappa * Rn - nu_1 * Rn
  dRvdt <- gamma_1 * Iv  - kappa * Rv - nu_1star * Rv
  dRddt <- gamma_1 * Id  - kappa * Rd
  ## combine results into a single vector
  dxdt <- c(dSndt,dSv1dt,dSv2dt,dSddt,dVdt,dIndt,dIv1dt,dIv2dt,dIddt,dRndt,dRvdt,dRddt)
  ## return result as a list!
  list(dxdt)
}

get.coords.for.label <- function(lookup, xmin, xmax, ymin, ymax, 
                     width = 50, height = 0.1 * (ymax - ymin)){
  if( max(xmin + width, lookup - width) ){
    
  }
  list(xmin = max(xmin + width/2, lookup - width),
       xmax = min(xmax - width/2, lookup + width),
       ymax = ymax - height, ymin = ymax - 2*height)
}


# Run the application 
shinyApp(ui = ui, server = server)

