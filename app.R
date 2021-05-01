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
library(shinyjs)
source("html_helpers.R")
source("initial_parameters.R")
source("model_description_functions.R")

plan(multisession)

# Define UI for application that draws a histogram

ui <- dashboardBody(
  tabBox(width = 12,
         tabPanel("VIP-SIRS",
                  #fluidPage(
                    fluidPage(
                      get_html_head() , 
                      shinyjs::useShinyjs(),
                      fluidRow(
                        column(6, 
                               titlePanel("VAP-SIRS: Vaccination Passes model"),
                               get_short_description()),
                        column(6, align="center", img(src = "model.png", height = 455, width = 490))
                      ),
                      fluidRow(
                      
                        #column(4,
                        #       mainPanel(h4("Vaccination:")),
                        #       checkboxInput("fast_vacc", "fast (100% in 200 days)", value = F),
                        #       checkboxInput("variant", "Fast loss of immunity (~200 days)", value = F)),
                        #column(4,
                        #       mainPanel(h4("Only two betas:")),
                        #       checkboxInput("two_betas", withMathJax("Use only $\\beta_{S_V,I_V}^{I_V}$ others are equal") , value = T),
                        #       checkboxInput("two_betas", withMathJax("Use only $\\beta$ and $\\beta_V$") , value = T),
              
                        #       checkboxInput("fix_gamma", withMathJax("Fix: $\\gamma = \\frac{1}{6}$") , value = T),
                        #       checkboxInput("fix_reduction", withMathJax("Fix: $r = 0.8$") , value = F),
                        #       checkboxInput("fix_nu1star", withMathJax("Fix: $\\upsilon = \\upsilon_r$") , value = T)) 
                       
                               
                        #column(4,
                               #mainPanel(h4("Contacts:")),
                               #checkboxInput("check1", "Use means", value = F),
                               #checkboxInput("go_crazy", withMathJax("Go crazy: $\\beta_{S_V,I_V}^{I_V} = 4 \\cdot \\beta_{S,I}^{I}$"), value = F))
                        ),
                    
                      #fluidRow(
                      #  mainPanel(h4("Selected parameters' values: "))
                      #  ),
                      #uiOutput("current_parameters"),
                      fluidRow( 
                        #column(4,
                        #       checkboxInput("log_scale", "Use log-scale on Y axis", value = F), 
                        #       checkboxInput("displayTable", "Display table values", value = F),
                        #       sliderInput("lookup", withMathJax("$t_i$"),  min = 0,  max = 2*365, value = 400, step=1 )),
                        
                               #checkboxInput("proportional_mixing", "Use proportional mixing", value = F)),
                        #column(4,
                        #       selectInput("plot_type", "Select plot type:",
                        #                   c("Daily incidence" = "daily"#,
                        #                     #"Infections per 1M" = "per1M"
                        #                     )))
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
                      
                      sidebarLayout(
                        sidebarPanel(
                          fluidRow(
                            column(6,
                                   h5("General properties"))
                          ),
                          fluidRow(
                            column(6,
                                   checkboxInput("two_betas", withMathJax("Use VAP-SIRS setting") , value = T)),
                            column(6,
                                   checkboxInput("fix_nu1star", withMathJax("Fix: $\\upsilon = \\upsilon_r$") , value = T))
                          ),
                          fluidRow(
                            column(6,
                                   h5("Initial values"))
                          ),
                          fluidRow(
                            column(6,
                                   sliderInput("S", withMathJax("$S = (1-d) S_N + d  S_D$"), min = 0,  max = 100000, value = S_INIT, step=1000)),
                            column(2,
                                   sliderInput("Sv1", withMathJax("$S_1$"), min = 0,  max = 100000, value = 0, step=1000)),
                            column(6,
                                   sliderInput("I", withMathJax("$I = (1-d) I_N + d  I_D$"),   min = 0, max = 20000,  value = I_INIT, step=1))
                          ), fluidRow(
                            column(6,
                                   sliderInput("R", withMathJax("$R = (1-d ) R_N + d  R_D$"),  min = 0,  max = 100000, value = R_INIT, step=100 )),
                            column(6,
                                   sliderInput("V", withMathJax("$V$"),  min = 0,  max = 100000, value = V_INIT , step=100 ))
                          ), fluidRow(
                            column(12,
                                   h5("Infection parameters setup"))
                          ), fluidRow(
                            column(4, sliderInput("infectivity",withMathJax("$\\beta_0$") , 
                                                  min = 0,max = 1,value = BETA_0, step = 0.001)),
                            column(4,  sliderInput("freedom",withMathJax("$f$") , 
                                                   min = 0,max = 1,value = F_GENERAL, step = 0.01)),
                            column(4, sliderInput("beta11",withMathJax("$\\beta$"),
                                                  min = 0,max = 1,value = BETA, step = 0.001))
                          ), 
                          fluidRow(
                            column(4, sliderInput("gamma", withMathJax("$\\gamma$"),  
                                                  min = 0,max = 1,value = GAMMA, step = 0.001)),
                            column(4, sliderInput("freedom_V",withMathJax("$f_V$") , 
                                                  min = 0,max = 1,value = F_V, step = 0.01)),
                            column(4, sliderInput("beta22",withMathJax("$\\beta_V$"), 
                                                  min = 0,max = 1,value = BETA, step = 0.001))
                          ), fluidRow(
                            column(12,
                                   h5("Vaccination parameters setup"))
                          ), fluidRow(
                            column(4, sliderInput("nu1",withMathJax("$\\upsilon$"), 
                                                  min = 0,max = 0.05, value = 1/UPSILON, step = 0.0005)),
                            column(4, sliderInput("nu2",withMathJax("$\\omega$"), 
                                                  min = 0,max = 0.05,value = 1/OMEGA, step = 0.0005)),
                            column(4, sliderInput("kappa", withMathJax("$\\kappa$"), 
                                                  min = 0,max = 0.05,value = 1/KAPPA, step = 0.0005))
                          ), 
                          fluidRow(
                            column(4, sliderInput("nu1star",withMathJax("$\\upsilon_r$"), 
                                                  min = 0, max = 0.05, value = 1/UPSILON_R, step = 0.0005)),
                            column(4, sliderInput("sceptics", withMathJax("$d$"),  
                                                  min = 0,max = 1,value = D_FRACTION, step = 0.01)),
                            column(4, sliderInput("alpha1", withMathJax("$a$"), 
                                                  min = 0,max = 1,value = A_1, step = 0.01))
                            #sliderInput("gamma", withMathJax("$\\gamma$"),  min = 0,max = 1,value = GAMMA, step = 0.001),
                            
                          )
                          ),
                        mainPanel(
                          fluidRow(
                            column(6,
                                   h5("Plotting Settings"))
                          ),
                          fluidRow(
                            column(3,
                                   checkboxInput("log_scale", "Use log-scale on Y axis", value = F)),
                            column(3,
                                   checkboxInput("displayTable", "Display table values", value = F)),
                            column(3,
                                   sliderInput("lookup", withMathJax("$t_i$"),  
                                               min = 1,  max = 2*365, value = 400, step=1 ))
                          ),
                          fluidRow(
                            column(6,
                                   h4("Proportional Mixing")),
                            column(6,
                                   h4("Preferential Mixing"))
                          ),
                          fluidRow(
                            plotOutput("extDistPlotv2", height = "600px")
                          ),
                          fluidRow(
                            column(3, 
                                   downloadButton("downloadPlot", "Download plot")),
                            column(3, 
                                   downloadButton("downloadData", "Download .csv")),
                            column(3,
                                   selectInput("plot_type", "Select plot type:",
                                               c("Daily incidence" = "daily"#,
                                                 #"Infections per 1M" = "per1M"
                                               )))
                          )
                        ) # mainPanel
                      ), #SidebarLayout
                        
                      # fluidRow(
                      #   column(6,
                      #        h4("Initial values"))
                      # ),
                      # fluidRow(
                      #   column(2,
                      #          sliderInput("S", withMathJax("$S = (1-d) \\cdot S_N + d \\cdot S_D$"), min = 0,  max = 100000, value = S_INIT, step=1000)),
                      #   column(2,
                      #          sliderInput("Sv1", withMathJax("$S_1$"), min = 0,  max = 100000, value = 0, step=1000)),
                      #   column(2,
                      #          sliderInput("I", withMathJax("$I = (1-d) \\cdot I_N + d \\cdot I_D$"),   min = 0, max = 20000,  value = I_INIT, step=1)),
                      #   column(2,
                      #          sliderInput("R", withMathJax("$R = (1-d ) \\cdot R_N + d \\cdot R_D$"),  min = 0,  max = 100000, value = R_INIT, step=100 )),
                      #   column(2,
                      #          sliderInput("V", withMathJax("$V$"),  min = 0,  max = 100000, value = V_INIT , step=100 ))
                      # ),
                      fluidRow(
                      #  column(6,
                      #         h4("Parameter setup"))
                      ),       
                      
                      fluidRow(
                        column(2,
                        #       sliderInput("infectivity",withMathJax("$\\beta_0$") , min = 0,max = 1,value = BETA_0, step = 0.001),
                        #       #sliderInput("beta11",label = withMathJax("$\\beta_{S_N,I_N}^{I_N}$"), min = 0,max = 1,value = BETA, step = 0.001)
                        #       sliderInput("beta11",withMathJax("$\\beta$"), min = 0,max = 1,value = BETA, step = 0.001),
                               sliderInput("beta12",withMathJax("$\\beta_{S_N,I_N}^{I_V}$") , min = 0,max = 1,value = BETA, step = 0.001),
                               sliderInput("beta13",withMathJax("$\\beta_{S_N,I_N}^{I_D}$") , min = 0,max = 1,value = BETA, step = 0.001)
                               
                               ),
                        column(2,
                               #sliderInput("freedom",withMathJax("$f$") , min = 0,max = 1,value = F_GENERAL, step = 0.01),
                               sliderInput("beta21",withMathJax("$\\beta_{S_V,I_V}^{I_N}$"), min = 0,max = 1,value = BETA, step = 0.001),
                               #sliderInput("beta22",withMathJax("$\\beta_V$"), min = 0,max = 1,value = BETA, step = 0.001),
                               #sliderInput("beta22",withMathJax("$\\beta_{S_V,I_V}^{I_V}$"), min = 0,max = 1,value = BETA_V, step = 0.001),
                               sliderInput("beta23",withMathJax("$\\beta_{S_V,I_V}^{I_D}$"), min = 0,max = 1,value = BETA, step = 0.001)
                               
                               ),
                        column(2,
                               #sliderInput("freedom_V",withMathJax("$f_V$") , min = 0,max = 1,value = F_V, step = 0.01),
                              
                               sliderInput("beta31",withMathJax("$\\beta_{S_D,I_D}^{I_N}$"), min = 0,max = 1,value = BETA, step = 0.001),
                               sliderInput("beta32",withMathJax("$\\beta_{S_D,I_D}^{I_V}$"), min = 0,max = 1,value = BETA, step = 0.001),
                               sliderInput("beta33",withMathJax("$\\beta_{S_D,I_D}^{I_D}$"), min = 0,max = 1,value = BETA, step = 0.001)
                               
                               ),
                        column(2,
                               #sliderInput("reduction",withMathJax("$r$") , min = 0,max = 1,value = 1, step = 0.001), 
                               #sliderInput("nu1",withMathJax("$\\upsilon$"), min = 0,max = 0.05, value = 1/UPSILON, step = 0.0005),
                               #sliderInput("nu1star",withMathJax("$\\upsilon_r$"), min = 0, max = 0.05, value = 1/UPSILON_R, step = 0.0005),
                               #sliderInput("nu2",withMathJax("$\\omega$"), min = 0,max = 0.05,value = 1/OMEGA, step = 0.0005),
                               #sliderInput("kappa", withMathJax("$\\kappa$"),  min = 0,max = 0.05,value = 1/KAPPA, step = 0.0005)
                               
                               ),
                        column(2,
                               #sliderInput("sceptics", withMathJax("$d$"),  min = 0,max = 1,value = D_FRACTION, step = 0.01),
                               #sliderInput("gamma", withMathJax("$\\gamma$"),  min = 0,max = 1,value = GAMMA, step = 0.001),
                               #sliderInput("alpha1", withMathJax("$a$"),  min = 0,max = 1,value = A_1, step = 0.01),
                               sliderInput("alpha2", withMathJax("$a_2$"),  min = 0,max = 1,value = A_2, step = 0.01)
                               )
                        #column(2,
                        #       sliderInput("S", withMathJax("$S = (1-d) \\cdot S_N + d \\cdot S_D$"), min = 0,  max = 100000, value = S_INIT, step=1000),
                        #       sliderInput("Sv1", withMathJax("$S_1$"), min = 0,  max = 100000, value = 0, step=1000),
                        #       sliderInput("I", withMathJax("$I = (1-d) \\cdot I_N + d \\cdot I_D$"),   min = 0, max = 20000,  value = I_INIT, step=1),
                        #       sliderInput("R", withMathJax("$R = (1-d ) \\cdot R_N + d \\cdot R_D$"),  min = 0,  max = 100000, value = R_INIT, step=100 ),
                        #       sliderInput("V", withMathJax("$V$"),  min = 0,  max = 100000, value = V_INIT , step=100 ),
                               
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
                    uiOutput("ode_model_description"),
                    
                    fluidRow(
                      mainPanel(h4("Parameters description: "))
                    ),
                    uiOutput("description")
                    )
         )
         )
    )


update_on_freedom <- function(input, output, session){
  updateSliderInput(session, "beta11", value =  (1 - input$freedom) * input$infectivity ) 
  updateSliderInput(session, "beta22", value =   (1 - input$freedom_V) * input$infectivity ) #* input$reduction 
  update_betas(beta_val = (1 - input$freedom) * input$infectivity, session)
}

update_betas <- function(beta_val, session){
  updateSliderInput(session, "beta12", value =  beta_val ) 
  updateSliderInput(session, "beta13", value =  beta_val ) 
  updateSliderInput(session, "beta21", value =  beta_val ) 
  updateSliderInput(session, "beta23", value =  beta_val ) 
  updateSliderInput(session, "beta31", value =  beta_val ) 
  updateSliderInput(session, "beta32", value =  beta_val ) 
  updateSliderInput(session, "beta33", value =  beta_val ) 
}

signatures_df <- data.frame(f =  c(0.8, 0.63, 0.3, 0.8, 0.5), 
                            fv =  c(0.05, 0.05, 0.05, 0.4, 0.4))
signatures_f <- lapply(1:5, function(i) 
  function(input, output, session){ 
    updateSliderInput(session, "freedom", value = signatures_df[i,"f"]); 
    updateSliderInput(session, "freedom_V", value = signatures_df[i,"fv"]);
    update_on_freedom(input, output, session);
})

names(signatures_f) <- paste("signature", as.roman(1:5), sep = "")


# Define server logic required to draw a histogram
server <- function(input, output, session) {
  
  #output$beta11_slider <- renderUI({
    #if(!input$two_betas){
    #  sliderInput("beta11",label = withMathJax("$\\beta_{S_N,I_N}^{I_N}$"), min = 0,max = 1,value = BETA, step = 0.001)
    #} else {
    #  sliderInput("beta11",label = list(withMathJax(), "$\\beta$"), min = 0,max = 1,value = BETA, step = 0.001)
    #}
  #})
  
  observeEvent(input$signatureI,{ if (input$signatureI) signatures_f$signatureI(input, output, session) })
  observeEvent(input$signatureII,{ if (input$signatureII) signatures_f$signatureII(input, output, session) })
  observeEvent(input$signatureIII,{ if (input$signatureIII) signatures_f$signatureIII(input, output, session) })
  observeEvent(input$signatureIV,{ if (input$signatureIV) signatures_f$signatureIV(input, output, session) })
  observeEvent(input$signatureV,{ if (input$signatureV) signatures_f$signatureV(input, output, session) })
  
  observeEvent(input$two_betas,{
    #if(input$two_betas){
    #  updateSliderInput(session, "beta11", label = withMathJax("$\\beta$"))
    #  updateSliderInput(session, "beta22", label = paste(withMathJax("$\\beta_V$")))
    #} else {
    #  updateSliderInput(session, "beta11", label = withMathJax("$\\beta_{S_N,I_N}^{I_V}$"))
    #  updateSliderInput(session, "beta22", label = withMathJax("$\\beta_{S_V,I_V}^{I_V}$"))    
    #}
    shinyjs::toggle(id="beta12", ! input$two_betas) 
    shinyjs::toggle(id="beta13", ! input$two_betas) 
    shinyjs::toggle(id="beta21", ! input$two_betas) 
    shinyjs::toggle(id="beta23", ! input$two_betas) 
    shinyjs::toggle(id="beta31", ! input$two_betas) 
    shinyjs::toggle(id="beta32", ! input$two_betas) 
    shinyjs::toggle(id="beta33", ! input$two_betas) 
    shinyjs::toggle(id="alpha2", ! input$two_betas) 
    shinyjs::toggle(id="Sv1", ! input$two_betas) 
    toggleState(id="beta12") 
    toggleState(id="beta13")
    toggleState(id="beta21")
    toggleState(id="beta23")
    toggleState(id="beta31")
    toggleState(id="beta32")
    toggleState(id="beta33")
    toggleState(id="alpha2")
    toggleState(id="Sv1")
  })

  observeEvent(input$fix_nu1star, {toggleState(id="nu1star", !input$fix_nu1star)})
  observeEvent(!input$two_betas, {toggleState(id="two_betas", !input$two_betas)})
  
  
  observeEvent(input$nu1,{
    #toggleState(id="nu1star")
    if (input$fix_nu1star) { 
      updateSliderInput(session, "nu1star", value = input$nu1 ) 
    } 
  })
  
  observeEvent(input$freedom | input$freedom_V | input$infectivity , { #} | input$reduction ,{
    #print("OBSERVED CHANGE")
    if (input$two_betas) { 
      update_on_freedom(input, output, session)
    } 
  })
  observeEvent(input$displayTable, { 
    toggle(id="lookup", input$displayTable)
  })
  

  output$ode_model_description <- render.model.ode.description(input, output, session)
  output$description <- render.model.params.description(input, output, session)
  output$current_parameters <- render.model.current.parameters(input, output, session)
  

  observeEvent(input$plot_type ,{
    #print("PLOT_TYPE")
    if (input$plot_type == "daily") { 
      output$extDistPlotv2 <- showDailyIncidents()
    } else if (input$plot_type == "per1M") {
      output$extDistPlotv2 <- showPer1MCases()
    }
  })
  
  showDailyIncidents <- function(){
    renderPlot({
      #plots <- plotInput()
      ## remember that here the naming is messed up proportional<->preferential
      ggarrange(plotInput_daily_proportional(), plotInput_daily_preferential(),
                plotInput_muller_proportional(), plotInput_muller_preferential(), 
                align = "hv", nrow = 2, ncol = 2) 
    })
  }
  
  showPer1MCases <- function(){ 
    renderPlot({
      plots <- plotInput()
      
      ggarrange(plots[["g.cum"]], plots[["g.cum"]],
                plotInput_muller_proportional(), plotInput_muller_preferential(), nrow = 2, ncol = 2) 
    })
  }
  
  plotInput_daily_proportional <-   reactive({
    create_daily_plot(proportional_mixing = F)
  })
  
  plotInput_daily_preferential <-   reactive({
    create_daily_plot(proportional_mixing = T)
  })
  
  plotInput_muller_proportional <-   reactive({
    create_muller_plot(proportional_mixing = F)
  })
  
  plotInput_muller_preferential <-   reactive({
    create_muller_plot(proportional_mixing = T)
  })
  
  
  plotInputDownload <- reactive({
    #plots <- plotInput()
    ggarrange(create_daily_plot(proportional_mixing = F), 
              create_muller_plot(proportional_mixing = F), 
              align = "hv", nrow = 2, ncol = 1) 
  })
  
  
  create_daily_plot <- function(proportional_mixing = F){
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
    
    parms1 <- c(beta11=input$beta11,beta12=input$beta12,beta13=input$beta13,
                beta21=input$beta21,beta22=input$beta22,beta23=input$beta23,
                beta31=input$beta31,beta32=input$beta32,beta33=input$beta33,
                gamma_1=input$gamma, gamma_2=input$gamma,
                nu_1=input$nu1,nu_2=input$nu2,nu_1star=input$nu1star,
                kappa=input$kappa,a=input$sceptics, proportional_mixing = F,
                alpha=input$alpha1)
    
    if(proportional_mixing){
      parms1["proportional_mixing"] <- T
    }
    
    ode(
      func=ext.closed.sir.model.v2,
      y=xstart,
      times=times,
      parms=parms1
    ) %>%
      as.data.frame() -> out
    
    tmp.dat <- cbind(out[, c("Iv1", "Iv2")], 
                     I = rowSums(out[, c("I", "Id")]), 
                     Isum = rowSums(out[, c("I", "Iv1", "Iv2", "Id")]))
    
    
    dat.diff <- rbind(
      cbind(time = out[-1, c("time")], 
            alpha = input$alpha1, 
            data.frame(diff(as.matrix(
              tmp.dat 
            ))) + (1/6) * tmp.dat[-1,] 
      )) %>%
      gather(variable,value,-c(time,alpha) )
    
    
    dat.diff$alpha <- factor(dat.diff$alpha)
    my.labels <- paste0(as.double(levels(dat.diff$alpha)) * 100, "%")
    levels(dat.diff$alpha) <- my.labels
    print(levels(dat.diff$variable))
    print(levels(factor(dat.diff$variable)))
    dat.diff$variable <- factor(dat.diff$variable, levels = c( "Isum", "I", "Iv1", "Iv2"))
    #dat.diff$variable <- (dat.diff$variable, levels = c("I", "Isum", "Iv1", "Iv2"))
    print(levels(factor(dat.diff$variable)))
    
    if( input$displayTable ) {
      my.label.coords <- get.coords.for.label(input$lookup, 
                                              xmin = min(dat.diff$time), xmax = max(dat.diff$time),
                                              ymin = 0, ymax = max(dat.diff$value) * 1e6)
      
      for.table <- dat.diff[dat.diff$time==input$lookup, ]
      print(for.table)
      mytable <- cbind(Group=c("All", "I", "I1", "I2"),
                       data.frame(Value=c(round(for.table[for.table$variable == "Isum",]$value * 1e6), 
                                      round(for.table[for.table$variable == "I",]$value * 1e6), 
                                      round(for.table[for.table$variable == "Iv1",]$value * 1e6), 
                                      round(for.table[for.table$variable == "Iv2",]$value * 1e6))))
      
      print(mytable)
    }
    
    # DAILY PROPORTIONAL 
    if(length(levels(dat.diff$alpha)) > 1){
      g.diff <- dat.diff  %>% 
        ggplot(aes(x=time, y=value * 1e6, 
                   color=variable, 
                   group = interaction(alpha, variable), 
                   linetype = alpha)) +
        scale_linetype_manual("Efficacy: ", values = c("solid", "dotted")) 
    } else {
      g.diff <- dat.diff  %>% 
        ggplot(aes(x=time, 
                   y=value * 1e6, 
                   color=variable, 
                   group = interaction(alpha, variable))) 
    }
    
    g.diff <- g.diff +
      scale_color_manual("Group:",
                         labels = expression("All",I,   I[1], I[2]),
                         values = c("#ff3747","#dee500",  "#ff458f", "#5e3c99",  "darkblue") ) + 
      geom_line(size = 1, alpha = 0.8, linetype = "dashed") +
      geom_line(size = 1, 
                data = subset(dat.diff, variable == "Isum"), 
                linetype = "solid", 
                color = "#ff3747") +
      theme_minimal() +
      
      guides(color = guide_legend(ncol = 1, override.aes = list(size = 2, linetype = "solid"))) + 
      theme(text = element_text(size=18), 
            legend.text = element_text(size=18),
            legend.key.width = unit(1,"cm"), 
            legend.position = "right", 
            legend.box = "vertical") + 
      labs(x = 'Time (days)', 
           y = "Daily incidence per 1M pop."
      )
    if(input$log_scale){
      g.diff <- g.diff + 
        scale_y_log10( breaks = trans_breaks("log10", function(x) 10^x),
                       labels = trans_format("log10", math_format(10^.x)), 
                       limits = c(10^-10, 10^5)) + 
        annotation_logticks(sides="l") +
        labs( y = "Daily incidence per 1M pop. [log]")
    } 
    if(input$displayTable) {
      g.diff <- g.diff + geom_vline(xintercept = input$lookup) + 
        annotation_custom(tableGrob(mytable, rows = NULL, 
                                    theme = ttheme_default(
                                      core=list(bg_params = list(fill = "lightgrey", col=NA, alpha = 0.8),
                                                fg_params=list(fontface=3)))), 
                          xmin = my.label.coords$xmin, 
                          xmax = my.label.coords$xmax,
                          ymax = my.label.coords$ymax, 
                          ymin = my.label.coords$ymin) 
    }
    
    g.diff
  }
  

  create_muller_plot <- function(proportional_mixing = F){
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
    
    parms1 <- c(beta11=input$beta11,beta12=input$beta12,beta13=input$beta13,
                beta21=input$beta21,beta22=input$beta22,beta23=input$beta23,
                beta31=input$beta31,beta32=input$beta32,beta33=input$beta33,
                gamma_1=input$gamma, gamma_2=input$gamma,
                nu_1=input$nu1,nu_2=input$nu2,nu_1star=input$nu1star,
                kappa=input$kappa,a=input$sceptics, proportional_mixing = F,
                alpha=input$alpha1)
    
    if(proportional_mixing){
      parms1["proportional_mixing"] <- T
    }
    
    ode(
      func=ext.closed.sir.model.v2,
      y=xstart,
      times=times,
      parms=parms1
    ) %>%
      as.data.frame() -> out
    
    tmp.dat <- cbind(out[, c("Iv1", "Iv2")], 
                     I = rowSums(out[, c("I", "Id")]), 
                     Isum = rowSums(out[, c("I", "Iv1", "Iv2", "Id")]))
    
    
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
    
    
    g.muller <- Muller_plot(melted.out, add_legend = TRUE) + 
      theme_minimal() + 
      guides(fill = guide_legend(ncol = 1)) + 
      scale_y_continuous(labels = scales::percent_format(accuracy = 1)) + 
      scale_fill_manual("Group:", values = c("#FFF97E", "#EFBB40", "#F29837", "#81D452", "#EE6D56", "#79BFF9"), 
                        labels = expression(S, S[1], S[2], V, I[Sigma], R[Sigma])) +
      theme(text = element_text(size=18),
            axis.title.x =  element_blank(),
            legend.text=element_text(size=18),
            legend.position="right", 
            legend.box="vertical") + 
      labs(x = "Time (days)", 
           y = "Population fraction")  
     
     g.muller
  }
  

  
  
  plotInput <- reactive({
    
    times <- seq(from=0,to=2*365,by=1)
    Npop <- sum(input$S + input$I + input$R + input$Sv1 + input$V)
    xstart <- c(S  = (1 - input$sceptics) * input$S/Npop, 
                Sv1 = input$Sv1/Npop, 
                Sv2 = 0,
                Sd = (input$sceptics) * input$S/Npop,
                V  = input$V/Npop,
                I  = (1 - input$sceptics) * input$I/Npop, 
                Iv1 = 0, Iv2 = 0, 
                Id = (input$sceptics) * input$I/Npop,
                R  = (1 - input$sceptics) * input$R/Npop, 
                Rv = 0,
                Rd = (input$sceptics) * input$R/Npop )
    
 


    parms1 <- c(beta11=input$beta11,beta12=input$beta12,beta13=input$beta13,
                beta21=input$beta21,beta22=input$beta22,beta23=input$beta23,
                beta31=input$beta31,beta32=input$beta32,beta33=input$beta33,
                gamma_1=input$gamma, gamma_2=input$gamma,
                nu_1=input$nu1,nu_2=input$nu2,nu_1star=input$nu1star,
                kappa=input$kappa,a=input$sceptics, proportional_mixing = F,
                alpha=input$alpha1)
    
    ode(
      func=ext.closed.sir.model.v2,
      y=xstart,
      times=times,
      parms=parms1
    ) %>%
      as.data.frame() -> out
    
    parms2<- parms1
    parms2["alpha"] = input$alpha2
    
    ode(
      func=ext.closed.sir.model.v2,
      y=xstart,
      times=times,
      parms=parms2
    ) %>%
      as.data.frame() -> out2
    
    
    parms3<- parms1
    parms3["proportional_mixing"] = T
    
    ode(
      func=ext.closed.sir.model.v2,
      y=xstart,
      times=times,
      parms=parms3
    ) %>%
      as.data.frame() -> out3
      
    dat3 <- 
      cbind(out3[, c("time", "Iv1", "Iv2")],
            I = rowSums(out3[, c("I", "Id")]), 
            Isum = rowSums(out3[, c("I", "Iv1", "Iv2", "Id")]), 
            alpha = parms3["alpha"]) %>%
      gather(variable,value,-c(time,alpha) ) 
    
    dat3$variable <- factor(dat3$variable, c("I", "Iv1", "Iv2", "Isum"))
    print(levels(dat3$variable))
    levels(dat3$variable) <- c("Not vaccinated (I)", "After vaccination (Iv1)", "After vaccination (Iv2)", "All")
      
    
    my.caption <- paste0("Scenario: ", 
                         "Go crazy: ", ifelse(input$go_crazy, "T", "F"), 
                         "; Vaccination: ", ifelse(input$fast_vacc, "fast", "slow"), 
                         "; Immunity: ", ifelse(input$variant, "200 days", "400 days"))
    
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
    
    
    muller.out.prop.mix   <- out3
    muller.out.prop.mix$S <- muller.out.prop.mix$S + muller.out.prop.mix$Sd
    muller.out.prop.mix$I <- muller.out.prop.mix$I + muller.out.prop.mix$Id + muller.out.prop.mix$Iv1 + muller.out.prop.mix$Iv2
    muller.out.prop.mix$R <- muller.out.prop.mix$R + muller.out.prop.mix$Rd + muller.out.prop.mix$Rv
    muller.out.prop.mix   <- muller.out.prop.mix[, c("time", "S", "Sv1", "Sv2", "V", "I", "R")]
    
    muller.out.prop.mix %>% gather(variable,value,-c(time) ) -> melted.out.prop.mix
    colnames(melted.out.prop.mix) <- c("Generation", "Group_id", "Frequency")
    melted.out.prop.mix$Unique_id <- paste0(melted.out.prop.mix$Group_id, "_", melted.out.prop.mix$Generation)
    melted.out.prop.mix$Identity <- factor(melted.out.prop.mix$Group_id, levels = c("S", "Sv1", "Sv2", "V", "I", "R"))
    melted.out.prop.mix$Population <- melted.out.prop.mix$Frequency * 1e6
    
    dat <- rbind(
      cbind(out[, c("time", "Iv1", "Iv2")], I = rowSums(out[, c("I", "Id")]), Isum = rowSums(out[, c("I", "Iv1", "Iv2", "Id")]), alpha = input$alpha1), 
      cbind(out2[, c("time", "Iv1", "Iv2")],I = rowSums(out2[, c("I", "Id")]), Isum = rowSums(out2[, c("I", "Iv1", "Iv2", "Id")]), alpha = input$alpha2)) %>%
      gather(variable,value,-c(time,alpha) ) 
    
    tmp.dat <- cbind(out[, c("Iv1", "Iv2")], 
                     I = rowSums(out[, c("I", "Id")]), 
                     Isum = rowSums(out[, c("I", "Iv1", "Iv2", "Id")]))
    
    tmp.dat2 <- cbind(out2[, c("Iv1", "Iv2")], 
                     I = rowSums(out2[, c("I", "Id")]), 
                     Isum = rowSums(out2[, c("I", "Iv1", "Iv2", "Id")]))
    
    tmp.dat3 <- cbind(out3[, c("Iv1", "Iv2")], 
                      I = rowSums(out3[, c("I", "Id")]), 
                      Isum = rowSums(out3[, c("I", "Iv1", "Iv2", "Id")]))
    
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

    
    dat.diff$alpha <- factor(dat.diff$alpha)
    my.labels <- paste0(as.double(levels(dat.diff$alpha)) * 100, "%")
    levels(dat.diff$alpha) <- my.labels
    dat.diff$variable <- factor(dat.diff$variable, levels = c("I", "Iv1", "Iv2", "Isum"))
    
    
    dat.diff.prop.mix <- rbind(
      cbind(time = out3[-1, c("time")], 
            alpha = input$alpha1, 
            data.frame(diff(as.matrix(
              tmp.dat3 
            ))) + (1/6) * tmp.dat3[-1,] 
      ),
      cbind(time = out3[-1, c("time")], 
            alpha = input$alpha2, 
            data.frame(diff(as.matrix(
              tmp.dat3 
            ))) + (1/6) * tmp.dat3[-1,] 
      )) %>%
      gather(variable,value,-c(time,alpha) )
    
    
    dat.diff.prop.mix$alpha <- factor(dat.diff.prop.mix$alpha)
    my.labels <- paste0(as.double(levels(dat.diff.prop.mix$alpha)) * 100, "%")
    levels(dat.diff.prop.mix$alpha) <- my.labels
    dat.diff.prop.mix$variable <- factor(dat.diff.prop.mix$variable, levels = c("I", "Iv1", "Iv2", "Isum"))
    
    
    dat$alpha <- factor(dat$alpha)
    my.labels <- paste0(as.double(levels(dat$alpha)) * 100, "%")
    levels(dat$alpha) <- my.labels
    dat$variable <- factor(dat$variable, levels = c("I", "Iv1", "Iv2", "Isum"))
    levels(dat$variable) <- c("Not vaccinated (I)", "After vaccination (Iv1)", "After vaccination (Iv2)", "All")
    
    dat3$alpha <- factor(dat3$alpha)
    my.labels <- paste0(as.double(levels(dat3$alpha)) * 100, "%")
    levels(dat3$alpha) <- my.labels
    dat3$variable <- factor(dat3$variable, levels = c("I", "Iv1", "Iv2", "Isum"))
    levels(dat3$variable) <- c("Not vaccinated (I)", "After vaccination (Iv1)", "After vaccination (Iv2)", "All")
    
    
    In.max <- 1e6*out[out$time==input$lookup,]$I ; In.max.2 <- 1e6*out2[out2$time==input$lookup,]$I 
    In.max <- ifelse (In.max < 1, round(In.max, 3), round(In.max))
    In.max.2 <- ifelse (In.max.2 < 1, round(In.max.2, 3), round(In.max.2))
    
    Iv1.max <- 1e6*out[out$time==input$lookup,]$Iv1; Iv1.max.2 <- 1e6*out2[out2$time==input$lookup,]$Iv1
    Iv1.max <- ifelse (Iv1.max < 1, round(Iv1.max, 3), round(Iv1.max))
    Iv1.max.2 <- ifelse (Iv1.max.2 < 1, round(Iv1.max.2, 3), round(Iv1.max.2))
    Iv2.max <- 1e6*out[out$time==input$lookup,]$Iv2; Iv2.max.2 <- 1e6*out2[out2$time==input$lookup,]$Iv2
    Iv2.max <- ifelse (Iv2.max < 1, round(Iv2.max, 3), round(Iv2.max))
    Iv2.max.2 <- ifelse (Iv2.max.2 < 1, round(Iv2.max.2, 3), round(Iv2.max.2))
    
    Isum.max   <- 1e6*rowSums(out[, c("I", "Iv1", "Iv2")])[input$lookup]; 
    Isum.max.2 <- 1e6*rowSums(out2[, c("I", "Iv1", "Iv2")])[input$lookup] 
    Isum.max   <- ifelse (Isum.max < 1, round(Isum.max, 3), round(Isum.max))
    Isum.max.2 <- ifelse (Isum.max.2 < 1, round(Isum.max.2, 3), round(Isum.max.2))
    

    mytable <- cbind(sites=c("Isum","I","Iv1","Iv2"),
                     data.frame(a=c(Isum.max, In.max, Iv1.max, Iv2.max), 
                                b=c(Isum.max.2, In.max.2, Iv1.max.2, Iv2.max.2)))

    names(mytable) <- c("group", my.labels)
    rownames(mytable) <- NULL
    
    cat(mytable)
    
    my.label.coords <- get.coords.for.label(input$lookup, 
                                  xmin = min(dat$time), xmax = max(dat$time),
                                  ymin = 0, ymax = max(dat$value) * 1e6)
    
    if(length(levels(dat$alpha)) > 1){
      g <- dat %>% 
        ggplot(aes(x=time,y=value * 1e6, color=variable, 
                   group = interaction(alpha, variable), linetype = alpha)) +
        scale_linetype_manual("Efficacy: ", values = c("solid", "dotted")) 
    } else {
      g <- dat %>% 
        ggplot(aes(x=time,y=value * 1e6, color=variable, 
                   group = interaction(alpha, variable))) 
    }
    
    g <- g + 
      scale_color_manual("Group: ",
                         labels = expression(I, I[1], I[2], "All"),
                         values = c("#e05f38", "#e08138", "#e0ae38", "#EE6D56", "darkblue") ) + 
      geom_line(size = 1) + 
      #### Annotation table with values
      geom_vline(xintercept = input$lookup) + 
      annotation_custom(tableGrob(mytable, rows = NULL, 
                                  theme = ttheme_default(
                                    core=list(bg_params = list(fill = "lightgrey", col=NA, alpha = 0.8),
                                              fg_params=list(fontface=3)))), 
                        xmin = my.label.coords$xmin, 
                        xmax = my.label.coords$xmax,
                        ymax = my.label.coords$ymax, 
                        ymin = my.label.coords$ymin) +
      theme_minimal() +
      guides(color = guide_legend(ncol = 1, override.aes = list(size = 2))) + 
      theme(text = element_text(size=18),
            legend.text = element_text(size=18),
            legend.key.width = unit(1,"cm"), 
            legend.position = "right", 
            legend.box = "vertical") + 
      labs(x = 'Time (days)', 
           y = "Infected in 1M pop.",
           caption = "Proportional mixing") 
    
    
    # DAILY PROPORTIONAL 
    if(length(levels(dat.diff$alpha)) > 1){
      g.diff <- dat.diff  %>% 
        ggplot(aes(x=time, y=value * 1e6, 
                   color=variable, 
                   group = interaction(alpha, variable), 
                   linetype = alpha)) +
        scale_linetype_manual("Efficacy: ", values = c("solid", "dotted")) 
    } else {
      g.diff <- dat.diff  %>% 
        ggplot(aes(x=time,y=value * 1e6, 
                   color=variable, 
                   group = interaction(alpha, variable))) 
    }

    g.diff <- g.diff +
      scale_color_manual("Group:",
                         labels = expression(I, I[1], I[2], "All"),
                         values = c("#dee500", "#ff458f", "#5e3c99", "#ff3747", "darkblue") ) + 
      geom_line(size = 1, alpha = 0.8, linetype = "dashed") +
      geom_line(size = 1, data = subset(dat.diff, variable == "Isum"), linetype = "solid", color = "#ff3747") +
      theme_minimal() +
      guides(color = guide_legend(ncol = 1, override.aes = list(size = 2))) + 
      theme(text = element_text(size=18), 
            legend.text = element_text(size=18),
            legend.key.width = unit(1,"cm"), 
            legend.position = "right", 
            legend.box = "vertical") + 
      labs(x = 'Time (days)', 
           y = "Daily incidence per 1M pop."
      )
    if(input$log_scale){
      g.diff <- g.diff + 
        scale_y_log10( breaks = trans_breaks("log10", function(x) 10^x),
                       labels = trans_format("log10", math_format(10^.x)), 
                       limits = c(10^-10, 10^5)) + 
        annotation_logticks(sides="l") +
        labs( y = "Daily incidence per 1M pop. [log]")
    } 
    # DAILY PREFERENTIAL 
    if(length(levels(dat.diff.prop.mix$alpha)) > 1){
      g.diff.prop.mix <- dat.diff.prop.mix  %>% 
        ggplot(aes(x=time, y=value * 1e6, 
                   color=variable, 
                   group = interaction(alpha, variable), 
                   linetype = alpha)) +
        scale_linetype_manual("Efficacy: ", values = c("solid", "dotted")) 
    } else {
      g.diff.prop.mix <- dat.diff.prop.mix  %>% 
        ggplot(aes(x=time, y=value * 1e6, 
                   color=variable, 
                   group = interaction(alpha, variable))) 
    }
    g.diff.prop.mix <- g.diff.prop.mix +
      scale_color_manual("Group:",
                         labels = expression(I, I[1], I[2], "All"),
                         values = c("#dee500", "#ff458f", "#5e3c99", "#ff3747", "darkblue") ) + 
      geom_line(size = 1, alpha = 0.8, linetype = "dashed") +
      geom_line(size = 1, data = subset(dat.diff.prop.mix, variable == "Isum"), linetype = "solid", color = "#ff3747") +
      theme_minimal() +
      guides(color = guide_legend(ncol = 1, override.aes = list(size = 2))) + 
      theme(text = element_text(size=18), 
            legend.text = element_text(size=18),
            legend.key.width = unit(1,"cm"), 
            legend.position = "right", 
            legend.box = "vertical") + 
      labs(x = 'Time (days)', 
           y = "Daily incidence per 1M pop."
      )
    if(input$log_scale){
      g.diff.prop.mix <- g.diff.prop.mix + 
        scale_y_log10( breaks = trans_breaks("log10", function(x) 10^x),
                       labels = trans_format("log10", math_format(10^.x)), 
                       limits = c(10^-10, 10^5)) + 
        annotation_logticks(sides="l") +
        labs( y = "Daily incidence per 1M pop. [log]")
    } 
   
    
    g.prop.mix <- dat3 %>% 
      ggplot(aes(x=time,y=value * 1e6, color=variable, group = variable)) +
      scale_color_manual("Group:",
                         labels = expression(I, I[1], I[2], "All"),
                         values = c("#e05f38", "#e08138", "#e0ae38", "#EE6D56", "darkblue") ) + 
      geom_line(size = 1) + 
      theme_minimal() +
      guides(color = guide_legend(ncol = 1, override.aes = list(size = 2))) + 
      theme(text = element_text(size=18), 
            legend.text=element_text(size=18),
            legend.key.width = unit(1,"cm"), 
            legend.position="right", 
            legend.box="vertical") + 
      labs(x='Time (days)', 
           y = "Infected in 1M pop.",
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

    
   g.muller <- Muller_plot(melted.out, add_legend = TRUE) + 
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
                     y = "Population fraction")  
   
   g.prop.mix.muller <- Muller_plot(melted.out.prop.mix, add_legend = TRUE) + 
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
          y = "Population fraction")  
   
   list(g.cum = g, g.muller = g.muller, g.daily = g.diff, g.daily.prop.mix = g.diff.prop.mix,
        g.prop.mix = g.prop.mix,  g.prop.mix.muller = g.prop.mix.muller)
   
                
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
  R <- Rn
  
  prop_mix_V <- 1
  prop_mix <- 1
  
  if(proportional_mixing){
    prop_mix_V <- beta22 / (beta11 * (S+I+R) + beta22 * (1 - (S+I+R)))
    prop_mix <- beta11 / (beta11 * (S+I+R) + beta22 * (1 -(S+I+R)))
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

