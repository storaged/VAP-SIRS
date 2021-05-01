## LATEX BASED DESCRIPTION OF THE ODE EQUATIONS OF THE MODEL 
render.model.ode.description <- function(input, output, session) {
  renderUI({
  fluidRow(
    column(6,
           withMathJax(
             sprintf('$$\\frac{dS_D}{dt} = -  ( \\beta I + \\beta I_V ) S_D + \\kappa R_D $$ \n 
                     $$\\frac{dS_N}{dt}  = -  ( \\beta I + \\beta I_V ) S_N - \\upsilon  S_N + \\kappa R_N$$ \n 
                     $$\\frac{dS_1}{dt}  = -  ( \\beta I + \\beta_V I_V ) S_1 + \\upsilon_{r}  (1 - \\alpha) S_2 + \\upsilon (1 - a)  S_N - \\omega  S_1 $$ \n 
                     $$\\frac{dS_2}{dt}  = -  ( \\beta I + \\beta_V I_V ) S_2 + \\upsilon_{r}   S_2 + \\omega  S_1 + \\omega V + \\kappa R_V$$ \n 
                     $$\\frac{dI_D}{dt} =  ( \\beta I + \\beta I_V ) S_D  - \\gamma  I_D$$ \n 
                     $$\\frac{dI_N}{dt}  =  ( \\beta I + \\beta I_V ) S_N  - \\gamma  I_N$$ \n 
                     $$\\frac{dI_1}{dt}  =  ( \\beta I + \\beta_V I_V ) S_1  - \\gamma  I_1$$ \n
                     $$\\frac{dI_2}{dt}  =  ( \\beta I + \\beta_V I_V ) S_2  - \\gamma  I_2$$ \n
                     ')
             )),
    column(6,
           withMathJax(
             sprintf('$$\\frac{dV}{dt}   =  \\upsilon  a  S_N + \\upsilon_{r} a S_2 - \\omega V + \\upsilon_{r} R_V + \\upsilon R_N$$ \n
                     $$\\frac{dR_D}{dt}   =  \\gamma  I_D - \\kappa  R_D$$ \n
                     $$\\frac{dR_N}{dt}   =  \\gamma  I_N - \\kappa  R_N - \\upsilon R_N$$ \n
                     $$\\frac{dR_V}{dt}   =  \\gamma  I_V - \\kappa  R_V - \\upsilon_{r} R_V$$ \n
                     $$\\quad \\hbox{where also the following relations hold:}$$ \n
                     $$S_V   = S_1 + S_2$$ \n
                     $$S   = S_N + S_D$$ \n
                     $$I_V   = I_1 + I_2$$ \n
                     $$I   = I_N + I_D$$ \n
                     $$R   = R_N + R_D$$ \n
                     $$d = S_D + I_D + R_D $$
                     ')
             )
             )
             )
  })
}


## PARAMETERS DESCRIOPTON
render.model.params.description <- function(input, output, session) {
  renderUI({
  fluidRow(
    column(6, 
           withMathJax(
             sprintf('$$\\beta_0 \\hbox{ - transmissability of SARS-CoV-2} $$ \n
                     $$f, f_V \\hbox{ - contact restrictions put on passport holders (IP) and others, respectively} $$ \n
                     $$\\beta, \\beta_V \\hbox{ - contact (infectivity) rate, where } \\beta=\\beta_0(1-f) \\hbox{ and } \\beta_V=\\beta_0(1-f_V) $$\n
                     $$\\gamma \\hbox{ - recovery rate } $$ \n
                     $$\\kappa \\hbox{ - becoming susceptible after recovery rate } $$ ')
             )
             ),
    column(6, 
           withMathJax(
             sprintf('
                     $$\\upsilon \\hbox{ - vaccination rate }$$ \n
                     $$\\upsilon_r \\hbox{ - re-vaccination rate }$$ \n
                     $$a \\hbox{ - vaccination efficacy}$$ \n
                     $$\\omega \\hbox{ - loss of immunity rate }$$ \n
                     $$d \\hbox{ - fraction of population will not get vaccinated } $$'
                      )
             )
             )
             )
})
}

  
render.model.current.parameters <- function(input, output, session){
  renderUI({
  fluidRow(
    
    column(6, 
           if(input$two_betas) {
             withMathJax(
               sprintf('$$\\beta_0 = %.03f$$ \n 
                       $$f = %.03f, f_V = %.03f$$ \n 
                       $$\\beta = %.03f, \\beta_V = %.03f$$', 
                       input$infectivity,
                       input$freedom, input$freedom_V,
                       input$beta11, input$beta22
               )
             )
           } else {
             withMathJax(
               sprintf('$$\\beta_{S_N,I_N}^{I_N} = %.03f,  \\beta_{S_N,I_N}^{I_V} = %.03f, \\beta_{S_N,I_N}^{I_D} = %.03f$$ \n 
                       $$ \\beta_{S_V,I_V}^{I_N} = %.03f,  \\beta_{S_V,I_V}^{I_V} = %.03f, \\beta_{S_V,I_V}^{I_D} = %.03f $$ \n
                       $$ \\beta_{S_D,I_D}^{I_N} = %.03f,  \\beta_{S_D,I_D}^{I_V} = %.03f, \\beta_{S_D,I_D}^{I_D} = %.03f $$', 
                       input$beta11, input$beta12, input$beta13, 
                       input$beta21, input$beta22, input$beta23, 
                       input$beta31, input$beta32, input$beta33 
               )
             )
           }
    ),
    column(6, 
           withMathJax(
             sprintf('$$ \\upsilon = %.04f, \\omega = %.04f; \\quad $$ \n
                     $$ \\gamma = %.04f, \\kappa = %.04f; \\quad $$ \n
                     $$ a_1 =  %.04f,  a_2 =  %.04f $$', 
                     input$nu1, input$nu2, 
                     input$gamma, input$kappa, 
                     input$alpha1, input$alpha2 )
           )
    )
  )
})
}
