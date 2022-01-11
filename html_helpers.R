# file contains only helper html headers/styles for the frontend off the app.

get_html_head <- function() {
  tags$head(
  tags$link(rel="stylesheet", 
            href="https://cdn.jsdelivr.net/npm/katex@0.10.1/dist/katex.min.css", 
            integrity="sha384-dbVIfZGuN1Yq7/1Ocstc1lUEm+AT+/rCkibIcC/OmWo5f0EA48Vf8CytHzGrSwbQ",
            crossorigin="anonymous"),
  tags$style(HTML("
                                      div.MathJax_Display{
                                      text-align: left !important;
                                      }
                                      ")),
  HTML('<script defer src="https://cdn.jsdelivr.net/npm/katex@0.10.1/dist/katex.min.js" integrity="sha384-2BKqo+exmr9su6dir+qCw08N2ZKRucY4PrGQPPWU1A7FtlCGjmEGFqXCv5nyM5Ij" crossorigin="anonymous"></script>'),
  HTML('<script defer src="https://cdn.jsdelivr.net/npm/katex@0.10.1/dist/contrib/auto-render.min.js" integrity="sha384-kWPLUVMOks5AQFrykwIup5lo0m3iMkkHrD0uJ4H5cjeGihAutqP0yW0J6dpFiVkI" crossorigin="anonymous"></script>'),
  HTML('
                           <script>
                           document.addEventListener("DOMContentLoaded", function(){
                           renderMathInElement(document.body, {
                           delimiters: [{left: "$", right: "$", display: false}]
                           });
                           })
                           </script>') 
  )
}

get_short_description <- function(){
    p(style="text-align: justify;", "The model is a part of the article:", 
      h4(style="text-align: justify;", "Unfavourable impact of COVID-19 vaccination passes on epidemic dynamics"),
      h5(style="text-align: justify;", "Tyll Krueger, Krzysztof Gogolewski, Marcin Bodych, 
          Anna Gambin, Giulia Giordano, Sarah Cuschieri , Thomas Czypionka,  Matjaz Perc, Elena Petelos,
          Magdalena Rosińska and Ewa Szczurek"),
      p(style="text-align: justify;", "Thanks to this app you can visualize possible COVID-19 epidemic dynamics 
        depending on crucial parameters, i.e. vaccine effectiveness, 
        vaccination rate, mixing type and the restrictions for VP holders 
        and for the rest of the population. For detailed model dynamics formulation and description of parameters see", 
        actionLink("link_to_tabpanel_a", "VAP-SIRS: formulation"), "tab."),
      h5("Quick tutorial"),
      p(style="text-align: justify;", "1. The Sidebar panel allows you to play with the initial conditions and parameters of the simulation. 
      On each change of any of the parameters the plots in the main panel update.
      The parameters setting you can see on the start is called in the article ", em("the reference setting"), "."),
      p(style="text-align: justify;", "2. There are", strong("five different settings"), "that correspond to the Figures 2. and 3. from the article. 
        These are the signatures of restrictions setting that introduce a specific epidemics dynamis."),
      p(style="text-align: justify;", "3. For some settings the values of Infected per 1M are low.
        For clear picture you can use the", strong("log-scale Y-axis"), ".",
        "Additionally, to see the specific values of Infected group, you can display a table values for selected time points."),
      br(),
      "The code of the shiny app and the model is available on-line:", 
        a("VAP-SIRS @ Github.", 
          href = "https://github.com/storaged/VAP-SIRS")
    )
}
