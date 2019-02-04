library("shiny")

# Define UI for application that plots random distributions
shinyUI(fluidPage(

  # Application title
  titlePanel("Simulation for Mediation Analysis with Outliers"),

  # First row for input and output
  fluidRow(

    # Side panel with a slider input for number of observations
    column(4,
      wellPanel(
        strong("Methods to compare:"),
        checkboxInput("standard", "Standard bootstrap", value = TRUE),
        checkboxInput("huberized", "Huberized bootstrap", value = FALSE),
        checkboxInput("median", "Bootstrap based on median regression",
                      value = FALSE),
        checkboxInput("robust", "ROBMED", value = TRUE),
        sliderInput("n_obs", "Number of observations:",
                    min = 50, max = 1000, value = 250, step = 1),
        sliderInput("a", "Path a:", min = -1, max = 1, value = 0.2, step = 0.001),
        sliderInput("b", "Path b:", min = -1, max = 1, value = 0.2, step = 0.001),
        sliderInput("c", "Path c:", min = -1, max = 1, value = 0.2, step = 0.001),
        sliderInput("n_out", "Number of outliers:",
                    min = 0, max = 10, value = 5, step = 1),
        sliderInput("d", "Expected distance of outliers:",
                    min = 0, max = 15, value = 9.37, step = 0.01),
        numericInput("seed", "Seed of random number generator:",
                     value = as.numeric(gsub("-", "", as.character(Sys.Date())))
        )
      )
    ),

    # Show a plot of the generated distribution
    column(8,
      plotOutput("scatterPlotMatrix"),
      plotOutput("densityPlot")
    )

  ),

  # Second row for some description and references
  fluidRow(
    column(12,
      h3("Description"),
      p("You can compare various bootstrap methods for mediation analysis on
         simulated data: the standard bootstrap test of Preacher & Hayes (2004,
         2008), the Huberized bootstrap test of Zu & Yuan (2010), the bootstrap
         test based on median regression of Yuan & MacKinnon (2014), and the
         fast and robust boostrap test (ROBMED) of Alfons, Ates & Groenen
         (2018)."),

      p("The default settings follow the simulation design of Zu & Yuan (2010).
         The good data points are generated via the mediation model",
         withMathJax("$$ \\begin{align}
                         M &= a X + e_{M}, \\\\
                         Y &= b M + c X + e_{Y},
                         \\end{align} $$"),
        "with standard normal error terms.  You can adjust the total number of
         observations, the values of the coefficients, the number of outliers,
         and the expected distance of the outliers from the main point cloud."),

      p("As this simulation is just for illustration, the bootstrap procedures
         use only 1000 replicates.  For each selected method, the bootstrap
         distribution of the indirect effect is shown together with a shaded
         area representing the 95% confidence interval."),
      h3("References"),
      HTML('<p>
              Alfons, A., Ates, N.Y. and Groenen, P.J.F. (2018) A robust
              bootstrap test for mediation analysis.
              <em>ERIM Report Series in Management</em>, Erasmus Research
              Institute of Management.  URL
              <a href="https://hdl.handle.net/1765/109594">https://hdl.handle.net/1765/109594</a>.
            </p>'),
      HTML("<p>
              Preacher, K.J. and Hayes, A.F. (2004) SPSS and SAS procedures for
              estimating indirect effects in simple mediation models.
              <em>Behavior Research Methods, Instruments, & Computers</em>,
              <strong>36</strong>(4), 717-731.
            </p>"),
      HTML("<p>
              Preacher, K.J. and Hayes, A.F. (2008) Asymptotic and resampling
              strategies for assessing and comparing indirect effects in multiple
              mediator models.
              <em>Behavior Research Methods</em>,
              <strong>40</strong>(3), 879-891.
            </p>"),
      HTML("<p>
              Yuan, Y. and MacKinnon, D.P. (2014) Robust mediation analysis
              based on median regression.
              <em>Psychological Methods</em>,
              <strong>19</strong>(1), 1-20.
            </p>"),
      HTML("<p>
              Zu, J. and Yuan, K.-H. (2010) Local influence and robust
              procedures for mediation analysis.
              <em>Multivariate Behavioral Research</em>,
              <strong>45</strong>(1), 1-44.
            </p>")
    )
  )

))
