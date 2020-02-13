library("robmed")
library("shiny")

# Define server logic required to generate and plot a random distribution
shinyServer(function(input, output) {

  # Expression that generates a plot of the distribution. The expression
  # is wrapped in a call to renderPlot to indicate that:
  #
  #  1) It is "reactive" and therefore should be automatically
  #     re-executed when inputs change
  #  2) Its output type is a plot

  ## colors for plot
  colors <- c("#F8766D", "#F564E3", "#619CFF", "#00BFC4")

  ## function to generate data
  generate_data <- function(n_obs, a, b, c, n_out, d, seed) {
    # set seed of random number generator
    set.seed(seed)
    # generate clean data
    x <- rnorm(n_obs)
    m <- a * x + rnorm(n_obs)
    y <- b * m + c * x + rnorm(n_obs)
    # compute covariance matrix
    s_M <- a^2 + 1
    s_Y <- b^2 * (a^2 + 1) + c^2 + 2 * a * b * c + 1
    s_MX <- a
    s_YX <- a * b + c
    s_YM <- b * (a^2 + 1) + a * c
    Sigma <- matrix(c(1, s_YX, s_MX, s_YX, s_Y, s_YM, s_MX, s_YM, s_M),
                    nrow = 3, ncol = 3)
    # contaminate the data
    if(n_out > 0) {
      # mean of clean data
      mu <-  c(0, 0, 0)
      # define means of contaminated data
      mu_out <- c(0, 1, -1)
      mu_out <- mu_out / sqrt(mahalanobis(mu_out, center = mu, cov = Sigma))
      # contaminate the first observations
      i <- seq_len(n_out)
      # x[i] <- x[i] + d * mu_out[1]  # not necessary (mu_out[1] is set to 0)
      y[i] <- y[i] + d * mu_out[2]
      m[i] <- m[i] + d * mu_out[3]
    }
    # return data frame
    data.frame(X = x, M = m, Y = y)
  }


  ## generate data
  df <- reactive({
    generate_data(input$n_obs, input$a, input$b, input$c,
                  input$n_out, input$d, input$seed)
  })

  ## scatter plot matrix of generated data
  output$scatterPlotMatrix <- renderPlot({
    col_points <- rep.int(colors, c(input$n_out, 0, 0, input$n_obs-input$n_out))
    plot(df(), pch = 16, cex = 2, col = col_points, cex.axis = 1.5,
         cex.labels = 2, las = 1, oma = c(2.5, 2.5, 2.5, 10.5))
    legend("right", inset = c(-0.03, 0), legend = c("Good point", "Outlier"),
           pch = 16, pt.cex = 1.4, col = colors[c(4, 1)], cex = 1.1,
           y.intersp = 0.75, bty = "n", xpd = TRUE)
  })


  ## density plot of bootstrap distributions
  output$densityPlot <- renderPlot({

    ## perform mediation analysis with nonrobust bootstrap test
    # With only 1000 bootstrap replicates, sometimes there is a warning that
    # extreme order statistics are used as endpoints of confidence intervals.
    # Since this shiny app is just illustrative, such warnings are ignored.
    suppressWarnings({
      if (input$standard) {
        standard_boot <- test_mediation(df(), x = "X", y = "Y", m = "M",
                                        test = "boot", R = 1000,
                                        method = "regression",
                                        robust = FALSE)
      }
      if (input$huberized) {
        huberized_boot <- test_mediation(df(), x = "X", y = "Y", m = "M",
                                         test = "boot", R = 1000,
                                         method = "covariance",
                                         robust = TRUE)
      }
      if (input$median) {
        median_boot <- test_mediation(df(), x = "X", y = "Y", m = "M",
                                      test = "boot", R = 1000,
                                      method = "regression",
                                      robust = "median")
      }
      if (input$robust) {
        robust_boot <- test_mediation(df(), x = "X", y = "Y", m = "M",
                                      test = "boot", R = 1000,
                                      method = "regression",
                                      robust = "MM")
      }
    })

    # create object containing selected tests and select colors accordingly
    tests <- list()
    if (input$standard) tests$Standard <- standard_boot
    if (input$huberized) tests$Huberized <- huberized_boot
    if (input$median) tests$Median <- median_boot
    if (input$robust) tests$ROBMED <- robust_boot

    # select colors
    select <- c(input$standard, input$huberized, input$median, input$robust)
    selected_colors <- colors[select]

    # plot the density of the bootstrap distribution
    if (length(tests) > 0) {
      dp <- density_plot(tests) +
        geom_vline(xintercept = input$a * input$b) +
        scale_color_manual(values = selected_colors) +
        scale_fill_manual(values = selected_colors) +
        theme(title = element_text(size = 15),
              axis.text = element_text(size = 13),
              axis.title = element_text(size = 15),
              legend.title = element_text(size = 15),
              legend.text = element_text(size = 13))
      print(dp)
    }
  })
})
