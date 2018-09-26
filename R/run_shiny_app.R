# --------------------------------------
# Author: Andreas Alfons
#         Erasmus Universiteit Rotterdam
# --------------------------------------

#' Shiny app: simulation for mediation analysis with outliers
#'
#' Compare various bootstrap methods for mediation analysis on simulated data.
#'
#' The default settings follow the simulation design of Zu & Yuan (2010).  You
#' can adjust the total number of observations, the values of the coefficients
#' in the mediation model, the number of outliers, as well as the expected
#' distance of the outliers from the main point cloud.
#'
#' As this simulation is just for illustration, the bootstrap procedures use
#' only 1000 replicates.  For each selected methods, the bootstrap distribution
#' of the indirect effect is shown together with a shaded area representing the
#' 95\% confidence interval.
#'
#' @author Andreas Alfons
#'
#' @seealso \code{\link{test_mediation}}
#'
#' @examples
#' \dontrun{
#' run_shiny_app()
#' }
#'
#' @keywords documentation
#'
#' @importFrom shiny runApp
#' @export

run_shiny_app <- function() {
  # find application folder
  folder <- system.file("shiny-examples", "simulation", package = "robmed")
  if (folder == "") {
    stop("Could not find shiny app.  Try re-installing package 'robmed'.")
  }
  # run shiny app
  shiny::runApp(folder, display.mode = "normal")
}
