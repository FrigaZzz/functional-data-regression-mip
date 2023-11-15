
#' Plots the functional predictors for each observation
#'
#' @param result A list containing the functional predictors and the response variable
#' @param n The number of observations
#' @param M The number of functional predictors
#' @param t_range The range of the time variable
#'
#' @return A plot of the functional predictors for each observation
#'
#' @export
plot_functional_predictors <- function(result, n, M, t_range) {
  for (i in 1:n) {
    for (m in 1:M) {
      # Plot one plot per combination (i, m)
      par(mfrow = c(1, 1))  # Reset to default, plotting one plot per device
      
      # Plot the functional predictor m for observation i
      plot(t_range, result$X[i, m, ], type = "l", lty = 1, col = rainbow(M)[m], 
           main = paste("Observation", i, "- Functional Predictor", m, " Y=", result$Y[i]), 
           xlab = "t", ylab = "X(t)")
      
      # Add a legend to the plot. Since there's one curve per plot, the legend will have only one entry
      legend("topright", legend = paste("Predictor", m), col = rainbow(M)[m], lty = 1, cex = 0.5)
    }
  }
}

#' Plots the design matrix
#'
#' @param result The design matrix
#' @param title The title of the plot
#' @param xlab The label for the x-axis
#' @param ylab The label for the y-axis
#'
#' @return A plot of the design matrix
#'
#' @export
plot_design_matrix <- function(result,title,xlab,ylab) {
  library(fields)
  par(mfrow = c(1, 1))
  # Create an image plot with a color bar
  image.plot(1:ncol(result), 1:nrow(result), t(result), 
             main = paste("Matrix ",title), xlab = xlab, ylab = ylab)
  axis(2)
  axis(1, at = 1:ncol(result))
}