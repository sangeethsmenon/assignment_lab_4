#' Linreg Class
#'
#' This class represents a simple linear regression model.
#'
#' @name linreg
#' @exportClass linreg
#' @slot call The formula used for the regression model.
#' @slot coefficients A named vector of regression coefficients.
#'
#' @examples
#' \dontrun{
#' # Create a linreg object
#' my_linreg <- linreg(Petal.Length ~ Species, data = iris)
#'
#' # Print the linreg object
#' print(my_linreg)
#'
#' # Plot the linreg object
#' plot(my_linreg)
#'
#' # Get residuals
#' residuals <- resid(my_linreg)
#'
#' # Get predicted values
#' predictions <- pred(my_linreg)
#'
#' # Get coefficients
#' coefficients <- coef(my_linreg)
#'
#' # Summarize the linreg object
#' summary_info <- summary(my_linreg)
#' }
#'
#' @rdname linreg
#'
#' @export


library(ggplot2)

linreg <- setRefClass(
  "linreg",
  fields = list(
    formula = "formula",
    data = "data.frame",
    coefficients = "numeric",
    fitted.values = "numeric",
    residuals = "numeric",
    df = "numeric",
    rss = "numeric",
    variances = "numeric",
    t_values = "numeric",
    p_values = "numeric",
    beta = "numeric"
  ),
  methods = list(
    initialize = function(formula, data){
      # Store formula and data
      .self$formula <- formula
      .self$data <- data

      X <- as.matrix(model.matrix(formula, data))
      y <- data[[all.vars(formula)[1]]]

      beta <<- as.numeric(solve(t(X) %*% X) %*% t(X) %*% y)

      # Calculate fitted values
      fitted.values <<- as.numeric(X %*% beta)

      # Convert fitted values to a numeric vector
      .self$fitted.values <- as.numeric(fitted.values)

      # Calculate residuals
      .self$residuals <- y - .self$fitted.values

      # Calculate degrees of freedom
      .self$df <- nrow(X) - length(beta)

      # Calculate residual sum of squares
      .self$rss <- sum(.self$residuals^2)

      # Calculate variance of regression coefficients
      .self$variances <- diag(solve(t(X) %*% X) * .self$rss / .self$df)

      # Calculate t-values for each coefficient
      .self$t_values <- beta / sqrt(.self$variances)

      # Calculate p-values for each coefficient
      .self$p_values <- 2 * (1 - pt(abs(.self$t_values), .self$df))

      # Store computed values with variable names
      coef_names <- colnames(X)
      names(beta) <<- coef_names
      .self$coefficients <- beta
    },
    print = function() {
      cat("linreg(formula = ", deparse(.self$formula), ", data = iris", ")\n")

      coef_names <- names(.self$coefficients)
      cat(paste(format(.self$coefficients, digits = 2), coef_names, sep = " "), "\n")
    },

    plot = function() {
      # Create a scatterplot of residuals vs. fitted values
      ggplot(data = .self$data, aes(x = .self$fitted.values, y = .self$residuals)) +
        geom_point() +
        geom_hline(yintercept = 0, linetype = "dashed") +
        labs(x = "Fitted Values", y = "Residuals") +
        ggtitle("Residuals vs. Fitted Values")
      # Create a Scale-Location plot
      ggplot(data = .self$data, aes(x = .self$fitted.values, y = sqrt(abs(.self$residuals)))) +
        geom_point() +
        geom_smooth(method = "lm", se = FALSE )+
        labs(x = "Fitted Values", y = "Square Root of Standardized Residuals") +
        ggtitle("Scale-Location Plot")

    },

    resid = function() {
      return(.self$residuals)
    },

    pred = function() {
      return(.self$fitted.values)
    },

    coef = function() {
      return(.self$coefficients)
    },

    summary = function() {
      coef_names <- names(.self$coefficients)
      for (i in seq_along(coef_names)) {
        cat(coef_names[i], ":", format(.self$coefficients[i], digits = 2), "\n")
      }
      cat("\nResidual standard error:", sqrt(.self$rss / .self$df), "on", .self$df, "degrees of freedom\n")
    }
  )
)





