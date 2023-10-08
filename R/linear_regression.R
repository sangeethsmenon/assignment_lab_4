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
    dname = "character",
    coefficients = "numeric",
    fitted.values = "numeric",
    residuals = "numeric",
    df = "numeric",
    rss = "numeric",
    variances = "numeric",
    t_values = "numeric",
    p_values = "numeric",
    std_error = "numeric"
  ),
  methods = list(
    initialize = function(formula, data) {
      # Store formula and data
      .self$formula <- formula
      .self$data <- data
      .self$dname <- "iris"  # Set the dataset name

      X <- as.matrix(model.matrix(formula, data))
      y <- data[[all.vars(formula)[1]]]

      beta <- as.numeric(solve(t(X) %*% X) %*% t(X) %*% y)

      # Calculate fitted values
      .self$fitted.values <- as.numeric(X %*% beta)

      # Calculate residuals
      .self$residuals <- y - .self$fitted.values

      # Calculate degrees of freedom
      .self$df <- nrow(X) - length(beta)

      # Calculate residual sum of squares
      .self$rss <- sum(.self$residuals^2)

      # Calculate variance of regression coefficients
      variances <<- diag(solve(t(X) %*% X) * .self$rss / .self$df)

      # Calculate t-values for each coefficient
      t_values <<- beta / sqrt(variances)

      # Calculate p-values for each coefficient
      p_values <<- 2 * (1 - pt(abs(t_values), .self$df))

      # Store computed values with variable names
      coef_names <- colnames(X)
      names(beta) <- coef_names
      .self$coefficients <- beta
      .self$std_error <- sqrt(variances)
      .self$t_values <- t_values
      .self$p_values <- p_values
    },

    print = function() {
      cat("Call:\n")
      cat(paste0("linreg(formula = ", deparse(.self$formula), ", data = ", .self$dname, ")\n"))
      cat("\nCoefficients:\n")
      #cat("linreg(formula = Petal.Length ~ Sepal.Width + Sepal.Length, data = iris)")
      # Print coefficients with proper formatting
      coef_names <- names(.self$coefficients)
      coef_values <- format(.self$coefficients, digits = 2)

      # Combine coefficient names and values into two separate lines
      coef_line <- paste(coef_names, collapse = " ")
      value_line <- paste(coef_values, collapse = " ")

      cat(coef_line, "\n")
      cat(value_line, "\n")
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
      cat("call:\n")
      cat(paste0("linreg(formula = ", deparse(.self$formula), ", data = ", .self$dname, ")\n"))
      cat("\nResiduals:\n")

      min_residual <- min(.self$residuals)
      q1_residual <- quantile(.self$residuals, 0.25)
      median_residual <- median(.self$residuals)
      q3_residual <- quantile(.self$residuals, 0.75)
      max_residual <- max(.self$residuals)
      resid_matrix <- matrix(0, nrow = 1, ncol = 5)
      colnames(resid_matrix) <- c("Min", "1Q", "Median", "3Q", "Max")
      resid_matrix[1, ] <- c(min_residual, q1_residual, median_residual, q3_residual, max_residual)
      base::print(resid_matrix, quote = FALSE)
      cat("\nCoefficients:\n")

      coef_matrix <- cbind(.self$coefficients, .self$std_error, .self$t_values, .self$p_values)

      colnames(coef_matrix) <- c("Estimate", "Std. Error", "t value", "Pr(>|t|)")

      coef_dataframe <- as.data.frame(coef_matrix)

      star <- ifelse(.self$p_values < 0.05, "***", "")
      coef_dataframe$`Pr(>|t|)` <- paste0(coef_dataframe$`Pr(>|t|)`, star)

      base::print(coef_dataframe)

      cat(paste0("\nResidual standard error: ", format(sqrt(.self$rss / .self$df), digits = 9),
                 " on ", .self$df, " degrees of freedom"))
    }
  )
)

