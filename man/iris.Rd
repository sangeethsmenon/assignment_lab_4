\name{iris}
\alias{iris}
\docType{data}
\title{
Iris Dataset - Description and Usage
}
\description{
This documentation provides a comprehensive description of the Iris dataset, including its origin, structure, and intended usage.
}
\usage{data("iris")}
\format{
  A data frame with 150 observations on the following 5 variables.
  \describe{
    \item{\code{Sepal.Length}}{a numeric vector}
    \item{\code{Sepal.Width}}{a numeric vector}
    \item{\code{Petal.Length}}{a numeric vector}
    \item{\code{Petal.Width}}{a numeric vector}
    \item{\code{Species}}{a factor with levels \code{setosa} \code{versicolor} \code{virginica}}
  }
}
\details{
  The iris dataset is a famous dataset in the field of statistics and machine learning.
  It contains measurements of sepal and petal length and width for three species of iris flowers: setosa, versicolor, and virginica.

  The dataset is often used as a benchmark for classification and clustering algorithms.
  Researchers and data scientists use it to demonstrate various data analysis and visualization techniques.

}

\examples{
data(iris)
## maybe str(iris) ; plot(iris) ...
}
\keyword{datasets}
