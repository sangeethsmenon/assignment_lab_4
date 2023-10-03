\name{mod_object}
\alias{mod_object}
\docType{data}
\title{
 Model Object - Description and Usage
}
\description{
The mod_object is a fundamental component of the lrmodels package, representing a linear regression model. This documentation provides comprehensive insights into the mod_object, including its attributes, structure, and practical utilization within the package. Gain a deeper understanding of how to create, manipulate, and employ linear regression models using this object.
}
\usage{data("mod_object")}
\format{
  The format is:
List of 13
 $ coefficients : Named num [1:3] 1.46 2.8 4.09
  ..- attr(*, "names")= chr [1:3] "(Intercept)" "Speciesversicolor" "Speciesvirginica"
 $ residuals    : Named num [1:150] -0.062 -0.062 -0.162 0.038 -0.062 ...
  ..- attr(*, "names")= chr [1:150] "1" "2" "3" "4" ...
 $ effects      : Named num [1:150] -46.026 4.347 20.45 0.058 -0.042 ...
  ..- attr(*, "names")= chr [1:150] "(Intercept)" "Speciesversicolor" "Speciesvirginica" "" ...
 $ rank         : int 3
 $ fitted.values: Named num [1:150] 1.46 1.46 1.46 1.46 1.46 ...
  ..- attr(*, "names")= chr [1:150] "1" "2" "3" "4" ...
 $ assign       : int [1:3] 0 1 1
 $ qr           :List of 5
  ..$ qr   : num [1:150, 1:3] -12.2474 0.0816 0.0816 0.0816 0.0816 ...
  .. ..- attr(*, "dimnames")=List of 2
  .. .. ..$ : chr [1:150] "1" "2" "3" "4" ...
  .. .. ..$ : chr [1:3] "(Intercept)" "Speciesversicolor" "Speciesvirginica"
  .. ..- attr(*, "assign")= int [1:3] 0 1 1
  .. ..- attr(*, "contrasts")=List of 1
  .. .. ..$ Species: chr "contr.treatment"
  ..$ qraux: num [1:3] 1.08 1.05 1.09
  ..$ pivot: int [1:3] 1 2 3
  ..$ tol  : num 1e-07
  ..$ rank : int 3
  ..- attr(*, "class")= chr "qr"
 $ df.residual  : int 147
 $ contrasts    :List of 1
  ..$ Species: chr "contr.treatment"
 $ xlevels      :List of 1
  ..$ Species: chr [1:3] "setosa" "versicolor" "virginica"
 $ call         : language lm(formula = Petal.Length ~ Species, data = iris)
 $ terms        :Classes 'terms', 'formula'  language Petal.Length ~ Species
  .. ..- attr(*, "variables")= language list(Petal.Length, Species)
  .. ..- attr(*, "factors")= int [1:2, 1] 0 1
  .. .. ..- attr(*, "dimnames")=List of 2
  .. .. .. ..$ : chr [1:2] "Petal.Length" "Species"
  .. .. .. ..$ : chr "Species"
  .. ..- attr(*, "term.labels")= chr "Species"
  .. ..- attr(*, "order")= int 1
  .. ..- attr(*, "intercept")= int 1
  .. ..- attr(*, "response")= int 1
  .. ..- attr(*, ".Environment")=<environment: R_GlobalEnv>
  .. ..- attr(*, "predvars")= language list(Petal.Length, Species)
  .. ..- attr(*, "dataClasses")= Named chr [1:2] "numeric" "factor"
  .. .. ..- attr(*, "names")= chr [1:2] "Petal.Length" "Species"
 $ model        :'data.frame':	150 obs. of  2 variables:
  ..$ Petal.Length: num [1:150] 1.4 1.4 1.3 1.5 1.4 1.7 1.4 1.5 1.4 1.5 ...
  ..$ Species     : Factor w/ 3 levels "setosa","versicolor",..: 1 1 1 1 1 1 1 1 1 1 ...
  ..- attr(*, "terms")=Classes 'terms', 'formula'  language Petal.Length ~ Species
  .. .. ..- attr(*, "variables")= language list(Petal.Length, Species)
  .. .. ..- attr(*, "factors")= int [1:2, 1] 0 1
  .. .. .. ..- attr(*, "dimnames")=List of 2
  .. .. .. .. ..$ : chr [1:2] "Petal.Length" "Species"
  .. .. .. .. ..$ : chr "Species"
  .. .. ..- attr(*, "term.labels")= chr "Species"
  .. .. ..- attr(*, "order")= int 1
  .. .. ..- attr(*, "intercept")= int 1
  .. .. ..- attr(*, "response")= int 1
  .. .. ..- attr(*, ".Environment")=<environment: R_GlobalEnv>
  .. .. ..- attr(*, "predvars")= language list(Petal.Length, Species)
  .. .. ..- attr(*, "dataClasses")= Named chr [1:2] "numeric" "factor"
  .. .. .. ..- attr(*, "names")= chr [1:2] "Petal.Length" "Species"
 - attr(*, "class")= chr "lm"
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
data(mod_object)
## maybe str(mod_object) ; plot(mod_object) ...
}
\keyword{datasets}
