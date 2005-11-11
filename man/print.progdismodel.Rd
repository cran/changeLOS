\name{print.progdismodel}
\alias{summary.progdismodel}
\alias{print.progdismodel}
\title{Print and summary method for objects of class 'progdismodel'}
\description{
  Print and summary method for objects of class \code{progdismodel}.
}
\usage{
print.progdismodel(x, ...)

summary.progdismodel(object, ...)

}
\arguments{
  \item{x}{An object of class \code{progdismodel}.}
  \item{object}{An object of class \code{progdismodel}.}
  \item{...}{other arguments} 
}

\author{Matthias Wangler \email{mw@imbi.uni-freiburg.de}}

\examples{
data(los.data)

p <- prepare.progdismodel(data=los.data)

pdm <- progdismodel(p$model, p$observ)

print(pdm)

summary(pdm)
}

\keyword{print}
