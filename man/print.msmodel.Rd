\name{print.msmodel}
\alias{summary.msmodel}
\alias{print.msmodel}
\title{Print and summary method for objects of class 'msmodel'}
\description{
  Print and summary method for objects of class \code{msmodel}.
}
\usage{
print.msmodel(x, ...)

summary.msmodel(object, ...)

}
\arguments{
  \item{x}{An object of class \code{msmodel}.}
  \item{object}{An object of class \code{msmodel}.}
  \item{...}{other arguments} 
}

\author{Matthias Wangler \email{mw@imbi.uni-freiburg.de}}

\examples{
my.model <- msmodel(c("0","1","2","3"),cens.name="cens")
summary(my.model)
print(my.model)
}

\keyword{print}
