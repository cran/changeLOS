\name{print.msmodel}
\alias{print.msmodel}
\title{Print method for objects of class 'msmodel'}
\description{
  Print method for objects of class \code{msmodel}.
}
\usage{
\method{print}{msmodel}(x, ...)
}
\arguments{
  \item{x}{An object of class \code{msmodel}.}
  \item{...}{other arguments} 
}

\author{Matthias Wangler \email{mw@imbi.uni-freiburg.de}}

\examples{
my.model <- msmodel(c("0","1","2","3"),cens.name="cens")
print(my.model)
}

\keyword{print}
