\name{print.trans}
\alias{print.trans}
\title{Print method for objects of class 'trans'}
\description{
  Print method for objects of class \code{\link{trans}}.
}
\usage{
print.trans(x, ...)
}
\arguments{
  \item{x}{An object of class \code{trans}.}
  \item{...}{other arguments}
}

\author{Matthias Wangler \email{mw@imbi.uni-freiburg.de}}

\examples{
data(los.data)
my.observ <- prepare.los.data(x=los.data)
my.model <- msmodel(c("0","1","2","3"),cens.name="cens")
my.trans <- trans(model=my.model,observ=my.observ)
summary(my.trans)
}

\keyword{print}
