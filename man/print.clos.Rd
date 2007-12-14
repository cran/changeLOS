\name{print.clos}
\alias{print.clos}
\title{Print method for objects of class 'clos'}
\description{
  Print method for objects of class \code{clos}.
}
\usage{
\method{print}{clos}(x, ...)
}
\arguments{
  \item{x}{An object of class \code{clos}.}
  \item{...}{other arguments}
}

\author{Matthias Wangler \email{mw@imbi.uni-freiburg.de}}
\seealso{  \code{\link{clos}}}
\examples{
data(los.data)
my.observ <- prepare.los.data(x=los.data)
trans <- matrix(FALSE,4,4)
diag(trans) <- TRUE
trans[1,] <- TRUE
trans[2,3:4] <- TRUE
my.model <- msmodel(c("0","1","2","3"),trans,cens.name="cens")
los <- clos(model=my.model,observ=my.observ)
print(los)
}

\keyword{print}
