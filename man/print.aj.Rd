\name{print.aj}
\alias{print.aj}
\alias{summary.aj}

\title{Print and summary method for objects of class 'aj'}
\description{
   Print and summary method for objects of class \code{aj}.
}
\usage{
print.aj(x, ...)

summary.aj(object, ...)
}

\arguments{
  \item{x}{ an object of class 'aj' }
  \item{object}{ an object of class 'aj' }
  \item{...}{other arguments}
}

\author{Matthias Wangler \email{mw@imbi.uni-freiburg.de}}

\seealso{ \code{\link{aj}} }
\examples{
data(los.data)
my.observ <- prepare.los.data(x=los.data)
my.model <- msmodel(c("0","1","2","3"),cens.name="cens")
my.trans <- trans(model=my.model,observ=my.observ)
my.aj <- aj(my.trans,s=0,t=80)
summary(my.aj)
print(my.aj)
}

\keyword{print}
