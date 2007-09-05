\name{print.progdismodel}
\alias{print.progdismodel}
\title{Print method for objects of class 'progdismodel'}
\description{
  Print method for objects of class \code{progdismodel}.
}
\usage{
print.progdismodel(x, ...)
}
\arguments{
  \item{x}{An object of class \code{progdismodel}.}
  \item{...}{other arguments} 
}

\author{Matthias Wangler \email{mw@imbi.uni-freiburg.de}}

\examples{
data(los.data)

p <- prepare.progdismodel(data=los.data)

pdm <- progdismodel(p$model, p$observ)

print(pdm)
}

\keyword{print}
