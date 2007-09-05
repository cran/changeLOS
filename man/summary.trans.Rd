\name{summary.trans}
\alias{summary.trans}
\title{Summary method for objects of class 'trans'}
\description{
  Summary method for objects of class \code{\link{trans}}.
}
\usage{
summary.trans(object, ...)
}
\arguments{
  \item{object}{An object of class \code{trans}.}
  \item{...}{other arguments}
}
\value{
  A list is returned by the function \code{summary.trans} with the
  following elements:
  \item{nrtransitions:}{ a matrix with
                  \describe{
                  \item{column 1}{the state from where a transition occurs}
                  \item{column 2}{the state to which a transition occurs}
                  \item{column 3}{the number of transitions}
               }    
  }
  \item{state.names:}{ vector with the names of the states}
  \item{times}{ the transition times}
  \item{nr.before}{ matrix with the number in each state just before the transition times}
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
