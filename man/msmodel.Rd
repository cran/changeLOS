\name{msmodel}
\alias{msmodel}
\title{ multi-state model }
\description{
  Makes a 'msmodel' - object to describe a 'multi-state model'
}
\usage{
msmodel(state.names, tra, cens.name)
}
\arguments{
  \item{state.names}{character vector of the state names}
  \item{tra}{ quadratic matrix of locigal values. TRUE/FALSE: transition is/is not possible}
  \item{cens.name}{character string, name of the censoring variable}
}
\value{
  An object of the class 'msmodel'. The object is a list of
  \item{msmodel$tra}{ quadratic matrix of locigal values. TRUE/FALSE: transition is/is not possible}
  \item{msmodel$states}{ numeric vector, the set of \emph{names} of the states. The internal
                         representation will be numbers 1, 2, 3, ...; if censoring occurrs,
                         the highest number will be the censoring code}
  \item{msmodel$state.names}{ character vector of the statenames}
  \item{msmodel$transitions}{matrix with two columns, 1.column: state 'from', 2.column: state 'to',the number of rows is the number of possibles transitions}
}

\author{ Matthias Wangler \email{mw@imbi.uni-freiburg.de}}
\examples{
my.model <- msmodel(c("0","1","2","3"),cens.name="cens")
}

\keyword{models}
