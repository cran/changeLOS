\name{exampleclos}
\alias{exampleclos}
\title{ example for clos }
\description{
  An example program estimating change in length of stay associated with
  an intermediate event.
}
\usage{
exampleclos()
}

\value{
  An object of class 'clos'. The object is a list of
  \item{cLOS      }{change in LOS}
  \item{trans     }{an object of class 'trans'}
  \item{e.given.1 }{estimates E(\eqn{T|X_s = }intermediate event) for
    all observed transition times \eqn{s}}
  \item{e.given.0 }{estimates E(\eqn{T|X_s = }initial state)  for
    all observed transition times \eqn{s}}
  \item{weights   }{weights for the weighted average }
  \item{called }{how the function was called}
  \item{patients}{total number of observed patients}
  \item{patients.discharge}{number of patients being discharged}
  \item{patients.death}{number of patients who die}
  \item{patients.cens}{number of patients being censored, i. e. for whom
  neither discharge or death was observed}
  \item{cases}{number of patients who experienced the IE}
  \item{cases.discharge}{number of patients who experienced the IE being discharged}
  \item{cases.death}{number of patients who experienced the IE and died}
  \item{cases.cens}{number of patients who experienced the IE and were censored}
}
\note{required packages: survival }
\author{ Matthias Wangler \email{mw@imbi.uni-freiburg.de}}
\seealso{ \code{\link{clos}} }

\keyword{programming}
