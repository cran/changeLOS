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
  An object of class c('clos'). The object is a list of
  \item{cLOS      }{change in LOS}
  \item{trans     }{an object of class 'trans'}
  \item{e.given.1 }{estimates E(LOS|\eqn{X_s = }intermediate event) for
    all observed transition times \eqn{s}, where \eqn{X_s} denotes the
    state by time \eqn{s}}
  \item{e.given.0 }{estimates E(LOS|\eqn{X_s = }initial state)  for
    all observed transition times \eqn{s}, where \eqn{X_s} denotes the
    state by time \eqn{s}}
  \item{phi2}{weighted average of the difference between phi2.case and phi2.control}
  \item{phi2.case}{estimates E(LOS 1(\eqn{X_LOS = }discharge)|\eqn{X_s = }intermediate event)}
  \item{phi2.control}{P(\eqn{X_LOS = }discharge|\eqn{X_s = }intermediate
    event) E(LOS|\eqn{X_s = }initial state)} 
  \item{phi3}{weighted average of the difference between phi3.case and phi3.control}
  \item{phi3.case}{estimates E(LOS 1(\eqn{X_LOS = }death)|\eqn{X_s = }intermediate event)}
  \item{phi3.control}{P(\eqn{X_LOS = }death|\eqn{X_s = }intermediate
    event) E(LOS|\eqn{X_s = }initial state)}  
  \item{empty.1}{event times with the group `intermediate, but no terminal event yet' being empty}
  \item{empty.0}{event times with the group `no intermediate or terminal event yet' being empty}
  \item{weights}{weights for the weighted average }
  \item{w.times}{time points corresponding to the weights}
  \item{called }{how the function was called}
  \item{patients}{total number of observed patients}
  \item{patients.discharge}{number of patients being discharged}
  \item{patients.death}{number of patients who die}
  \item{patients.cens}{number of patients being censored, i. e. for whom
  neither discharge or death was observed}
  \item{cases}{number of patients who experienced the IE}
  \item{cases.discharge}{number of patients who experienced the IE being discharged}
  \item{cases.death}{number of patients who experienced the IE and died}
  \item{cases.cens}{number of patients who experienced the IE and were
    censored}
}
\note{required packages: survival }
\author{ Matthias Wangler \email{mw@imbi.uni-freiburg.de}}
\seealso{ \code{\link{clos}} }

\keyword{programming}
