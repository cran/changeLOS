\name{aj}
\alias{aj}

\title{ Aalen-Johansen estimator}
\description{
  computes the Aalen-Johansen estimator for the matrix of transition
  probabilities \eqn{P(s,t)}. The entry \eqn{(l,m)} of the matrix
  denotes the estimated probability that state \eqn{m} has been
  reached by time \eqn{t} given state \eqn{l} has been occupied by
  time \eqn{s}. 
}
\usage{
aj(tr, s, t)
}

\arguments{
  \item{tr}{ an object of 'trans'}
  \item{s}{ begin of the time interval}
  \item{t}{ end of the time interval}
}
\details{
  The Aalen-Johansen estimator is considered in detail by Andersen et
  al. (1993). Usually, the process describing movements between states
  is considered to be Markovian; this assumption may be relaxed,
  cf. Datta and Satten (2001).

  The Aalen-Johansen estimator is a finite matrix product, one matrix
  for every observed transition time in the time interval \eqn{(s,t]}. These
  matrices can be obtained from \code{\link{trans}}.

}
\value{
  An object of class 'aj'. The object is a list of:
  \item{matrix}{ the matrix of transition probabilities \eqn{P(s,t)}}
  \item{start}{ the beginning \eqn{s} of the time interval \eqn{(s,t]}}
  \item{end}{ the end \eqn{t} of the time interval \eqn{(s,t]}}
  \item{times}{ the transition times in the interval \eqn{(s,t]}}
  \item{matrices}{ array of estimators for \eqn{P(s,u)} for all transition
    times \eqn{u} in \eqn{(s,t]}}
}
\references{
  P Andersen, O Borgan, R Gill, and N Keiding (1993). \emph{Statistical
    models based on counting processes.} New York: Springer

  S Datta, and G Satten (2001). Validity of the Aalen-Johansen estimators of stage occupation
  probabilities and Nelson-Aalen estimators of integrated transition
  hazards for non-Markov models. \emph{Statistics and Probability
    Letters} 55 (4), 403--411.

  Andersen and Keiding (2002). Multi-state models for event history
  analysis. \emph{Statistical Methods in Medical Research} 11 (2),
  91--115.
  }
\author{ Matthias Wangler \email{mw@imbi.uni-freiburg.de} }

\seealso{ \code{\link{trans} } }
\examples{
data(los.data)
my.observ <- prepare.los.data(x=los.data)
my.model <- msmodel(c("0","1","2","3"),cens.name="cens")
my.trans <- trans(model=my.model,observ=my.observ)
my.aj <- aj(my.trans,s=0,t=80)
}

\keyword{models}
