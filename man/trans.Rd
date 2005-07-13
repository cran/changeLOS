\name{trans}
\alias{trans}
\title{ transition matrices}
\description{
  computes the Aalen-Johansen estimator for the matrix of transition
  probabilities \eqn{P(u-, u)} for all observed transition times \eqn{u}. The
  entry \eqn{(l,m)} of the matrix
  denotes the estimated probability that state \eqn{m} has been
  reached by time \eqn{u} given state \eqn{l} has been occupied
  just before time \eqn{u}.
}
\usage{
trans(model, observ)
}
\arguments{
  \item{model}{ an object of the class 'msmodel' which describes the multi-state model }
  \item{observ}{ a data.frame of the form data.frame( id, from, to, time, oid )
    (see also \code{\link{prepare.los.data}}):
            \describe{
                \item{id}{id (patient id, admision id)}
                \item{from}{the state from where a transition occurs}
                \item{to}{the state to which a transition occurs}
                \item{time}{the time a transition occurs}
		\item{oid}{the observation id}		
            }    
  }
}
\details{
  The estimator for \eqn{P(u-, u)} is described by Andersen et
  al. (1993) at the bottom of p. 288. Non-diagonal entries \eqn{(h, j)}
  are given as the number of observed transitions from state\eqn{h}
  to state \eqn{j}, divided by the number of individuals in state
  \eqn{h} just
  prior to time \eqn{u}. The diagonal elements are
  chosen such that the sum of each row equals 1.

  The Aalen-Johansen estimator for \eqn{P(s, t)} can then be computed as
  matrix product of all matrices \eqn{P(u-, u)} for all transition times
  \eqn{u} in \eqn{(s,t]}, see \code{\link{aj}}.

}
\value{
  An object of the class 'trans'. The object is a list of:
  \item{matrices}{ array of matrices \eqn{P(u-, u)} for every transition
    time \eqn{u}}
  \item{times}{ the transition times}
  \item{nrtransitions}{ a matrix with
                           \describe{
                           \item{column 1}{the state from where a transition occurs}
                           \item{column 2}{the state to which a transition occurs}
                           \item{column 3}{the number of transitions}
                           }
                        }  
  \item{state.names}{ vector with the names of the states}
  \item{nr.before}{ matrix with the number in each state just before the transition times}
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
\author{ Matthias Wangler \email{mw@imbi.uni-freiburg.de}}

\seealso{ \code{\link{msmodel}, \link{clos}, \link{aj}}}
\examples{
data(los.data)
my.observ <- prepare.los.data(x=los.data)
my.model <- msmodel(c("0","1","2","3"),cens.name="cens")
my.trans <- trans(model=my.model,observ=my.observ)
my.matrices <- my.trans$matrices
my.times <- my.trans$times
}

\keyword{models}
