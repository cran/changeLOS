\name{progdismodel}
\alias{progdismodel}

\title{progressive disability model}

\description{
  Change in LOS and impact of an intermediate event on mortality can
  also be investigated in a so-called progressive disability model.
  This multi-state model can be described and the Aalen-Johansen
  estimator for transition probabilities can be computed.
  Further will be computed the mortality, attributable mortality
  and the population attributable mortality.
}
\usage{
progdismodel(model,observ, max.time)
}

\arguments{
  \item{observ}{a data.frame of the form data.frame(id, from, to, time,
    oid):
    \describe{
      \item{id:}{ id (patient id, admision id)}
      \item{from:}{ the state from where a transition occurs}
      \item{to:}{ the state to which a transition occurs}
      \item{time:}{ time of the transition}
      \item{oid:}{ the observation id}
    }
  }
  \item{model}{an object of class 'msmodel', which describes the
    progressive disability model with two transient and four absorbing states.}
  \item{max.time}{the last time point of interest}
}

\details{
  mortality:
  \eqn{P( death, t )}, \eqn{P(death | risk factor absent,t)} and
  \eqn{P(death | risk factor present,t)}.

  attributable mortality:
  \eqn{P(death | risk factor present,t) - P(death | risk factor absent,t)}.

  population attributable mortality:
  \eqn{(P(death,t) - P(death | risk factor absent,t))/P(death,t)}. 
}

\value{
  A list of
  \item{trans}{an object of class \code{\link{trans}}}
  \item{aj}{an object of class \code{\link{aj}}}
  \item{times.par}{the time points for which the following values are computed.}
  \item{PAR}{population attributable mortality}
  \item{AR}{attributable mortality}
  \item{death}{\eqn{P(death,t)}}
  \item{death.given.rfa}{\eqn{P(death | risk factor absent,t)}}
  \item{death.given.rfp}{\eqn{P(death | risk factor present,t)}}
}

\author{Matthias Wangler \email{mw@imbi.uni-freiburg.de}}

\examples{
data(los.data)

p <- prepare.progdismodel(data=los.data)

pdm <- progdismodel(p$model, p$observ)
}

\keyword{models}
