\name{clos}
\alias{clos}
\title{ Change in LOS }
%\encoding{latin1}
\description{
  estimates the expected change in length of stay (LOS) associated with
  an intermediate event (IE). In order to account for the timing of events (an
  IE can only have an effect once it has occurred), a
  multi-state model is used. The data are stratified into `cases' (the
  IE has occurred) and `controls' (the IE has
  \emph{not yet} occurred) on a daywise basis. (`Daywise' if time is in
  days.) The expected change in LOS is estimated for each day and a
  weighted average is computed.
}
\usage{
clos(model, observ, aw=FALSE)
}
\arguments{
  \item{model}{ an object of the class 'msmodel' which describes the
    multi-state model }
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
  \item{aw}{locical, with (TRUE) or without(FALSE) alternative weighting}	  
	 
}
\details{
  \code{clos} is based on `approach B' in Schulgen and Schumacher
  (1996), however with some modifications: We use a multi-state model
  with four states. All individuals start in an initial state. They may
  pass through the intermediate state. LOS is determined by reaching one
  of two competing, absorbing states. In clinical research these latter
  two states are typically `discharge' and `death',
  respectively. The IE is often interpreted as a complication. Backward
  transitions are not possible. Schulgen and
  Schumacher
  (1996) modelled a second intermediate state, motivated by their
  application. However, this model has some individuals `drop out' in
  the sense that they become neither `case' nor `control'. In addition,
  \code{clos} computes an weighted average,
  where the weights are given by the weighting time distribution in the
  intial state. Schulgen and Schumacher (1996) used a conditional version
  of it, given one reaches the IE. In the case of no censoring, the
  former weights have every individual contribute to the weighting,
  whereas the latter has not. One can interpret the conditional
  weighting in Schulgen and Schumacher (1996) as assuming a
  patient's viewpoint who experiences the IE (and analogously for
  patients uninfected,
  with weights given one is directly discharged/dies without prior IE.)
  Finally, \code{clos} treats the daywise
  comparison as being zero if there are not both patients with and
  without IE on that day.
  See also Beyersmann et al. (2005) for a non-technical explanation of
  these methods.
  
  See the examples for special features:
  \itemize{
    \item bootstrapping
    \item use of the alternative weight
    \item distinguishing between patients discharged and patients deceased    
    }
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
  \item{phi2}{weighted average of the difference between phi2.case and
    phi2.control, this quantity can be interpreted as the contribution
    to the expected change in LOS at time \eqn{s} by patients infected
    at time \eqn{s} who eventually discharge
  }  
  \item{phi2.case}{estimates E(LOS \strong{1}(\eqn{X_LOS = }discharge)|\eqn{X_s = }intermediate event),
    where \strong{1} denotes the indicator function}
  \item{phi2.control}{P(\eqn{X_LOS = }discharge|\eqn{X_s = }intermediate
    event) E(LOS|\eqn{X_s = }initial state)} 
  \item{phi3}{weighted average of the difference between phi3.case and
    phi3.control, this quantity can be interpreted as the contribution
    to the expected change in LOS at time \eqn{s} by patients infected
    at time \eqn{s} who eventually die
  }
  \item{phi3.case}{estimates E(LOS \strong{1}(\eqn{X_LOS = }death)|\eqn{X_s = }intermediate event),
    where \strong{1} denotes the indicator function}
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
\section{Value 'alternative weighting'}{
  When aw=TRUE an object of class c('clos','closa'). The object is a
  list like above with the following extra items:

  \describe{
  \item{weights.1}{weights corresponding to the conditional weighting
    time in the intial state given one reaches the IE.}
  \item{weights.23}{weights corresponding to the conditional weighting
    time in the intial state given one \emph{does not} reach the IE.}
  \item{given.1}{change in LOS corresponding to the alternative weights 'weights.1'}
  \item{given.23}{change in LOS corresponding to the alternative weights
    'weights.2'}
  }
}
\note{required packages: survival }
\references{
  G Schulgen and M Schumacher (1996). Estimation of prolongation of
  hospital stay attributable to nosocomial infections.  \emph{Lifetime
  Data Analysis} 2, 219-240.

  J Beyersmann, P Gastmeier, H Grundmann, S Baerwolf, C Geffers,
  M Behnke, H Rueden, and M Schumacher (2006). Use of Multistate
  Models to Assess Prolongation of Intensive Care Unit Stay Due to
  Nosocomial Infection. \emph{Infection Control and Hospital
    Epidemiology} 27, 493-499.

  M Wrangler, J Beyersmann and M Schumacher (2006). changeLOS: An
  R-package for change in length of hospital stay based on the
  Aalen-Johansen estimator. \emph{R News} 6(2), 31--35.
  
  A reference on multi-state models:

  P Andersen and N Keiding (2002). Multi-state models for event history
  analysis. \emph{Statistical Methods in Medical Research} 11, 91-115}
\author{
  Jan Beyersmann \email{jan@fdm.uni-freiburg.de}
  Matthias Wangler \email{mw@imbi.uni-freiburg.de}
}
\seealso{ \code{\link{survival}}, \link{exampleclos} }
\examples{

## run clos
data(los.data)
my.observ <- prepare.los.data(x=los.data)
tra <- matrix(FALSE,4,4)
diag(tra) <- TRUE
tra[1,] <- TRUE
tra[2,3:4] <- TRUE
my.model <- msmodel(c("0","1","2","3"),tra,cens.name="cens")
los <- clos(model=my.model,observ=my.observ)
summary(los)
plot(los)

## expected change in LOS due to an IE
phi <- los$e.given.1 - los$e.given.0

## distinguishing between patients discharged
phi2 <- los$phi2.case - los$phi2.control
## and patients deceased
phi3 <- los$phi3.case - los$phi3.control 

## we have phi = phi2 + phi3
all(round(phi, digits=10) == round(phi2+phi3, digits=10),na.rm=TRUE)


## compute bootstrap SE with function boot() from library boot

## first we need a statistic, which boot takes as an argument
"clos.for.bstrap" <- function(data, index, mod){
  my.observ <- prepare.los.data(x=data[index,])
  return(clos(model=mod,observ=my.observ)$cLOS)
}

## our estimate is
clos.for.bstrap(data=los.data, index=1:length(los.data[,1]), mod=my.model)

## now bootstrap 
library(boot)
nb <- 20 ## nb <- 2000
los.bootstrap <- boot(los.data, clos.for.bstrap, nb, mod=my.model)
sqrt(var(los.bootstrap$t))


## compute change in LOS `by hand'; also works for alternative weights

## first, get pure event times
my.evtimes <- sort(unique(my.observ$time[my.observ$to !="cens"]))
## compute daywise difference between expected LOS
my.diffs <- los$e.given.1 - los$e.given.0
## restrict to those days when a transition out of the initial state was observed
my.diffs <- my.diffs[is.element(my.evtimes, los$w.times)]
## compute weighted average, but don't sum up over
## days where one of the groups (no) IE (yet) experienced was empty.
sum((my.diffs * los$weights)[!(is.element(los$w.times, c(los$empty.0, los$empty.1)))])

}

\keyword{survival}
\keyword{models}
