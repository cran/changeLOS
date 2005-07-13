\name{progdismodel}
\alias{progdismodel}

\title{progressive disability model}

\description{
  Change in LOS and impact of an intermediate event on mortality can
  also be investigated in a so-called progressive disability model.
  This multi-state model can be described and the Aalen-Johansen
  estimator for transition probabilities can be computed.
  Further will be computed and plotted the mortality, attributable mortality
  and the population attributable mortality.
}
\usage{
progdismodel(data, plot=FALSE, file1, file2, file3, max.time,
             lwd=2, cex=1.2, lty1=1, lty2=2, lty3=3,
             color1=4, color2=1, color3=2)
}

\arguments{
  \item{data}{data.frame of the form data.frame(id, j.01, j.02, j.03, j.12, j.13, cens):
    \describe{
    \item{id:}{id (patient id, admision id)}      
    \item{j.01:}{observed time for jump from 0 to 1}
    \item{j.02:}{observed time for jump from 0 to 2}
    \item{j.03:}{observed time for jump from 0 to 3}
    \item{j.12:}{observed time for jump from 1 to 2}
    \item{j.13:}{observed time for jump from 1 to 3}
    \item{cens:}{censoring time (either in initial or intermediate state)}
    }
  }
  \item{plot}{TRUE/FALSE}
  \item{file1}{ps file name for the first plot (mortality)}
  \item{file2}{ps file name for the second plot(attributable mortality)}
  \item{file3}{ps file name for the third plot (population attributable
    mortality)}
  \item{max.time}{the last time point of interest}
  \item{lwd}{the line widths}
  \item{cex}{character expansion factor}
  \item{lty1}{the line type of line 1 of the first plot}
  \item{lty2}{the line type of line 2 of the second plot}
  \item{lty3}{the line type of line 3 of the third plot}
  \item{color1}{the color of line 1 of the first plot}
  \item{color2}{the color of line 2 of the second plot}
  \item{color3}{the color of line 3 of the third plot}
}

\details{
  The first plot shows \eqn{P( death, t )}, \eqn{P(death | risk factor absent,t)} and
  \eqn{P(death | risk factor present,t)}: mortality.

  The second plot shows \eqn{P(death | risk factor present,t) - P(death
    | risk factor absent,t)}: attributable mortality.

  The third plot shows \eqn{(P(death,t) - P(death | risk factor
    absent,t))/P(death,t)}: population attributable mortality.
  
}

\value{
  A list of
  \item{trans}{an object of class \code{\link{trans}}}
  \item{aj}{an object of class \code{\link{aj}}}
  \item{times.par}{}
  \item{PAR}{population attributable mortality}
  \item{AR}{attributable mortality}
  \item{death}{\eqn{P(death,t)}}
  \item{death.given.rfa}{\eqn{P(death | risk factor absent,t)}}
  \item{death.given.rfp}{\eqn{P(death | risk factor present,t)}}
  \item{tab}{table with the values at the last time point. if max.time
    is not missing, this will be the last time point.}}
}

\references{Matthias Wangler \email{mw@imbi.uni-freiburg.de}}


\examples{
\dontrun{
 data(los.data)
 res <- progdismodel(los.data)
}
}

\keyword{models}
