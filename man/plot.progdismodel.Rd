\name{plot.progdismodel}
\alias{plot.progdismodel}

\title{Plot function for objects of class 'progdismodel'}

\description{
  Draws three plots for an object of class 'progdismodel' in the active graphics
  devices or to the passed postscript files.
}

\usage{
plot.progdismodel(x,
                  file1, file2, file3, 
                  lwd=2, cex=1.2,
                  lty1=1, lty2=1, lty3=1,
                  color1=4, color2=1, color3=2, ... )
}

\arguments{
  \item{x}{ an object of the class 'aj'}
  \item{file1}{ps file name for the first plot (mortality)}
  \item{file2}{ps file name for the second plot(attributable mortality)}
  \item{file3}{ps file name for the third plot (population attributable
    mortality)}
  \item{lwd}{the line widths}
  \item{cex}{character expansion factor}
  \item{lty1}{the line type of line 1 of the first plot}
  \item{lty2}{the line type of line 2 of the second plot}
  \item{lty3}{the line type of line 3 of the third plot}
  \item{color1}{the color of line 1 of the first plot}
  \item{color2}{the color of line 2 of the second plot}
  \item{color3}{the color of line 3 of the third plot}
  \item{...}{other graphical arguments}
}

  
\details{
  The first plot shows \eqn{P( death, t )}, \eqn{P(death | risk factor absent,t)} and
  \eqn{P(death | risk factor present,t)}: mortality.

  The second plot shows \eqn{P(death | risk factor present,t) - P(death
    | risk factor absent,t)}: attributable mortality.

  The third plot shows \eqn{(P(death,t) - P(death | risk factor
    absent,t))/P(death,t)}: population attributable mortality.
  
}

\author{ Matthias Wangler \email{mw@imbi.uni-freiburg.de}}

\examples{
data(los.data)

p <- prepare.progdismodel(data=los.data)

pdm <- progdismodel(p$model, p$observ)

plot(pdm)
}

\keyword{hplot}
