\name{plot.clos}
\alias{plot.clos}
\title{ Plot function for objects of class 'clos'}
\description{
  Draws a plot for an object of class 'clos' in the active graphics device:
  It illustrates the estimated expected change in LOS associated with an
  intermediate event.
}
\usage{
plot.clos(x,opt=0,  
          xlab=expression(paste(Time, " ", italic(t))),
          ylab.1="Expected LOS",ylab.2="Weights",
          xlim = c(0,max(x$trans$times[!is.na(x$e.given.1) | !is.na(x$e.given.0)])),
          ylim.1=c(0,max(x$e.given.1,x$e.given.0,na.rm=TRUE)),
          ylim.2=c(0,max(x$weights,na.rm=TRUE)),
          col1=c(1,1),col2=c(1),
          lab.1=c(10,10,7), lab.2=c(10,3,7),
          lgd=expression(paste(Intermediate, " ", event, " ", by, " ", time, " ",italic(t)),
                paste(No, " ", intermediate, " ", event, " ", by, " ", time, " ",italic(t))),
          x.lgd=0,y.lgd=ylim.1[2]*0.9, bty.lgd="n", ...)
}

\arguments{
  \item{x}{  an object of class 'clos'}
  \item{opt}{ numeric with the following valid values:
              \describe{
                \item{0}{plots weights and expected LOS (the default)}
                \item{1}{plots expected LOS}
                \item{2}{plots weights}
              }
	    }	    
  \item{xlab}{  a title for the x axis}
  \item{ylab.1}{ a title for the y axis in the plot of the expected LOS}
  \item{ylab.2}{ a title for the y axis in the plot of the weights}
  \item{xlim}{ the x limits (min,max) of the plot.}
  \item{ylim.1}{ the y limits of the plot of the expected LOS}
  \item{ylim.2}{ the y limits of the plot of the weights}
  \item{col1}{ the color of the lines in the plot of the expected LOS}
  \item{col2}{ the color of the line in the plot of the weights}
  \item{lab.1}{ A numerical vector of the form 'c(x, y, len)' which
    modifies the way that axes are annotated. The values of 'x' and 'y'
    give the (approximate) number of tickmarks on the x and y axes and
    'len' specifies the label size. Plot of the expected LOS}
  \item{lab.2}{ like lab.1, but for the plot of the weights}
  \item{lgd}{ a vector of text values or an 'expression' to appear in
    the legend of the plot of the expected LOS}
  \item{x.lgd}{ the x co-ordinate to be used to position the legend}
  \item{y.lgd}{ the y co-ordinate to be used to position the legend}
  \item{bty.lgd}{ the type of box to be drawn around the legend. The
    allowed values are 'n' (the default) and 'o'.}
  \item{...}{other graphical arguments}
}
\details{By default, two graphs are drawn. In a lower graph, the
  expected LOS given the intermediate event has (not yet) occurred by
  the time running on the x-axis is plotted. The difference between
  these two curves equals the change in LOS associated with the
  intermediate event \emph{by the time} running on the x-axis. I. e. the
  difference between the two curves is the estimator of E(LOS|\eqn{X_s =
  }intermediate event) - E(LOS|\eqn{X_s = }initial state), where where
  \eqn{X_s} denotes the state by time \eqn{s}.

  The estimated change in LOS is then computed as an weighted average
  over all these differences; the weights are illustrated in the upper
  plot. The weights derive from the estimated waiting time distribution
  in the initial state.
  
}
\author{Matthias Wangler \email{mw@imbi.uni-freiburg.de}}
\seealso{ \code{\link{clos}}}
\examples{
data(los.data)
my.observ <- prepare.los.data(x = los.data)
tra <- matrix(FALSE, 4, 4)
diag(tra) <- TRUE
tra[1, ] <- TRUE
tra[2, 3:4] <- TRUE
my.model <- msmodel(c("0", "1", "2", "3"), tra, cens.name = "cens")
los <- clos(model = my.model, observ = my.observ)
plot(los, xlim=c(0,80), ylim.1=c(0,120))
## same plot with nicer y-axis in the lower graph
plot(los, xlim=c(0,80), ylim.1=c(0,120), lab.1=c(8,6,7))
}

\keyword{hplot}
