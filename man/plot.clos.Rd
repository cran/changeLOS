\name{plot.clos}
\alias{plot.clos}
\title{ Plot function for objects of class 'clos'}
\description{
  Draws a plot for an object of class 'clos' in the active graphics device:
  It illustrates the estimated expected change in LOS associated with an
  intermediate event.
}
\usage{
\method{plot}{clos}(x, opt=0,
          xlab=expression(paste(Time, " ", italic(t))),
          ylab.1="Expected LOS",ylab.2="Weights",
          xlim = c(0,max(x$trans$times[!is.na(x$e.given.1) |
                                       !is.na(x$e.given.0)])),
          xlim.2 = c(0,max(x$trans$times[!is.na(x$phi2.case) |
                                         !is.na(x$phi2.control) |
                                         !is.na(x$phi3.case) |
                                         !is.na(x$phi3.control)])),
          ylim.0=c(0,max(x$weights,na.rm=TRUE)),
          ylim.1=c(0,max(x$e.given.1,x$e.given.0,na.rm=TRUE)),
          ylim.2=c(0,max(x$x$phi2.case,x$phi2.control,
                         x$phi3.case,x$phi3.control, na.rm=TRUE)),
          col1=c(1,2),col2=c(1), lty1=c(1,1), lty2=c(1), lwd1=c(2,2), lwd2=c(2),
          lab.1=c(10,10,7), lab.2=c(10,3,7),
          lgd.1=expression(
              paste(Intermediate, " ", event, " ", by, " ", time, " ",italic(t)),
              paste(No, " ", intermediate, " ", event, " ", by, " ",
                    time, " ",italic(t))),
          lgd.2=expression(
              paste(Case, " ", term, " ", by, " ",  time, " ",
                    italic(t),", ",  patients, " ", discharged),
              paste(Control, " ", term, " ", by, " ", time, " ",
                    italic(t),", ",  patients, " ", discharged)),
          lgd.3=expression(
              paste(Case, " ", term, " ", by, " ", time, " "
                    ,italic(t),", ",  patients, " ", deceased),
              paste(Control, " ", term, " ", by, " ", time, " ",
                    italic(t),", ",  patients, " ", deceased)),             
          x.lgd=0,
          y.lgd.1=ylim.1[2]*0.9,
          y.lgd.2=ylim.2[2]*0.9,
          bty.lgd="n",
          cexlab=1,
          cexleg=1, ...)
}

\arguments{
  \item{x}{  an object of class 'clos'}
  \item{opt}{ numeric with the following valid values:
              \describe{
                \item{0}{plots weights and expected LOS (the default)}
                \item{1}{plots expected LOS}
                \item{2}{plots weights}
		\item{3}{plots phi2 and phi3}
		\item{4}{plots phi2}
		\item{5}{plots phi3}
              }
	    }	    
  \item{xlab}{  a title for the x axis}
  \item{ylab.1}{ a title for the y axis in the plot of the expected LOS}
  \item{ylab.2}{ a title for the y axis in the plot of the weights}
  \item{xlim}{ the x limits (min,max) of the plot of the expected LOS
    and of the weights}
  \item{xlim.2}{ the x limits (min,max) of the plot of phi2 and phi3}
  \item{ylim.0}{ the y limits of the plot of the weights}
  \item{ylim.1}{ the y limits of the plot of the expected LOS}
  \item{ylim.2}{ the y limits of the plot of phi2 and phi3}
  \item{col1}{ the color of the lines in the plot of the expected LOS,
    phi2 and phi3}
  \item{col2}{ the color of the line in the plot of the weights}
  \item{lty1}{ the line type of the lines in the plot of the expected
    LOS, phi2 and phi3}
  \item{lty2}{ the line type of the line in the plot of the weights}
  \item{lwd1}{the positive line width of the lines in the plot of the
    expected LOS, phi2 and phi3}
  \item{lwd2}{the positive line width of the line in the plot of the weights}  
  \item{lab.1}{ A numerical vector of the form 'c(x, y, len)' which
    modifies the way that axes are annotated. The values of 'x' and 'y'
    give the (approximate) number of tickmarks on the x and y axes and
    'len' specifies the label size. Plot of the expected LOS, phi2 and phi3}
  \item{lab.2}{ like lab.1, but for the plot of the weights}
  \item{lgd.1}{ a vector of text values or an 'expression' to appear in
    the legend of the plot of the expected LOS}
  \item{lgd.2}{ a vector of text values or an 'expression' to appear in
    the legend of the plot of phi2}
  \item{lgd.3}{ a vector of text values or an 'expression' to appear in
    the legend of the plot of phi3}  
  \item{x.lgd}{ the x co-ordinate to be used to position the legend of
    the plot of the expected LOS, phi2 and phi3}
  \item{y.lgd.1}{ the y co-ordinate to be used to position the legend of
    the the plot of the expected LOS}
  \item{y.lgd.2}{ the y co-ordinate to be used to position the legend of
    the the plot of phi2 and phi3}  
  \item{bty.lgd}{ the type of box to be drawn around the legend. The
    allowed values are 'n' (the default) and 'o'.}
  \item{cexlab}{The magnification to be used for x and y labels relative
    to the current.}
   \item{cexleg}{ character expansion factor used for the legend}
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
