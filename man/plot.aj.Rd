\name{plot.aj}
\alias{plot.aj}
\title{  Plot function for objects of class 'aj'}
\description{
  Draws a plot for an object of class 'aj' in the active graphics
  device: Transition probabilities are plotted as a function of
  time. The time origin is determined by the object of the class \code{\link{aj}}. 
}
\usage{
\method{plot}{aj}(x,
        from,
        to,
        xlab=expression(paste(Time, " ", italic(t))),
        ylab= eval(substitute(expression(paste("Estimate of ",
                   P[{a}][{b}], "(", italic(s), ",",  italic(t), ")")),
                   list(a=from[1],b=to[1],s=x$start))),
        xlim = c(x$start,max(x$times)),
        ylim=c(0,1),
        lab=c(10,10,7),
        txt=eval(substitute(expression(paste(hat(P)[{a}][{b}],
                 "(",italic(s), ",", italic(t), ")")),
                 list(a=from[1],b=to[1],s=x$start))),
        x.txt=(xlim[2]+xlim[1])/2,
        y.txt=ylim[2]*0.9,
        col=1, ...)   
}
\arguments{
  \item{x}{ an object of the class 'aj'}
  \item{from}{ a character vector naming the states 'from'}
  \item{to}{ a character vector naming the states 'to'}
  \item{xlab}{ a title for the x axis}
  \item{ylab}{ a title for the y axis}
  \item{xlim}{ the x limits (min,max) of the plot}
  \item{ylim}{ the y limits (min,max) of the plot}
  \item{lab}{ A numerical vector of the form 'c(x, y, len)' which modifies
              the way that axes are annotated.  The values of 'x' and 'y'
              give the (approximate) number of tickmarks on the x and y
              axes and 'len' specifies the label size.}
  \item{txt}{ one or more character strings or expressions specifying a text to be written.}
  \item{x.txt}{ the x co-ordinates to be used to position the text}
  \item{y.txt}{ the y co-ordinate to be used to position the text}
  \item{col}{ the color of the line(s)}
  \item{...}{other graphical arguments}
}
\value{
  A matrix with ncol = 1 + length(from):
  \item{column 1:}{vector of the timepoints, the x-cooridinates}
  \item{column 2 to column ncol:}{vector of the estimated transition probabilities, the y-coordinates}  
}

\author{ Matthias Wangler \email{mw@imbi.uni-freiburg.de}}
\seealso{ \code{\link{aj}}}
\examples{
data(los.data)
my.observ <- prepare.los.data(x=los.data)

my.model <- msmodel(c("0","1","2","3"),cens.name="cens")
my.trans <- trans(model=my.model,observ=my.observ)
my.aj <- aj(my.trans, s=0, t=80)  
plot(my.aj,c("0","0","0","0"),c("0","1","2","3"))
}
\keyword{hplot}

