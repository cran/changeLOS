plot.clos <- function(x,opt=0, xlab=expression(paste(Time, " ", italic(t))),
                      ylab.1="Expected LOS",ylab.2="Weights",
                      xlim = c(0,max(x$trans$times[!is.na(x$e.given.1) | !is.na(x$e.given.0)])),
                      ylim.1=c(0,max(x$e.given.1,x$e.given.0,na.rm=TRUE)),
                      ylim.2=c(0,max(x$weights,na.rm=TRUE)),
                      col1=c(1,1),col2=c(1),
                      lab.1=c(10,10,7), lab.2=c(10,3,7),
                      lgd=expression(paste(Intermediate, " ", event, " ", by, " ", time, " ",italic(t)),
                                      paste(No, " ", intermediate, " ", event, " ", by, " ", time, " ",italic(t))),
                      x.lgd=0,y.lgd=ylim.1[2]*0.9, bty.lgd="n", ...) {
## ----------------------------------------------------------------------------
## Title: plot.clos
## ----------------------------------------------------------------------------
## Author: Matthias Wangler, mw@imbi.uni-freiburg.de
## Institute of Med. Biometry and Med. Computer Science
## Stefan-Meier-Strasse 26, D-79104 Freiburg,
## http://www.imbi.uni-freiburg.de
## ----------------------------------------------------------------------------
## Description: plot function for an object of class 'clos'
## ----------------------------------------------------------------------------
## Required Packages: 
## ----------------------------------------------------------------------------
## Usage: plot.clos(x)
##
## x:      an object of class 'clos'
## opt:    0: plots weights and expected LOS
##         1: plots expected LOS
##         2: plots weights
## xlab:   a title for the x axis
## ylab.1: a title for the y axis in the plot of the expected LOS
## ylab.2: a title for the y axis in the plot of the weights
## xlim:   the x limits (min,max) of the plot.
## ylim.1: the y limits of the plot of the expected LOS
## ylim.2: the y limits of the plot of the wights
## col.1:  the line color of the plot of the expected LOS
## col.2:  the line color of the plot of the wights
## lab.1:   A numerical vector of the form 'c(x, y, len)' which modifies
##          the way that axes are annotated.  The values of 'x' and 'y'
##          give the (approximate) number of tickmarks on the x and y
##          axes and 'len' specifies the label size.
##          Plot of the expected LOS
## lab.2:   like lab.1, but for the plot of the weights
## lgd:     a vector of text values or an 'expression' to appear in the legend
##          of the plot of the expected LOS
## x.lgd,y.lgd: the x and y co-ordinates to be used to position the legend
## bty.lgd: the type of box to be drawn around the legend. The allowed values are
##          "n" (the default) and "o".
## ----------------------------------------------------------------------------
## Value: 
## ----------------------------------------------------------------------------
## Notes: -
## ----------------------------------------------------------------------------
## Example:  > data(los.data)
##           > my.observ <- prepare.los.data(x=los.data)  
##           > tra <- matrix(F, 4, 4)
##           > diag(tra) <- T
##           > tra[1, ] <- T
##           > tra[2, 3:4] <- T
##           > my.model <- msmodel(c("0", "1", "2", "3"), tra, cens.name = "cens")
##           > los <- clos(model=my.model,observ=my.observ)
##           > plot(los, xlim=c(0,80), ylim.1=c(0,120))
## ----------------------------------------------------------------------------
## License: GPL 2
##-----------------------------------------------------------------------------
## History: 09.08.2004, Matthias Wangler
##                      the first version
## ----------------------------------------------------------------------------
  if (!inherits(x, "clos")) {
    stop("Argument 'x' must be an object of class \"clos\".")
  }

  if(!(opt %in% c(0,1,2))) {
    stop("Allowed values for argument 'opt' are 0, 1 and 2.")
  }
    
  ## save default, for resetting
  def.par <- par(no.readonly = TRUE)
  on.exit(par(def.par))
  
  if( opt == 0 ) {   
    split.screen(figs=matrix(c(rep(0,2), rep(1,2), c(0, 0.6), c(0.7, 1)), ncol=4))
  }

  if( opt == 0 || opt == 2 ) {
    par(lab=lab.2)
    
    w.lab = xlab
    
    if( opt == 0 ) {
      screen(2)
      w.lab=""
      plot(x=c(0,x$w.times),y=c(0,x$weights), type ="s",axes=FALSE, lty=1, lwd=2, xlim=xlim,ylim=ylim.2,xlab=w.lab,ylab=ylab.2,col=col2)
      axis(side=2)
      box()
    }
    else {
      plot(x=c(0,x$w.times),y=c(0,x$weights), type ="s", lty=1, lwd=2, xlim=xlim,ylim=ylim.2,xlab=w.lab,ylab=ylab.2,col=col2)
    }
    
  }
  
  if( opt == 0 || opt == 1 ) {
    if( opt == 0 ) {
      screen(1)
    }
    
    par(lab=lab.1)
    
    plot(x=x$trans$times,y=x$e.given.1, type ="s", lty=1, lwd=2, xlim=xlim,ylim=ylim.1,xlab=xlab,ylab=ylab.1,col=col1[1])

    lines(x=x$trans$times,y=x$e.given.0, type ="s", lty=3, lwd=2,col=col1[2])

    legend(x.lgd, y.lgd, lgd , lty=c(1,3), bty=bty.lgd, col=col1)
    
  }  

  if( opt == 0 ) {
    close.screen(all = TRUE)
  }
} 
