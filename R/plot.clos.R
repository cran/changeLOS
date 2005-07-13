plot.clos <- function(x,opt=0,
                      xlab=expression(paste(Time, " ", italic(t))),
                      ylab.1="Expected LOS",ylab.2="Weights",
                      xlim = c(0,max(x$trans$times[!is.na(x$e.given.1) | !is.na(x$e.given.0)])),
                      xlim.2 = c(0,max(x$trans$times[!is.na(x$phi2.case) | !is.na(x$phi2.control) |
                                                     !is.na(x$phi3.case) | !is.na(x$phi3.control) ])),                     
                      ylim.0=c(0,max(x$weights,na.rm=TRUE)),
                      ylim.1=c(0,max(x$e.given.1,x$e.given.0,na.rm=TRUE)),
                      ylim.2=c(0,max(x$x$phi2.case,x$phi2.control,x$phi3.case,x$phi3.control, na.rm=TRUE)),                
                      col1=c(1,2),col2=c(1), lty1=c(1,1), lty2=c(1), lwd1=c(2,2), lwd2=c(2),
                      lab.1=c(10,10,7), lab.2=c(10,3,7),
                      lgd.1=expression(paste(Intermediate, " ", event, " ", by, " ", time, " ",italic(t)),
                          paste(No, " ", intermediate, " ", event, " ", by, " ", time, " ",italic(t))),
                      lgd.2=expression(
                          paste(Case, " ", term, " ", by, " ",  time, " ",italic(t),", ",  patients, " ", discharged),
                          paste(Control, " ", term, " ", by, " ", time, " ",italic(t),", ",  patients, " ", discharged)),
                      lgd.3=expression(
                          paste(Case, " ", term, " ", by, " ", time, " ",italic(t),", ",  patients, " ", deceased),
                          paste(Control, " ", term, " ", by, " ", time, " ",italic(t),", ",  patients, " ", deceased)),             
                      x.lgd=0,
                      y.lgd.1=ylim.1[2]*0.9,
                      y.lgd.2=ylim.2[2]*0.9,
                      bty.lgd="n",
                      cexlab=1,
                      cexleg=1, ...) {
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
## xlim.2:   
## ylim.0: the y limits of the plot of the wights 
## ylim.1: the y limits of the plot of the expected LOS
## ylim.2: 
## col.1:  the line color of the plot of the expected LOS
## col.2:  the line color of the plot of the wights
## lab.1:   A numerical vector of the form 'c(x, y, len)' which modifies
##          the way that axes are annotated.  The values of 'x' and 'y'
##          give the (approximate) number of tickmarks on the x and y
##          axes and 'len' specifies the label size.
##          Plot of the expected LOS
## lab.2:   like lab.1, but for the plot of the weights
## lgd.1:   a vector of text values or an 'expression' to appear in the legend
##          of the plot of the expected LOS
## lgd.2:
## lgd.3:
## x.lgd:  the x co-ordinate to be used to position the legends
## y.lgd.1: the y co-ordinate to be used to position the legend of the plot of the expected LOS
## y.lgd.2
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

  if(!(opt %in% 0:5)) {
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
      op <- par(mar=c(2, 5, 2, 1))
      w.lab=""
      plot(x=c(0,x$w.times),y=c(0,x$weights), type ="s",axes=FALSE, lty=lty2, lwd=lwd2, xlim=xlim,
           ylim=ylim.0,xlab=w.lab,ylab=ylab.2,col=col2, cex.lab=cexlab, ...)
      axis(side=2)
      box()
      par(op)
    }
    else {
      plot(x=c(0,x$w.times),y=c(0,x$weights), type ="s", lty=lty2, lwd=lwd2, xlim=xlim,ylim=ylim.0,
           xlab=w.lab,ylab=ylab.2,col=col2, ...)
    }    
  }
  
  if( opt == 0 || opt == 1 ) {
    if( opt == 0 ) {            
      screen(1)      
      op <- par(mar=c(5, 5, 4, 1))  
    }
    
    par(lab=lab.1)
    
    plot(x=x$trans$times,y=x$e.given.1, type ="s", lty=lty1[1], lwd=lwd1[1], xlim=xlim,ylim=ylim.1,
         xlab=xlab,ylab=ylab.1,col=col1[1], cex.lab=cexlab, ...)

    lines(x=x$trans$times,y=x$e.given.0, type ="s", lty=lty1[2], lwd=lwd1[2],col=col1[2], ...)

    legend(x.lgd, y.lgd.1, legend=lgd.1 , lty=lty1, bty=bty.lgd, col=col1, lwd=lwd1,cex=cexleg, ...)
    
  }  



  if(opt==3) {
    split.screen(c(2,1))
  }
  
  if(opt==3 || opt==4) {
    ax <- TRUE
    
    w.lab = xlab
    
    if(opt==3) {
      screen(1)
      ##op <- par(mar=c(5, 5, 4, 1))
      op <- par(mar=c(1, 5, 2, 1))
      w.lab=""
      ax <- FALSE
    }
    
    par(lab=lab.1)
    
    plot(x=x$trans$times,y=x$phi2.case, type ="s", axes=ax, lty=lty1[1], lwd=lwd1[1], xlim=xlim.2,ylim=ylim.2,
         xlab=w.lab,ylab=ylab.1,col=col1[1], cex.lab=cexlab, ...)
    axis(side=2)
    box()
    lines(x=x$trans$times,y=x$phi2.control, type ="s", lty=lty1[2], lwd=lwd1[2],col=col1[2], ...)
    legend(x.lgd, y.lgd.2, legend=lgd.2 , lty=lty1, bty=bty.lgd, col=col1, lwd=lwd1,cex=cexleg, ...)

  }
      
    
  if(opt==3 || opt==5) {
    if(opt==3) {
      screen(2)

      op <- par(mar=c(5, 5, 1, 1))  
    }

    par(lab=lab.1)
    
    plot(x=x$trans$times,y=x$phi3.case, type ="s", lty=lty1[1], lwd=lwd1[1], xlim=xlim.2,ylim=ylim.2,
         xlab=xlab,ylab=ylab.1,col=col1[1], cex.lab=cexlab, ...)

    lines(x=x$trans$times,y=x$phi3.control, type ="s", lty=lty1[2], lwd=lwd1[2],col=col1[2], ...)

    legend(x.lgd, y.lgd.2, legend=lgd.3 , lty=lty1, bty=bty.lgd, col=col1, lwd=lwd1,cex=cexleg, ...)
    
  }

 if( opt == 0 || opt==3 ) {
    par(op)
    close.screen(all = TRUE)
  } 
} 
