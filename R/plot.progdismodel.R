"plot.progdismodel" <-
function(x,
                                file1, file2, file3, 
                                lwd=2, cex=1.2,
                                lty1=1, lty2=1, lty3=1,
                                color1=4, color2=1, color3=2, ... ) {
## ----------------------------------------------------------------------------
## Title: plot.progdismodel
## ----------------------------------------------------------------------------
## Author: Matthias Wangler, <mw@imbi.uni-freiburg.de>
## Institute of Med. Biometry and Med. Computer Science 
## Stefan-Meier-Strasse 26, D-79104 Freiburg 
## http://www.imbi.uni-freiburg.de 
## ----------------------------------------------------------------------------
## Description: plot function for objects of class 'progdismodel'
## ----------------------------------------------------------------------------
## Required Packages: 
## ----------------------------------------------------------------------------
## Usage: plot.progdismodel(x,
##                          file1, file2, file3, 
##                          lwd=2, cex=1.2,
##                          lty1=1, lty2=2, lty3=3,
##                          color1=4, color2=1, color3=2, ... )
## ----------------------------------------------------------------------------
## Value:
## ----------------------------------------------------------------------------
## Notes: -
## ----------------------------------------------------------------------------
## Example:
## ----------------------------------------------------------------------------
## License: GPL 2
##-----------------------------------------------------------------------------
## History:
## ----------------------------------------------------------------------------
  if (!inherits(x, "progdismodel")) {
    stop("Argument 'x' must be an object of class \"progdismodel\".")
  }

  if( missing(file1) ) {
    ##x11()
    get(getOption("device"))()
  }
  else {
    postscript(file=file1)
  }          
  
  plot(y=x$death, x=x$times.par, type="l",
       xlab="Time t", ylab="mortality",
       main="Mortality",
       lwd=lwd,
       lty=lty1,
       col=color1,
       ylim=c(min(x$death, x$death.given.rfa , x$death.given.rfp),
              max(x$death, x$death.given.rfa , x$death.given.rfp)) )
  
  lines(y=x$death.given.rfa, x=x$times.par, type="l", lwd=lwd, lty=lty2, col=color2)
  lines(y=x$death.given.rfp, x=x$times.par, type="l", lwd=lwd, lty=lty3, col=color3)
  
  legend(x=max(x$times.par)/2,
         y=max(x$death, x$death.given.rfa , x$death.given.rfp)/4,
         legend=c("P( death )","P( death | risk factor absent)","P( death | risk factor present)") ,
         bty="o", pch=c(15,15,15), cex=cex,
         lty=c(lty1,lty2,lty3), col=c(color1,color2,color3))
  
  
  if( !missing(file1) ) {
    dev.off()
  }
  
  
  if( missing(file2) ) {    
    ##x11()
    get(getOption("device"))()
  }
  else {
    postscript(file=file2)
  }
  
  plot(y=x$AR, x=x$times.par, type="l",
       xlab="Time t", ylab="AR",
       main="Attributable Mortality",
       lwd=lwd,
       col=color1)
  
  lines(y=rep(0,length(x$times.par)), x=x$times.par, type="l",lty=2)
  
  if( !missing(file2) ) {
    dev.off()
  }
  
  if( missing(file3) ) {
    ##x11()
    get(getOption("device"))()
  }
  else {
    postscript(file=file3)
  }
  
  plot(y=x$PAR, x=x$times.par, type="l",
       xlab="Time t", ylab="PAR",
       main="Population Attributable Mortality",
       lwd=lwd,
       col=color1)
  
  lines(y=rep(0,length(x$times.par)), x=x$times.par, type="l",lty=2)
  
  if( !missing(file3) ) {
    dev.off()
  }
  
} ## end of function plot.progdismodel

