summary.aj <- function(object, ...) {
## ----------------------------------------------------------------------------
## Title: summary.aj
## ----------------------------------------------------------------------------
## Author: Matthias Wangler, mw@imbi.uni-freiburg.de
## Institute of Med. Biometry and Med. Computer Science
## Stefan-Meier-Strasse 26, D-79104 Freiburg,
## http://www.imbi.uni-freiburg.de
## ----------------------------------------------------------------------------
## Description: summary for an object of class 'aj'
## ----------------------------------------------------------------------------
## Required Packages: 
## ----------------------------------------------------------------------------
## Usage: summary(object)
##
## object: an object of class 'aj'
## ----------------------------------------------------------------------------
## Value: x$matrix
## ----------------------------------------------------------------------------
## Notes: -
## ----------------------------------------------------------------------------
## Example: > data(los.data)
##          > my.observ <- prepare.los.data(x=los.data)
##          > my.model <- msmodel(c("0","1","2","3"),cens.name="cens")
##          > my.trans <- trans(model=my.model,observ=my.observ)
##          > my.aj <- aj(my.trans,s=0,t=80)  
##          > summary(my.aj)
## ----------------------------------------------------------------------------
## License: GPL 2
##-----------------------------------------------------------------------------
## History: 03.08.2004, Matthias Wangler
##                      the first version
## ----------------------------------------------------------------------------
  if (!inherits(object, "aj")) {
    stop("Argument 'object' must be an object of class 'aj'")
  }

  cat("\n")
  
  cat(paste("The Aalen-Johansen estimator for the transition matrix P(", object$start, ",", object$end, ")", sep=""))
  cat("\n")
  
  print(object$matrix)

  cat("\n")
    
} ## end of function print.aj
