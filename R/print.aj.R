print.aj <- function(x, ...) {
## ----------------------------------------------------------------------------
## Title: print.aj
## ----------------------------------------------------------------------------
## Author: Matthias Wangler, mw@imbi.uni-freiburg.de
## Institute of Med. Biometry and Med. Computer Science
## Stefan-Meier-Strasse 26, D-79104 Freiburg,
## http://www.imbi.uni-freiburg.de
## ----------------------------------------------------------------------------
## Description: prints an object of class 'aj'
## ----------------------------------------------------------------------------
## Required Packages: 
## ----------------------------------------------------------------------------
## Usage: print.aj(x)
##
## x: an object of class 'aj'
## ----------------------------------------------------------------------------
## Value:
## ----------------------------------------------------------------------------
## Notes: -
## ----------------------------------------------------------------------------
## Example: > data(los.data)
##          > my.observ <- prepare.los.data(x=los.data)
##          > my.model <- msmodel(c("0","1","2","3"),cens.name="cens")
##          > my.trans <- trans(model=my.model,observ=my.observ)
##          > my.aj <- aj(my.trans,s=0,t=80)  
##          > print(my.aj)
## ----------------------------------------------------------------------------
## License: GPL 2
##-----------------------------------------------------------------------------
## History: 03.08.2004, Matthias Wangler
##                      the first version
## ----------------------------------------------------------------------------
  if (!inherits(x, "aj")) {
    stop("Argument 'x' must be an object of class 'aj'")
  }

  summary(x)

  cat(paste("Time points in the interval (", x$start, ",", x$end, "]", sep=""))
  cat("\n")
  print(x$times)
  cat("\n")
  
  for(i in 1:length(x$times)) {
    cat(paste("Estimate of P(", x$start, ",", x$times[i], ")", sep=""))
    cat("\n")
    print(x$matrices[,,i])
    cat("\n")
  }
  
  cat("\n")
  
  
} ## end of function print.aj
