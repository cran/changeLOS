print.trans <- function(x, ...) {
## ----------------------------------------------------------------------------
## Title: print.trans
## ----------------------------------------------------------------------------
## Author: Matthias Wangler, mw@imbi.uni-freiburg.de
## Institute of Med. Biometry and Med. Computer Science
## Stefan-Meier-Strasse 26, D-79104 Freiburg,
## http://www.imbi.uni-freiburg.de
## ----------------------------------------------------------------------------
## Description: prints an object of class 'trans'
## ----------------------------------------------------------------------------
## Required Packages: 
## ----------------------------------------------------------------------------
## Usage: print.trans(x)
##
## x: an object of class 'trans'
## ----------------------------------------------------------------------------
## Value:
## ----------------------------------------------------------------------------
## Notes: -
## ----------------------------------------------------------------------------
## Example: > data(los.data)
##          > my.observ <- prepare.los.data(x=los.data)
##          > nr.0 <- length(my.observ$from[my.observ$from=="0"])
##          > my.model <- msmodel(c("0","1","2","3"),c(nr.0,0,0,0),cens=TRUE,cens.state.name="cens")
##          > my.trans <- trans(model=my.model,observ=my.observ)
##          > print(my.trans)
## ----------------------------------------------------------------------------
## License: GPL 2
##-----------------------------------------------------------------------------
## History: 03.08.2004, Matthias Wangler
##                      the first version
## ----------------------------------------------------------------------------
  if (!inherits(x, "trans")) {
    stop("Argument 'x' must be an object of class 'trans'")
  }

  cat("\n")

  print(summary(x))

  cat("Transition matrices:\n\n")
  
  for(i in 1:length(x$times)) {
    cat("Estimated transition matrix for P(", x$times[i], "-,",x$times[i],"):", sep="")
    cat("\n")
    print(x$matrices[,,i])
    cat("\n")
  }
  
  cat("\n")
  
} ## end of function print.trans
