summary.msmodel <- function(object, ...) {
## ----------------------------------------------------------------------------
## Title: summary.msmodel
## ----------------------------------------------------------------------------
## Author: Matthias Wangler, mw@imbi.uni-freiburg.de
## Institute of Med. Biometry and Med. Computer Science
## Stefan-Meier-Strasse 26, D-79104 Freiburg,
## http://www.imbi.uni-freiburg.de
## ----------------------------------------------------------------------------
## Description: summary for an object of class 'msmodel'
## ----------------------------------------------------------------------------
## Required Packages: 
## ----------------------------------------------------------------------------
## Usage: summary.msmodel(object, ...)
##
## object: an object of class 'msmodel'
## ----------------------------------------------------------------------------
## Value: 
## ----------------------------------------------------------------------------
## Notes: -
## ----------------------------------------------------------------------------
## Example:
##          > my.model <- msmodel(c("0","1","2","3"),cens.name="cens")
##          > summary(my.model)
## ----------------------------------------------------------------------------
## License: GPL 2
##-----------------------------------------------------------------------------
## History: 03.08.2004, Matthias Wangler
##                      the first version
## ----------------------------------------------------------------------------
  if (!inherits(object, "msmodel")) {
    stop("Argument 'object' must be an object of class 'msmodel'")
  }

  len <- length(object$state.names)
    
  cat("\n")
  cat(paste(nrow(object$tra), "-state model ",sep=""), "\n\n")

  cat("states (internal representation): ")
  cat(object$states[1:len-1], sep=", ")
  cat("\n\n")

  cat("names of the states: ")
  for(i in 1:(len-1)) {
    cat( paste("'", object$state.names[i], "'", if(i<len-1) ", ", sep=""))
  }
  cat("\n\n")
    
} ## end of function summary.msmodel
