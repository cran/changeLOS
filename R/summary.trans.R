"summary.trans" <-
function(object, ...) {
## ----------------------------------------------------------------------------
## Title: summary.trans
## ----------------------------------------------------------------------------
## Author: Matthias Wangler, mw@imbi.uni-freiburg.de
## Institute of Med. Biometry and Med. Computer Science
## Stefan-Meier-Strasse 26, D-79104 Freiburg,
## http://www.imbi.uni-freiburg.de
## ----------------------------------------------------------------------------
## Description: summary for an object of class 'trans'
## ----------------------------------------------------------------------------
## Required Packages: 
## ----------------------------------------------------------------------------
## Usage: summary.trans(object)
##
## object: an object of class 'trans'
## ----------------------------------------------------------------------------
## Value: an object of class 'summary.trans'. The object is a list of:
##
## nrtransitions:  a matrix
##           column 1: state 'from'
##           column 2: state 'to'
##           column 3: number of transitions from state 'from' to state 'to'
##
## state.names: vector of the statenames
## ----------------------------------------------------------------------------
## Notes: -
## ----------------------------------------------------------------------------
## Example: > data(los.data)
##          > my.observ <- prepare.los.data(x=los.data)
##          > my.model <- msmodel(c("0","1","2","3"),cens.name="cens")
##          > my.trans <- trans(model=my.model,observ=my.observ)
##          > summary(my.trans)  
## ----------------------------------------------------------------------------
## License: GPL 2
##-----------------------------------------------------------------------------
## History: 03.08.2004, Matthias Wangler
##                      the first version
## ----------------------------------------------------------------------------
  if (!inherits(object, "trans")) {
    stop("Argument 'object' must be an object of class 'trans'")
  }

  summary <- list()

  summary$nrtransitions <- object$nrtransitions
  
  summary$state.names <- object$state.names

  summary$times <- object$times

  summary$nr.before <- object$nr.before
      
  class(summary) <- "summary.trans"

  summary
    
} ## end of function summary.trans

