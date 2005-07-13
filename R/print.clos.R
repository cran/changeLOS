print.clos <- function(x, ...) {
## ----------------------------------------------------------------------------
## Title: print.clos
## ----------------------------------------------------------------------------
## Author: Matthias Wangler, mw@imbi.uni-freiburg.de
## Institute of Med. Biometry and Med. Computer Science
## Stefan-Meier-Strasse 26, D-79104 Freiburg,
## http://www.imbi.uni-freiburg.de
## ----------------------------------------------------------------------------
## Description: print function for an object of class 'clos'
## ----------------------------------------------------------------------------
## Required Packages: 
## ----------------------------------------------------------------------------
## Usage: print.clos(x)
##
## x: an object of class 'clos'
## ----------------------------------------------------------------------------
## Value: 
## ----------------------------------------------------------------------------
## Notes: -
## ----------------------------------------------------------------------------
## Example: > data(los.data)
##          > my.observ <- prepare.los.data(x=los.data)
##          > trans <- matrix(F,4,4)
##          > diag(trans) <- T
##          > trans[1,] <- T
##          > trans[2,3:4] <- T  
##          > my.model <- msmodel(c("0","1","2","3"),trans,cens.name="cens")
##          > los <- clos(model=my.model,observ=my.observ)
##          > print(los)
## ----------------------------------------------------------------------------
## License: GPL 2
##-----------------------------------------------------------------------------
## History: 03.08.2004, Matthias Wangler
##                      the first version
## ----------------------------------------------------------------------------
  if (!inherits(x, "clos")) {
    stop("Argument 'x' must be an object of class \"clos\".")
  }

  print(summary(x))

  cat("\n")

  cat("Time points:\n")
  print(x$trans$times)
  cat("\n")

  cat("Estimates given in state 1:\n")
  print(x$e.given.1)
  cat("\n")

  cat("Estimates given in state 0:\n")
  print(x$e.given.0)
  cat("\n")  

  cat("Distinguishing  between patients discharged and patients deceased:\n")
  cat("Discharged, case term:\n")
  print(x$phi2.case)
  cat("\n")
  
  cat("Discharged, control term:\n")
  print(x$phi2.control)
  cat("\n")  

  cat("Deceased, case term:\n")
  print(x$phi3.case)
  cat("\n")
  
 cat("Deceased, control term:\n")
  print(x$phi3.control)
  cat("\n")   
  
  cat("the group `intermediate, but no terminal event yet'\n")
  cat("was empty for the following event times: ")
  print(x$empty.1)
  cat("\n")
  cat("the group `no intermediate or terminal event yet'\n")
  cat("was empty for the following event times: ")
  print(x$empty.0)
  cat("\n")
  cat("No comparison between groups was possible at these time points.\n")
  cat("Change in LOS associated with the intermediate event acquired \n")
  cat("up to such a time point was set to 0.\n")
  cat("\n")
  
  cat("Weights for the weighted average:\n")
  print(x$weights)
  cat("\n")    

  for(i in 1:length(x$trans$times)) {
    cat("Transition matrix for time point nr. ", i, ": ", x$trans$times[i], sep="")
    cat("\n")
    print(x$trans$matrices[,,i])
    cat("\n")
  }

  cat("The function 'clos' was called:\n")
  print(x$called)
  cat("\n")
  
} ## end of function print.clos
