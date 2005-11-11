"print.progdismodel" <- function(x, ...) {
## ----------------------------------------------------------------------------
## Title: print.progdismodel
## ----------------------------------------------------------------------------
## Author: Matthias Wangler, <mw@imbi.uni-freiburg.de>
## Institute of Med. Biometry and Med. Computer Science 
## Stefan-Meier-Strasse 26, D-79104 Freiburg 
## http://www.imbi.uni-freiburg.de 
## ----------------------------------------------------------------------------
## Description:  print function for an object of class 'progdismodel'
## ----------------------------------------------------------------------------
## Required Packages: 
## ----------------------------------------------------------------------------
## Usage: print.progdismodel(x)
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

  cat("The time points: \n")
  print(x$times.par)
  cat("\n")
  
  cat("P(death, t):\n")
  print(x$death)
  cat("\n")

  cat("P(death | risk factor absent, t):\n")
  print(x$death.given.rfa)
  cat("\n")

  cat("P(death | risk factor present, t):\n")
  print(x$death.given.rfp)
  cat("\n")

  cat("Attriburable Mortality: \n")
  print(x$AR)
  cat("\n")
  
  cat("Population Attriburable Mortality: \n")
  print(x$PAR)
  cat("\n")
  
} ## end of function print.progdismodel
