"summary.progdismodel" <- function(object, ...) {
## ----------------------------------------------------------------------------
## Title: summary.progdismodel
## ----------------------------------------------------------------------------
## Author: Matthias Wangler, <mw@imbi.uni-freiburg.de>
## Institute of Med. Biometry and Med. Computer Science 
## Stefan-Meier-Strasse 26, D-79104 Freiburg 
## http://www.imbi.uni-freiburg.de 
## ----------------------------------------------------------------------------
## Description: summary function for objects of class 'progdismodel'
## ----------------------------------------------------------------------------
## Required Packages: 
## ----------------------------------------------------------------------------
## Usage: summary.progdismodel(object)
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
  if (!inherits(object, "progdismodel")) {
    stop("Argument 'object' must be an object of class \"progdismodel\".")
  }

  
  nt <- length(object$times.par)

  a <- object$death[nt]
  b <- object$death.given.rfa[nt]
  c <- object$death.given.rfp[nt]
  d <- object$AR[nt]
  e <- object$PAR[nt]

  s <- matrix( c(a,b,c,d,e), ncol=1)

  rownames(s) <- c("P(death)", "P(death | risk factor absent)",
                   "P(death | risk factor present)", "Attriburable Mortality",
                   "Population Attrributable Mortality")
    
  colnames(s) <- ""

  cat(paste("At time point ", object$times.par[nt], ": \n", sep=""))

  s

} ## end of function summary.progdismodel
