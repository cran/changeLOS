summary.clos <- function(object, ...) {
## ----------------------------------------------------------------------------
## Title: summary.clos
## ----------------------------------------------------------------------------
## Author: Matthias Wangler, mw@imbi.uni-freiburg.de
## Institute of Med. Biometry and Med. Computer Science
## Stefan-Meier-Strasse 26, D-79104 Freiburg,
## http://www.imbi.uni-freiburg.de
## ----------------------------------------------------------------------------
## Description: summary for an object of class 'clos'
## ----------------------------------------------------------------------------
## Required Packages: 
## ----------------------------------------------------------------------------
## Usage: summary.clos(object)
##
## object: an object of class 'clos'
## ----------------------------------------------------------------------------
## Value: an object of class 'summary.clos'. The object is a list of:
##
## cLOS: change in LOS
## patients: numeric, total number of observed patients (admissions)
## patients.discharge: numeric, number of patients being discharged
## patients.death: numeric, number of patients being death
## patients.cens: numeric, number of patients being censored
## cases: numeric, number of patients with complication
## cases.discharge: numeric, number of patients with complication being discharged
## cases.death: numeric, number of patients with complication being death
## cases.cens: numeric, number of patients with complication being censored
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
##          > summary(los)  
## ----------------------------------------------------------------------------
## License: GPL 2
##-----------------------------------------------------------------------------
## History: 03.08.2004, Matthias Wangler
##                      the first version
## ----------------------------------------------------------------------------
  if (!inherits(object, "clos")) {
    stop("Argument 'object' must be an object of class \"clos\".")
  }

  summary <- list( cLOS=object$cLOS, patients=object$patients, patients.discharge=object$patients.discharge,
                   patients.death=object$patients.death, patients.cens=object$patients.cens,
                   cases=object$cases, cases.discharge=object$cases.discharge, cases.death=object$cases.death,
                   cases.cens=object$cases.cens )

  class(summary) <- "summary.clos"
  
  summary
  
} ## end of function summary.clos
