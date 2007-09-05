"exampleclos" <-
function() {
## ----------------------------------------------------------------------------
## Title: exampleclos
## ----------------------------------------------------------------------------
## Author: Matthias Wangler
##         mw@imbi.uni-freiburg.de
##         Institute of Med. Biometry and Med. Computer Science
##         Stefan-Meier-Strasse 26, D-79104 Freiburg,
##         http://www.imbi.uni-freiburg.de  
## ----------------------------------------------------------------------------
## Description: An example program, to compute change in length of stay due to
##              an intermediate event using event history analysis.
## ----------------------------------------------------------------------------
## Required Packages: survival
## ----------------------------------------------------------------------------
## Usage: exampleclos()
## ----------------------------------------------------------------------------
## Value: list with the following arguments:
##
## cLOS:      change in LOS
## times:     time points t1, ...,tn
## e.given.1  estimates E(T|Xti = 1), i = 1,..,n
## e.given.0  estimates E(T|Xti = 0), i = 1,..,n
## weights    weights for the weighted average
## matrices   jump matrices (transition matrices)
## called     how  the function clos() was called
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
## Example:
## ----------------------------------------------------------------------------
## License: GPL 2
##-----------------------------------------------------------------------------
## History: 28.07.2004, Matthias Wangler
## ----------------------------------------------------------------------------
  data(los.data)
  my.observ <- prepare.los.data(x=los.data)
  tra <- matrix(F,4,4)
  diag(tra) <- T
  tra[1,] <- T
  tra[2,3:4] <- T  
  my.model <- msmodel(c("0","1","2","3"),tra,cens.name="cens")
  los <- clos(model=my.model,observ=my.observ)
  return(los)
} ## end of function exampleclos

