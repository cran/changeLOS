print.summary.clos <- function(x, ...) {
## ----------------------------------------------------------------------------
## Title: print.summary.clos
## ----------------------------------------------------------------------------
## Author: Matthias Wangler, mw@imbi.uni-freiburg.de
## Institute of Med. Biometry and Med. Computer Science
## Stefan-Meier-Strasse 26, D-79104 Freiburg,
## http://www.imbi.uni-freiburg.de
## ----------------------------------------------------------------------------
## Description: prints a summary for an object of class 'clos'
## ----------------------------------------------------------------------------
## Required Packages: 
## ----------------------------------------------------------------------------
## Usage: print.summary.clos(x)
##
## x: an object of class 'summary.clos'
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
##          > print(summary(los))
## ----------------------------------------------------------------------------
## License: GPL 2
##-----------------------------------------------------------------------------
## History: 03.08.2004, Matthias Wangler
##                      the first version
## ----------------------------------------------------------------------------
  if (!inherits(x, "summary.clos")) {
    stop("Argument 'x' must be an object of class \"summary.clos\".")
  }

  cat("\n")

  cat("Change in LOS: ", x$cLOS, sep="")
  cat("\n\n")
  
  cat("Change in LOS, patients discharged: ", x$phi2, sep="")
  cat("\n\n")
  
  cat("Change in LOS, patients deceased: ", x$phi3, sep="")
  cat("\n\n")
 
  a <- c("Number of observed patients", "Number of patients being discharged", "Number of patients who die",
         "Number of patients being censored", "Number of patients who experienced the intermediate event(IE)",
         "Number of patients who experienced the IE being discharged",
         "Number of patients who experienced the IE and died",
         "Number of patients who experienced the IE and were censored")
  
  b <- c(x$patients,x$patients.discharge,x$patients.death,
         x$patients.cens,x$cases,x$cases.discharge,
         x$cases.death, x$cases.cens)
  
  c <- c(round(100,2), round(x$patients.discharge/x$patients*100,2), round(x$patients.death/x$patients*100,2),
         round(x$patients.cens/x$patients*100,2), round(x$cases/x$patients*100,2),
         round(x$cases.discharge/x$patients*100,2),
         round(x$cases.death/x$patients*100,2), round(x$cases.cens/x$patients*100,2))

  colnames <- c("Total", "%")
  
  ms <- matrix(c(b,c), nrow = 8, ncol = 2, dimnames=list(a, colnames))

  print(ms)

  cat("\n")
  
} ## end of function print.summary.clos
