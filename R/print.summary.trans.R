print.summary.trans <- function(x, ...) {
## ----------------------------------------------------------------------------
## Title: print.summary.trans
## ----------------------------------------------------------------------------
## Author: Matthias Wangler, mw@imbi.uni-freiburg.de
## Institute of Med. Biometry and Med. Computer Science
## Stefan-Meier-Strasse 26, D-79104 Freiburg,
## http://www.imbi.uni-freiburg.de
## ----------------------------------------------------------------------------
## Description: prints a summary for an object of class 'trans'
## ----------------------------------------------------------------------------
## Required Packages: 
## ----------------------------------------------------------------------------
## Usage: print.summary.trans(x, ...)
##
## x: an object of class 'summary.trans'
## ----------------------------------------------------------------------------
## Value: 
## ----------------------------------------------------------------------------
## Notes: -
## ----------------------------------------------------------------------------
## Example: > data(los.data)
##          > my.observ <- prepare.los.data(x=los.data)
##          > my.model <- msmodel(c("0","1","2","3"),cens.name="cens")
##          > my.trans <- trans(model=my.model,observ=my.observ)
##          > print(summary(my.trans))
## ----------------------------------------------------------------------------
## License: GPL 2
##-----------------------------------------------------------------------------
## History: 03.08.2004, Matthias Wangler
##                      the first version
## ----------------------------------------------------------------------------
  if (!inherits(x, "summary.trans")) {
    stop("Argument 'x' must be an object of class 'summary.trans'")
  }

  cat("\n")

  len <- length(x$state.names)
    
  cat("Total number of transitions:\n\n")

  from <- x$nrtransitions[,1]+1
  to <- x$nrtransitions[,2] + 1
  count <- x$nrtransitions[,3]

  mj <-  matrix( c( x$state.names[from[to!=len]],
                    x$state.names[to[to!=len]],
                    as.character(count[to!=len]) ),
                 nrow = length(x$state.names[from[to!=len]]),
                 ncol = 3, byrow = FALSE ) 
  
  prmatrix(mj, rowlab=rep("",nrow(mj)), collab=c("  from", "    to", "  transitions"), quote = FALSE, right = TRUE)
  cat("\n")
  
  if( length(count[to==len & count > 0]) > 0 ) {
    for( i in 1:length(from[to==len]) ) {
      cat(paste("censored in state ", x$state.names[from[to==len][i]], ": ",
                count[to==len][i], sep=""))
      cat("\n")
    }
  }
  
  cat("\n")

  cat("the initial distribution:\n")
  for(i in 1:(len-1)) {
    cat(x$nr.before[1,i], " in state '", as.character(x$state.names[i]), "'", sep= "")
    cat("\n")
  }
  cat("\n")

  cat("the number in the states just before the transition times:\n")
  m <- cbind(x$times, x$nr.before)
  colnames(m) <- c("time", colnames(x$nr.before))
  rownames(m) <- rep("",nrow(m))
  print(m)

  cat("\n")
 
  
  
} ## end of function print.summary.trans
