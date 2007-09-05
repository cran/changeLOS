"print.msmodel" <-
function(x, ...) {
## ----------------------------------------------------------------------------
## Title: print.msmodel
## ----------------------------------------------------------------------------
## Author: Mustermann, <yourfault@somewhere.net>
## Institute of Med. Biometry and Med. Computer Science 
## Stefan-Meier-Strasse 26, D-79104 Freiburg 
## http://www.imbi.uni-freiburg.de 
## ----------------------------------------------------------------------------
## Description: Muster 
## ----------------------------------------------------------------------------
## Required Packages: 
## ----------------------------------------------------------------------------
## Usage: print.msmodel(x, ...)
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
  if (!inherits(x, "msmodel")) {
    stop("Argument 'x' must be an object of class 'msmodel'")
  }

  summary(x)

  len <- length(x$state.names)
    
  from <-  x$state.names[][x$transitions[,1][x$transition[,2] != len]]
  to <- x$state.names[][x$transitions[,2][x$transition[,2] != len]]

  ufrom <- unique(from)
  uto   <- unique(to)
  absorb   <- setdiff(uto,ufrom)
  
  cat("transient states: ")
  for(i in 1:length(ufrom)) {
    cat( paste("'", ufrom[i], "'", if(i<length(ufrom)) ", ", sep=""))
  }
  cat("\n\n")
  
  cat("absorbing states: ")
  if( length(absorb) == 0 ) {
    cat("none")
  }
  else {
    for(i in 1:length(absorb)) {
      cat( paste("'", absorb[i], "'", if(i<length(absorb)) ", ", sep=""))
    }
  }  
  cat("\n\n")
  
  fromtocens <-  x$state.names[][x$transitions[,1][x$transition[,2] == len]]
  cat("censored observations may occurr in states: ")
  for(i in 1:length(fromtocens)) {
    cat( paste("'", fromtocens[i], "'", if(i<length(fromtocens)) ", ", sep=""))
  }
  cat("\n\n")
    
  cat( paste("name of censoring code: '", x$state.names[len], "' " ,
             "(internal representation is: ",x$states[len],")", sep=""))
  cat("\n\n")

  cat("the possibles transitions:\n")
  fromto <- c("from", "to")
  nr <- 1:length(from)
  transitions <- matrix(c(from, to), nrow = length(from), ncol = 2, byrow = FALSE, dimnames=list(nr,fromto))  
  print(transitions)
  cat("\n")  
  
} ## end of function print.msmodel

