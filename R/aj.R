"aj" <-
function(tr,s=0,t=tr$times[length(tr$times)]) {
## ----------------------------------------------------------------------------
## Title:  Aalen-Johansen estimator
## ----------------------------------------------------------------------------
## Author: Matthias Wangler
##         mw@imbi.uni-freiburg.de
##         Institute of Med. Biometry and Med. Computer Science
##         Stefan-Meier-Strasse 26, D-79104 Freiburg,
##         http://www.imbi.uni-freiburg.de
## ----------------------------------------------------------------------------
## Description: computes the Aalen-Johansen estimator for the transition matrix
##              P(s,t). The estimator is a finite matrix product, one matrix for
##              every observed event time in the time interval (s,t].
## ----------------------------------------------------------------------------
## Required Packages: 
## ----------------------------------------------------------------------------
## Usage: aj(tr, s, t)
##
## tr: an object of 'trans'
## s:  begin of the time interval
## t:  end of the time interval
## ----------------------------------------------------------------------------
## Value: An object of class 'aj'. The object is a list of:
##
## matrix:   estimator for the transition matrix P(s,t)
## start:    the begin s of the time interval (s,t]
## end:      the end t of the time interval (s,t]
## times:    the jumptimes in the interval (s,t]
## matrices: array of estimators for P(s,j) for all jumptimes j in (s,t]
## ----------------------------------------------------------------------------
## Notes: -
## ----------------------------------------------------------------------------
## Example: > data(los.data)
##          > my.observ <- prepare.los.data(x=los.data)
##          > my.model <- msmodel(c("0","1","2","3"),cens.name="cens")
##          > my.trans <- trans(model=my.model,observ=my.observ)
##          > my.aj <- aj(my.trans,s=0,t=80) 
## ----------------------------------------------------------------------------
## License: GPL 2
##-----------------------------------------------------------------------------
## History: 02.08.2004, Matthias Wangler
##                      the first version
## ----------------------------------------------------------------------------

    if( !inherits(tr,"trans") )
    {
      stop("Arguemnt 'tr' must be an object of class 'trans'.") 
    }

    if( !is.numeric(t) )
    {
      stop("Argument 't' must be numeric.")
    }
    
    if( !is.numeric(s) )
    {
      stop("Argument 's' must be numeric.")
    }
    
    if( s >= t )
    {
      stop("Time 's' must be before time 't'.")
    }

    if( s >=  tr$times[length(tr$times)] )
    {
      stop("Invalid time 's'.")
    }

    if( t < tr$times[1])
    {
      stop("Invalid time 't'.")
    }    
    
    dim <- dim(tr$matrices)[1]

    first <- length(tr$times[tr$times <= s]) + 1
    last <- length(tr$times[tr$times <= t])

    if(!all((first:last) %in% (1:dim(tr$matrices)[3])))
    {
      stop("Not for all time points exists a transition matrix.")
    }
    
    ajm <- diag(1,dim,dim)

    ajt <- array(diag(1,dim,dim), c(dim, dim,1))  

    if( first <= last )
    {
      j <- 1
      
      for(i in first:last) {
        ajm <- ajm %*% tr$matrices[,,i]

       ajt[,,j] <- ajm

       j <- j + 1

       if( i < last ) {
         ajt <- array(c(ajt,diag(1, dim, dim)), c(dim, dim, (dim(ajt)[3] + 1)))
       }
      }
    }

    dimnames(ajm) <- rep(dimnames(tr$matrices)[1],2)
    
    dimnames(ajt) <- list(dimnames(tr$matrices)[[1]], dimnames(tr$matrices)[[1]], paste("Estimate of P(", s, ",", tr$times[first:last], ")", sep=""))

    res <- list(matrix=ajm,start=s,end=t,times=tr$times[first:last],matrices=ajt)
    
    class(res) <- "aj"
      
    return(res)

} ## end of function aj

