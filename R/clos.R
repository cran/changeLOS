"clos" <- function(model,observ,aw=FALSE) {
## ---------------------------------------------------------------------------------
## Title: R-function clos()
## ---------------------------------------------------------------------------------
## Author: Jan Beyersmann
##         jan@fdm.uni-freiburg.de
##         1.)  Institute of Med. Biometry and Med. Computer Science
##              Stefan-Meier-Strasse 26, D-79104 Freiburg,
##              http://www.imbi.uni-freiburg.de
##         2.)  Freiburg Centre for Data Analysis and Modelling
##              Eckerstrasse 1, D-79104 Freiburg
##              http://www.fdm.uni-freiburg.de
##
##         Matthias Wangler
##         mw@imbi.uni-freiburg.de
##         Institute of Med. Biometry and Med. Computer Science
##         Stefan-Meier-Strasse 26, D-79104 Freiburg,
##         http://www.imbi.uni-freiburg.de  
## --------------------------------------------------------------------------------
## Description: Compute change in length of stay due to an intermediate event using
##              event history analysis.
## ---------------------------------------------------------------------------------
## Required Packages: survival
## ---------------------------------------------------------------------------------
## Usage: clos(model,observ)
##
## model:    an object of the class 'msmodel' which describes the multi-state model
##
## observ:   data.frame with of the form data.frame( id, from, to, time ):
##
##           id:   id (patient id, admision id ...)
##           from: the state from where a transition occurs
##           to:   the state to which a transition occurs
##           time: the time a transition occurs
##           oid:  the observation id   
##
## aw.       locical, alternative weighting TRUE/FALSE
## ---------------------------------------------------------------------------------
## Value: an object of class 'clos'. The object is a list of:
##
## cLOS:      change in LOS
## trans:     an object of class 'trans'
## e.given.1  estimates E(T|Xti = 1), i = 1,..,n
## e.given.0  estimates E(T|Xti = 0), i = 1,..,n
## empty.1    event times: the group `intermediate, but no terminal event yet' was empty
## empty.0    event times: the group `no intermediate or terminal event yet' was empty
## weights    weights for the weighted average
## called     how  the function clos() was called
## patients: numeric, total number of observed patients (admissions)
## patients.discharge: numeric, number of patients being discharged
## patients.death: numeric, number of patients being death
## patients.cens: numeric, number of patients being censored
## cases: numeric, number of patients with complication
## cases.discharge: numeric, number of patients with complication being discharged
## cases.death: numeric, number of patients with complication being death
## cases.cens: numeric, number of patients with complication being censored  
## --------------------------------------------------------------------------------
## Notes: It's possible that the same patient, person or object was observed several
##        times (e.g. bootstrap).
##        So for each observation the same id recieves different observation id's. 
## --------------------------------------------------------------------------------
## Example: > data(los.data)
##          > my.observ <- prepare.los.data(x=los.data)
##          > trans <- matrix(F,4,4)
##          > diag(trans) <- T
##          > trans[1,] <- T
##          > trans[2,3:4] <- T  
##          > my.model <- msmodel(c("0","1","2","3"),trans,cens.name="cens")
##          > los <- clos(model=my.model,observ=my.observ)
## ---------------------------------------------------------------------------------
## License: GPL 2
##----------------------------------------------------------------------------------
## History:   01/10/2003, Jan Beyersmann
##                        the first version
##            20/06/2004, Matthias Wangler
##                        1. error handling: checking the passed arguments
##                        2. no fixed column names in the passed data frame
##                        3. function for computing transition matrices for any finite set of states
## ---------------------------------------------------------------------------------
  
    ## save the call of the function
    fcall <- sys.call()
    
    ## check arguments
    if( nargs() > 3 || nargs() < 2 )
    {
      stop("Invalid number of arguments")
    }
    else
    {
      ## argument 'model' must be passed
      if( missing(model) )
      {
        stop("Argument 'model' is missing, with no defauls")
      }
      
      ## argument 'observ' must be passed
      if( missing(observ) )
      {
        stop("Argument 'observ' is missing, with no defauls")
      }
      
      if( !inherits(model,"msmodel") )
      {
        stop("Arguemnt 'model' must be an object of class 'msmodel'.") 
      }
      
      ## argument 'observ' must be data.frame
      if( !is.data.frame( observ ) )
      {
        stop("Argument 'observ' is not data.frame")
      }
    }

    ## check the model
    if( nrow(model$tra) != 4 )
    {
      stop("Argument 'model' must be a 'four-sate model'.")
    }
    len <- length(model$state.names)    
    from <-  model$state.names[][model$transitions[,1][model$transition[,2] != len]]
    to <- model$state.names[][model$transitions[,2][model$transition[,2] != len]]
    ufrom <- unique(from)
    uto   <- unique(to)
    absorb   <- setdiff(uto,ufrom)
    if( length(ufrom) != 2 ) {
       stop("Argument 'model' must be a 'four-sate model' with 2 transient states.")  
    }    
    if( length(absorb) != 2 ) {
      stop("Argument 'model' must be a 'four-sate model' with 2 absorbing states.")  
    }
    
    ## check the number of columns of the passed data.frame observ
    if( dim(observ)[2] != 5 )
    {
      stop("The passed data.frame 'observ' doesn't include 5 columns.")      
    }

    ## check the column names of the passed data.frame observ
    if( names(observ)[1] != "id" || names(observ)[2] != "from" || names(observ)[3] != "to"
       || names(observ)[4] != "time"|| names(observ)[5] != "oid")
    {
      stop("The passed data.frame 'observ' must have the columns 'id', 'from', 'to', 'time' and 'oid'.")
    }
    
    ## check the number of rows of the passed data.frame observ
    if( dim(observ)[1] == 0 )
    {
      stop("The passed data.frame 'observ' doesn't contain rows. There is nothing to do")
    }
    
    ## need package survival
    require(survival)

    ## name of the censoring variable
    cens.state <- model$state.names[length(model$states)]
    
    ## compute the transition times and the `transition matrix' for every transition time
    my.trans <-  trans(model, observ)

    if( length(unique(observ$oid)) != my.trans$nr.start[1] ) {
      stop(paste("Not all individuals start in the initial state '", model$state.names[1], "'.", sep=""))
    }
    
    ## compute expected LOS given the state at every observed transition time _except for_
    ## the greatest observed time (which may be a censoring time)

    ## is there a censoring time greater than the last observed transition time?
    my.times <- sort(unique(c(my.trans$times, max(observ$time[observ$to==cens.state], my.trans$times))))
    
    ## find `transition matrices' that corresponds to my.times[i + 1]...
    my.matrices <- array(0, c( nrow(model$tra), ncol(model$tra), length(my.times)))
    
    if( length(my.times) > 2 ) {     
      for(i in (length(my.times)-2):1)
        {
          my.matrices[,,i] <- my.trans$matrices[,,length(my.trans$times[my.trans$times <= my.times[i + 1]])]
        }    
    }
    
    los <- matrix(data=rep(my.times,3), ncol=3, byrow=FALSE, dimnames=list(NULL, c("Time", "Given in state 1", "Given in state 0")))

    dim <- dim(my.matrices)[1]
  
    ## will need to temporarily store Aalen-Johansen estimates
    aj <- array(diag(1,dim,dim), c(dim, dim,1))    
        
    ## will need function that does matrix multiplication running
    ## thru the `slices' of array aj
    "my.function" <- function(x,y){ x%*%aj[,,y] }
    
    if( length(my.times) > 2 ) {
      ## last two rows in los already correct. compute the rest, starting with the 
      ## last but two row (which corresponds to the third greatest time in my.times)
      los[length(my.times)-1,2:3] <- rep(max(my.times), 2)
      
      for(i in (length(my.times)-2):1) {
        ## compute time differences
        diffs <- diff(my.times[(i+1):length(my.times)])
        
        ## multiply `transition matrix' with Aalen-Johansen estimates of the previous loop
        aj <- array(apply(X=diag(1:dim(aj)[3]), 1, my.function, x=my.matrices[,,i]), c(dim ,dim, dim(aj)[3]))
        
        ## LOS given in state 1 at time my.times[i]
        los[,2][i] <- my.times[i+1] + matrix(diffs, nrow=1) %*% matrix(aj[2,2,],ncol=1)
        
        ## LOS given in state 0 at time my.times[i]
        los[,3][i] <- my.times[i+1] + matrix(diffs, nrow=1) %*% matrix((aj[1,1,] + aj[1,2,]),ncol=1)
        
        ## stack identity matrix on top for the next loop
        aj <- array(c(diag(1, dim, dim), aj), c(dim, dim, (dim(aj)[3] + 1)))
      }
    }
    
    ## if one of the two groups of patients ('cases' and 'controls') is empty,
    ## the expected change in LOS at the questioned time must be set to 'NA'
    times.empty.0 <- c()
    times.empty.1 <- c()
    for(i in 1:length(my.times))
    {
      group.0 <- sum(observ$from == model$state.names[1], na.rm=TRUE) -
                 sum(observ$from ==  model$state.names[1] & observ$time <= my.times[i], na.rm=TRUE)
      
      group.1 <- sum(observ$to ==  model$state.names[2] & observ$time <= my.times[i] , na.rm=TRUE) -
                 sum(observ$from ==  model$state.names[2] & observ$time <= my.times[i], na.rm=TRUE)

      if( group.0 == 0 )
      {
        los[,3][i] <- NA
        times.empty.0 <- c( times.empty.0, my.times[i] )
      }

      if( group.1 == 0 )
      {
        los[,2][i] <- NA
        times.empty.1 <- c( times.empty.1, my.times[i] )
      }
    }
    
    ## compute distribution to weight differences in LOS
    ## need waiting time distribution in initial state 0.
    ## create a survival object and fit it. event: left state 0.
    T0 <- Surv(observ$time[observ$from== model$state.names[1]], 1 - (observ$to==cens.state)[observ$from== model$state.names[1]])
    T0.fit <- survfit(T0)
    ## only need `time' and `surv' for non-censoring events
    T0.fit$time <- T0.fit$time[T0.fit$n.event!=0]
    T0.fit$surv <- T0.fit$surv[T0.fit$n.event!=0]

    ## weight by waiting time distribution in initial state 0
    ## need mass for each event time
    ## my.weights[i] corresponds to T0.fit$time[i]
    my.weights <- diff(c(0,1-T0.fit$surv))

    ## compute estimate
    los.diff <- los[,2]-los[,3]
    los.diff[is.na(los.diff)] <- 0    
    estimate <- matrix(los.diff[is.element(los[,1], T0.fit$time)], nrow=1)  %*% matrix(my.weights, ncol=1)    
    
    ## some results for a summary:
    
    ## number of patients
    my.patients <- sum(my.trans$nrtransitions[,3][my.trans$nrtransitions[,1]==1])
    ## number of patients being discharged
    my.patients.discharge <- sum(my.trans$nrtransitions[,3][my.trans$nrtransitions[,2]==3])
    ## number of patients being death
    my.patients.death <- sum(my.trans$nrtransitions[,3][my.trans$nrtransitions[,2]==4])
    ## number of patients being censored
    my.patients.cens <- sum(my.trans$nrtransitions[,3][my.trans$nrtransitions[,2]==5]) 

    ## number of patients with complication
    my.cases <- sum(my.trans$nrtransitions[,3][my.trans$nrtransitions[,2]==2]) 
    ## number of patients with complication being discharged
    my.cases.discharge <- sum(my.trans$nrtransitions[,3][my.trans$nrtransitions[,1]==2 & my.trans$nrtransitions[,2]==3])
    ## number of patients with complication being death
    my.cases.death <- sum(my.trans$nrtransitions[,3][my.trans$nrtransitions[,1]==2 & my.trans$nrtransitions[,2]==4])
    ## number of patients with complication being censored
    my.cases.cens <- sum(my.trans$nrtransitions[,3][my.trans$nrtransitions[,1]==2 & my.trans$nrtransitions[,2]==5])
        
    ## return results
    res <- list(cLOS=estimate, trans=my.trans, e.given.1=c(los[,2]), e.given.0=c(los[,3]),
                w.times=T0.fit$time,weights=my.weights,called=fcall,
                patients=my.patients, patients.discharge=my.patients.discharge,
                patients.death=my.patients.death, patients.cens=my.patients.cens,
                cases=my.cases, cases.discharge=my.cases.discharge,
                cases.death=my.cases.death, cases.cens=my.cases.cens,
                empty.0=times.empty.0, empty.1=times.empty.1)

    class(res) <- "clos"

    if( aw == FALSE ) {
      return(res)
    }
    
    ## alternative weighting

    ## 2nd: waiting time distribution in initial state 0 given cause for leaving is state 1

    pr.cause1 <- matrix(c(1, T0.fit$surv[1:length(T0.fit$surv)-1]), nrow=1) %*%
      matrix(my.trans$matrices[1,2,][is.element(my.trans$times, T0.fit$time)], ncol=1)

    ## my.weights.1[i] corresponds to T0.fit$time[i]
    my.weights.1 <- diag(diag(c(1, T0.fit$surv[1:length(T0.fit$surv)-1])) %*%
                         diag(my.trans$matrices[1,2,][is.element(my.trans$times, T0.fit$time)])) / pr.cause1

    ##estimate.1 <- matrix((los[,2]-los[,3])[is.element(los[,1], T0.fit$time)], nrow=1)  %*% matrix(my.weights.1, ncol=1)
    estimate.1 <- matrix(los.diff[is.element(los[,1], T0.fit$time)], nrow=1)  %*% matrix(my.weights.1, ncol=1)
    
    ## 3rd: waiting time distribution in initial state 0 given cause for leaving are states 2 or 3

    pr.cause23 <- matrix(c(1, T0.fit$surv[1:length(T0.fit$surv)-1]), nrow=1) %*%
      matrix((my.trans$matrices[1,3,] + my.trans$matrices[1,4,])[is.element(my.trans$times, T0.fit$time)], ncol=1)

    ## my.weights.23[i] corresponds to T0.fit$time[i]
    my.weights.23 <- diag(diag(c(1, T0.fit$surv[1:length(T0.fit$surv)-1])) %*%
                          diag((my.trans$matrices[1,3,] + my.trans$matrices[1,4,])[is.element(my.trans$times, T0.fit$time)])) / pr.cause23

    ##estimate.23 <- matrix((los[,2]-los[,3])[is.element(los[,1], T0.fit$time)], nrow=1)  %*% matrix(my.weights.23, ncol=1)
    estimate.23 <- matrix(los.diff[is.element(los[,1], T0.fit$time)], nrow=1)  %*% matrix(my.weights.23, ncol=1)

    ## return results
    res$weights.1=my.weights.1
    res$weights.23=my.weights.23
    res$given.1=estimate.1
    res$given.23=estimate.23

    class(res) <- c( "clos", "closa")

    return(res)      
    
  }## end of function
