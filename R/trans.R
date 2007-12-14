"trans" <-
function(model, observ) {
## ----------------------------------------------------------------------------
## Title: trans
## ----------------------------------------------------------------------------
## Author: Matthias Wangler
##         mw@imbi.uni-freiburg.de
##         Institute of Med. Biometry and Med. Computer Science
##         Stefan-Meier-Strasse 26, D-79104 Freiburg,
##         http://www.imbi.uni-freiburg.de
## ----------------------------------------------------------------------------
## Description: Estimates a matrix of transition probabilities for every transition time 
## ----------------------------------------------------------------------------
## Required Packages: -
## ----------------------------------------------------------------------------
## Usage: trans(model=my.model,observ=my.observ)
##
## model:     an object of the class 'msmodel' which describes the multi-state model
##
## observ:   a data frame of the form data.frame( id, from, to, time )
##           id   : id (patient id, admision id, ...)
##           from : the state from where the transition occurs
##           to   : the state to whiche the transition occurs
##           time : the time the transition occurs
##           oid:  the observation id     
## ----------------------------------------------------------------------------
## Value: an object 'trans'
##
##        trans$matrices: array of transition matrices for every transition time
##        trans$times: the transition times
##        trans$nrtransitions: a matrix holding the number of transitions
##        trans$state.names: vector withe the names of the states
## ----------------------------------------------------------------------------  
## Notes: It's possible that the same patient, person or object was observed several
##        times (e.g. bootstrap).
##        So for each observation the same id recieves different observation id's. 
## ----------------------------------------------------------------------------
## Example: > data(los.data)
##          > my.observ <- prepare.los.data(x=los.data)
##          > my.model <- msmodel(c("0","1","2","3"),cens.name="cens")
##          > my.trans <- trans(model=my.model,observ=my.observ)  
## ----------------------------------------------------------------------------
## License: GPL 2
##-----------------------------------------------------------------------------
## History: 26.07.2004, Matthias Wangler
##                      the first version
## ----------------------------------------------------------------------------

  ## check the passed parameters
  if( missing(model) )
  {
    stop("Argument 'model' is missing, with no defauls.")
  }
  
  if( missing(observ) )
  {
    stop("Argument 'data' is missing, with no defauls.")
  }

  if( !inherits(model,"msmodel") )
  {
    stop("Arguemnt 'model' must be an object of class 'msmodel'.") 

  }

  if( !is.data.frame(observ) )
  {
    stop("Argument 'observ' must be a 'data.frame'.")
  }
    
  ## check the number of columns of the passed data.frame observ
  if( dim(observ)[2] != 5 )
  {
    stop("The passed data.frame 'observ' doesn't include 5 columns.")      
  }

  ## check the column names of the passed data.frame observ
  if( names (observ)[1] != "id" || names (observ)[2] != "from"
     || names(observ)[3] != "to" || names(observ)[4] != "time" || names(observ)[5] != "oid")
  {
    stop("The passed data.frame 'observ' must have the columns 'id', 'from', 'to', 'time' and 'oid'.")
  }
    
  ## check the number of rows of the passed data.frame observ
  if( dim(observ)[1] == 0 )
  {
    stop("The passed data.frame 'observ' doesn't contain rows. There is nothing to do")
  }
  
  states <- model$states - 1
  
  state.names <- model$state.names
  
  len <- length(model$state.names)
  
  ## check the observation data for undefined states
  stn <- unique( c(as.character(observ$from), as.character(observ$to) ) )
  if( length(stn[!(stn %in% state.names)]) > 0 ) {
    stop("Undefined states in the observation.")
  }
    
  ## check the observations if censoring TRUE or FALSE
  censoring <- FALSE
  if( state.names[len] %in% observ$to ) {
    censoring <- TRUE
  }
  
  ## check the obsrvation data for undefined transitions
  observ.transitions <- unique(observ[,2:3])    
  a <- paste(observ.transitions[,1],observ.transitions[,2],sep="")
  b <- paste( state.names[model$transitions[,1]], state.names[model$transitions[,2]],sep="")
  if( length(a[!(a %in% b)]) > 0 )
  {
    stop("Undefined transitions 'from'-'to' in the observation.")
  }

  observ <- observ[order(observ[,"time"]),]
  observ.time <- observ[,"time"]
  observ.from <- observ[,"from"]
  observ.to <- observ[,"to"]
    
  observ.from <- as.integer(factor(observ.from,levels=state.names)) - 1
  observ.from <- as.numeric(observ.from)
    
  observ.to <- as.integer(factor(observ.to, levels=state.names)) - 1
  observ.to <- as.numeric(observ.to)

  observ.event <- rep(1,length(observ.to))
  observ.event[observ.to==states[len]] <- 0
    
  nr.observ  <- length(observ.time)

  ## time points
  times <- sort(unique(observ.time[observ.to != states[len]]))

  ## array for the transition matrices
  matrices <- array(0, c( nrow(model$tra), ncol(model$tra),length(times)))

  nr.start <- rep(0, ifelse(censoring==TRUE,len,len-1))
  dd <- split(observ,observ$oid)
  count <- NULL
  for( i in 1:length(unique(observ$oid)) )
    {
      dd[[i]]$from <- as.factor(dd[[i]]$from)
      dd[[i]]$from2 <- as.numeric(levels(dd[[i]]$from)[dd[[i]]$from])
      count[i] <- dd[[i]][order(dd[[i]][,"time"])[1],"from2"]
      nr.start[count[i]+1] <- nr.start[count[i]+1]+1
    }
  
  ## matrix to store the number in each state just before the transition times (risk)
  nr.before <- matrix(0, nrow=length(times), ncol=ifelse(censoring==TRUE,len,len-1))
  nr.before[1,] <- nr.start

  ## matrix for storing the total number of transitions for all possible transitions
  nrtransitions <-  matrix(c(model$transitions[,1]-1, model$transitions[,2]-1, rep(0,length(model$transitions[,1]))),
                     nrow = length(model$transitions[,1]), ncol = 3, byrow = FALSE)       
	
  out <- .C("trans",
             as.integer(length(states)-1),
             as.integer(length(times)),
             as.integer(length(observ.time)),
             as.integer(model$transitions[,1]-1),
             as.integer(model$transitions[,2]-1),
             as.integer(nrow(model$transitions)),
             nj = as.integer(nrtransitions[,3]),
             as.double(observ.time),
             as.integer(observ.event),
             as.integer(observ.from),
             as.integer(observ.to),
             nr = as.double(nr.before),
             ma = as.double(matrices),
             PACKAGE="changeLOS")

  nrtransitions[,3] <- out$nj
  
  matrices <- array(out$ma, c( nrow(model$tra), ncol(model$tra), length(times)))

  dimnames(matrices) <- list(dimnames(model$tra)[[1]],dimnames(model$tra)[[2]] ,
                             paste("Time Point ",1:length(times), ": ", times, sep=""))

  nr.before <- matrix(out$nr, nrow=length(times), ncol=ifelse(censoring==TRUE,len,len-1))
  if(censoring==TRUE) co <- state.names[1:len]
  else                co <- state.names[1:(len-1)]
  colnames(nr.before) <- co
  rownames(nr.before) <- times

  res <- list(matrices, times, nrtransitions, model$state.names, nr.before)
  names(res) <- c("matrices", "times", "nrtransitions", "state.names", "nr.before")
  
  class(res) <- "trans"
  
  return(res)
}

