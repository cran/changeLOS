"progdismodel" <-
function(model,observ, max.time) {
## --------------------------------------------------------------------------------
## Title: R-function progdismodel
## ---------------------------------------------------------------------------------
## Author: Matthias Wangler
##         mw@imbi.uni-freiburg.de
##         Institute of Med. Biometry and Med. Computer Science
##         Stefan-Meier-Strasse 26, D-79104 Freiburg,
##         http://www.imbi.uni-freiburg.de
## ---------------------------------------------------------------------------------
## Description: 
## ---------------------------------------------------------------------------------
## Required Packages: -
## ---------------------------------------------------------------------------------
## Usage: progdismodel(model,observ)
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
## max.time:
## ---------------------------------------------------------------------------------
## Value: 
## ---------------------------------------------------------------------------------
## License: GPL 2
##----------------------------------------------------------------------------------
## History:    27.09.2005, Matthias Wangler
##                         first version
## ---------------------------------------------------------------------------------

  ## argument 'model' must be passed
  if( missing(model) ) {
    stop("Argument 'model' is missing, with no defauls")
  }
  
  ## argument 'observ' must be passed
  if( missing(observ) ) {
    stop("Argument 'observ' is missing, with no defauls")
  }
  
  if( !inherits(model,"msmodel") ) {
    stop("Arguemnt 'model' must be an object of class 'msmodel'.") 
  }
      
  ## argument 'observ' must be data.frame
  if( !is.data.frame( observ ) ) {
    stop("Argument 'observ' is not data.frame")
  }

  ## check the model
  if( nrow(model$tra) != 6 ){
    stop("Argument 'model' must be a 'six-states model'.")
  }

  len <- length(model$state.names)    
  from <-  model$state.names[][model$transitions[,1][model$transition[,2] != len]]
  to <- model$state.names[][model$transitions[,2][model$transition[,2] != len]]
  ufrom <- unique(from)
  uto   <- unique(to)
  absorb   <- setdiff(uto,ufrom)

  if( length(ufrom) != 2 ) {
    stop("Argument 'model' must be a 'six-states model' with 2 transient states.")  
  }    
  if( length(absorb) != 4 ) {
    stop("Argument 'model' must be a 'six-states model' with 4 absorbing states.")  
  }
  
  ## check the number of columns of the passed data.frame observ
  if( dim(observ)[2] != 5 ){
    stop("The passed data.frame 'observ' doesn't include 5 columns.")      
  }

  ## check the column names of the passed data.frame observ
  if( names(observ)[1] != "id" || names(observ)[2] != "from" || names(observ)[3] != "to"
     || names(observ)[4] != "time"|| names(observ)[5] != "oid") {
    stop("The passed data.frame 'observ' must have the columns 'id', 'from', 'to', 'time' and 'oid'.")
  }
    
  ## check the number of rows of the passed data.frame observ
  if( dim(observ)[1] == 0 ) {
    stop("The passed data.frame 'observ' doesn't contain rows. There is nothing to do")
  }
     

  ## compute the Aalen-Johansen estimator for the matrices of transition posibilities
  my.trans <- trans(model, observ)
  my.aj <- aj(my.trans, s=0, t=max(my.trans$times))

  ## P(X(0)=0)
  PX0 <- my.trans$nr.before[1,1] / sum(my.trans$nr.before[1,])

  ## P(X(0)=1)
  PX1 <-  my.trans$nr.before[1,2] / sum(my.trans$nr.before[1,])

    
  p00 <- my.aj$matrices[1,1,]
  p01 <- my.aj$matrices[1,2,]
  p02 <- my.aj$matrices[1,3,]
  p03 <- my.aj$matrices[1,4,]
  p04 <- my.aj$matrices[1,5,]
  p05 <- my.aj$matrices[1,6,]
  p11 <- my.aj$matrices[2,2,]
  p14 <- my.aj$matrices[2,5,]
  p15 <- my.aj$matrices[2,6,]
  
  my.times.par <- my.aj$times

  if( !missing(max.time) ) {
    my.times.par <-  my.times.par[my.times.par <= max.time]
  }
  
  p00 <- p00[my.aj$times %in% my.times.par]
  p01 <- p01[my.aj$times %in% my.times.par]
  p02 <- p02[my.aj$times %in% my.times.par]
  p03 <- p03[my.aj$times %in% my.times.par]
  p04 <- p04[my.aj$times %in% my.times.par]
  p05 <- p05[my.aj$times %in% my.times.par]
  p11 <- p11[my.aj$times %in% my.times.par]
  p14 <- p14[my.aj$times %in% my.times.par]
  p15 <- p15[my.aj$times %in% my.times.par]
  
  prob.death <- PX0*(p03 + p05) + PX1*p15   ## death
  
  prob.death.rfa <- p03 / (p00+p02+p03) ## death, given risk factor absent
  prob.death.rfa[is.na(prob.death.rfa)] <- 0
  
  prob.death.rfp <- (PX0*p05 + PX1*p15) / (PX0*(p01+p04+p05) + PX1*(p11+p14+p15)) ## death, given risk factor present
  prob.death.rfp[is.na(prob.death.rfp)] <- 0

  my.ar <- prob.death.rfp - prob.death.rfa ## attributable risk
  
  my.par <- (prob.death - prob.death.rfa) / prob.death ## population attributable risk

  names(my.par) <- NULL
  names(my.ar) <- NULL
  names(prob.death) <- NULL
  names( prob.death.rfa) <- NULL
  names( prob.death.rfp) <- NULL

    
  ## return the result
  res <- list(trans=my.trans,
              aj=my.aj,
              times.par=my.times.par,
              PAR=my.par,
              AR=my.ar,
              death=prob.death,
              death.given.rfa=prob.death.rfa,
              death.given.rfp=prob.death.rfp)

  class(res) <- "progdismodel"

  return(res)
}

