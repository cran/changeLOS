"msmodel" <-
function(state.names,tra=matrix(TRUE,length(state.names),length(state.names)),cens.name="cens") {
## ----------------------------------------------------------------------------
## Title: multi-state model
## ----------------------------------------------------------------------------
## Author: Matthias Wangler
##         mw@imbi.uni-freiburg.de
##         Institute of Med. Biometry and Med. Computer Science
##         Stefan-Meier-Strasse 26, D-79104 Freiburg,
##         http://www.imbi.uni-freiburg.de
## ----------------------------------------------------------------------------
## Description: Makes a 'msmodel' - object for describing a 'multi-sate model'
## ----------------------------------------------------------------------------
## Required Packages: 
## ----------------------------------------------------------------------------
## Usage: msmodel(state.names,tra,cens.name)
##
##        state.names: character vector of the statenames
##  
##        tra: matrix of locigal values, which describes the posibles transitions
##
##        cens.name: character string (name of the censoring variable)
## ----------------------------------------------------------------------------
## Value: an object 'msmodel':
##
##        msmodel$tra: quadratic matrix of locigal values. TRUE/FALSE: transition is/is not possible
##        msmodel$states: numeric vector of the states  
##        msmodel$state.names: character vector of the statenames
##        msmodel$transitions: matrix with two columns, 1.column: state 'from', 2.column: state 'to'
##                             the number of rows is the number of possibles transitions
## ----------------------------------------------------------------------------
## Notes: The following example, a four-state model, is used in 'clos':
##
##        msmodel$tra:         0       1        2        3
##                         0   TRUE    TRUE     TRUE     TRUE
##                         1   FALSE   TRUE     TRUE     TRUE
##                         2   FALSE   FALSE    TRUE     FALSE
##                         3   FALSE   FALSE    FALSE    TRUE
##
##        msmodel$states:      1 (no complication),         2 (complication), 3 (discharge), 4 (death), 5(censoring state)
##
##        msmodel$state.names: 0 (no complication),         1 (complication), 2 (discharge), 3 (death), cens (censoring state)
##
##        msmodel$transitions: from   to
##                             1      2
##                             1      3
##                             1      4
##                             1      5
##                             2      3
##                             2      4
##                             2      5
##
## ----------------------------------------------------------------------------
## Example: > my.model <- msmodel(state.names=c("0","1","2","3"))
## ----------------------------------------------------------------------------
## License: GPL 2
##-----------------------------------------------------------------------------
## History: 26.07.2004, Matthias Wangler
##                      the firs version
## ----------------------------------------------------------------------------

  if( !is.character(state.names)) {
    stop("Argument 'state.names' must be character vector.")
  }
     
  if(!is.matrix(tra) ) {
    stop("'tra' must be a matrix, which describes the possible transitions.")
  }
    
  if(nrow(tra) != ncol(tra)) {
    stop("Argument 'tra' must be a quadratic  matrix.")    
  }

  if( nrow(tra) != length(state.names) ) {
    stop("The row number of 'tra' must be equal to the number of states.")
  }
    
  if(!is.logical(tra)) {
    stop("'tra' must be a matrix of logical values, which describes the possible transitions.")
  }
    
  if( length(state.names) != length(unique(state.names)) ) {
    stop("The statenames must be unique.")
  }

  if( cens.name %in% state.names ) {
    stop("The name of the censoring variable just is a name of the model states.")
  }  
  
  states <- 1:length(state.names)

  dimnames(tra) <- list(state.names,state.names)
     
  transitions.from <- row(tra)[tra[,]==TRUE][row(tra)[tra[,]==TRUE] != col(tra)[tra[,]==TRUE]]
  transitions.to   <- col(tra)[tra[,]==TRUE][row(tra)[tra[,]==TRUE] != col(tra)[tra[,]==TRUE]]

  states = c(states, length(states)+1)
  state.names = c(state.names, cens.name)
    
  transitions.to   <- c(transitions.to, rep(length(states),length(unique(transitions.from))))
  transitions.from <- c(transitions.from,unique(transitions.from))

  transitions <-  matrix(c(transitions.from,transitions.to), nrow = length(transitions.from), ncol = 2, byrow = FALSE)

  transitions <- transitions[order(transitions[,1], transitions[,2]),]

  msm <- list(tra, states, state.names, transitions)

  names(msm) <- c("tra", "states", "state.names", "transitions")

  class(msm) <- "msmodel"

  return(msm)
  
} ## end of function msmodel

