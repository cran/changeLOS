"prepare.progdismodel" <- function(data) {
## ----------------------------------------------------------------------------
## Title: prepare.progdismodel
## ----------------------------------------------------------------------------
## Author: Matthias Wangler, <mw@imbi.uni-freiburg.de>
## Institute of Med. Biometry and Med. Computer Science 
## Stefan-Meier-Strasse 26, D-79104 Freiburg 
## http://www.imbi.uni-freiburg.de 
## ----------------------------------------------------------------------------
## Description:
## ----------------------------------------------------------------------------
## Required Packages: 
## ----------------------------------------------------------------------------
## Usage: prepare.progdismodel(data)
##
## data: data.frame of the form data.frame( id, j.01, j.02, j.03, j.12, j.13, cens):
##
## id:      id (patient id, admision id, ...)
## j.01:    observed time for jump from "0" to "1"
## j.02:    observed time for jump from "0" to "2"
## j.03:    observed time for jump from "0" to "3"
## j.12:    observed time for jump from "1" to "2"
## j.13:    observed time for jump from "1" to "3"
## cens:    observed time for censoring
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

  ## prepare the data
  ## --> one row for each transition
  observ <- prepare.los.data(x = data)
  
  ## describe the model

  ## 6x6 matrix to describe the possible transitions
  tra <- rbind(c(TRUE,TRUE,TRUE,TRUE,FALSE,FALSE),
               c(FALSE,TRUE,FALSE,FALSE,TRUE,TRUE),
               c(FALSE,FALSE,TRUE,FALSE,FALSE,FALSE),
               c(FALSE,FALSE,FALSE,TRUE,FALSE,FALSE),
               c(FALSE,FALSE,FALSE,FALSE,TRUE,FALSE),
               c(FALSE,FALSE,FALSE,FALSE,FALSE,TRUE))
  
  model <- msmodel(c("0", "1", "2", "3", "4", "5"), tra, cens.name = "cens")

  
  ## prepare the data for the 'progressive disability model':

  observ$to <- factor(observ$to, levels=model$state.names)

  ## set 'to' = 4 if 'from' = 1 and 'to' = 2
  observ[,3][observ[,2]==1 & observ[,3]==2] <- 4

  ## set 'to' = 5 if 'from' = 1 and 'to' = 3
  observ[,3][observ[,2]==1 & observ[,3]==3] <- 5


  res <- list(observ=observ, model=model)

  return(res)
  
} ## end of function prepare.progdismodel
