## ------------------------------------------------------------------------------------------------------------
## Title: Datensatz SIR 3-Studie
## ------------------------------------------------------------------------------------------------------------
## Author: Matthias Wangler
## Institute of Med. Biometry and Med. Computer Science
## Stefan-Meier-Strasse 26, D-79104 Freiburg,
## http://www.imbi.uni-freiburg.de
## ------------------------------------------------------------------------------------------------------------
## Description: Einlesen des Datensatzes der SIR 3-Studie
## ------------------------------------------------------------------------------------------------------------
## Required Packages: -
## ------------------------------------------------------------------------------------------------------------
## Usage: los.data
## ------------------------------------------------------------------------------------------------------------
## Value:
## ------------------------------------------------------------------------------------------------------------
## Notes: -
## ------------------------------------------------------------------------------------------------------------
## Example:
## ------------------------------------------------------------------------------------------------------------
## License: GPL 2
##-------------------------------------------------------------------------------------------------------------
## History:
## ------------------------------------------------------------------------------------------------------------
  ##file <- "los.data.txt"
  
  ##sep <- ";"
 
  ##header <- TRUE

  ##row.names <- NULL

  ##pos.id <- 1

  ##pos.columns <- c(2,3,4,5,6,7)

  ## read the file
  los.data <- read.table(file="los.data.txt",sep=";",header=TRUE,row.names=NULL)
  
  ## check the number of columns of the passed data.frame
  if( dim(los.data)[2] < 7 )
  {
    stop("The passed data.frame doesn't include 7 columns for holding the observation id, observed times and censoring times.")
  }
  
  ## check the observation id
  if( length(los.data[, 1]) !=  length(unique(los.data[, 1])) ) {
    stop("The values of the observation id in the passed data.frame must be unique")    
  }
  
  ## cut the irrelevant columns, which doesn't hold informations for the computation
  los.data <- los.data[,c(1, c(2,3,4,5,6,7))]

  ## set my own columnnames
  names(los.data) <- c("id", "j.01", "j.02", "j.03", "j.12", "j.13", "cens")
  
  ## compute variables cens.0 for admissions censored in the initial state 0
  ## and cens.1 for admissions censored in state 1

  los.data$cens.0 <- los.data$cens
  los.data$cens.0[is.finite(los.data[,2])] <- Inf
    
  los.data$cens.1 <- los.data$cens
  los.data$cens.1[is.infinite(los.data[,2])] <- Inf

    
  los.data <- los.data[,c(1,2,3,4,5,6,8,9)]

 ## los.data <- read.los.data("los.data.txt",pos.id=1,pos.columns=c(2,3,4,5,6,7))
