"read.los.data" <-
function( file, sep=";", header=TRUE, row.names=NULL, pos.id=1, pos.columns) {
## --------------------------------------------------------------------------------
## Title: R-function read.los.data()
## ---------------------------------------------------------------------------------
## Author: Matthias Wangler
##         mw@imbi.uni-freiburg.de
##         Institute of Med. Biometry and Med. Computer Science
##         Stefan-Meier-Strasse 26, D-79104 Freiburg,
##         http://www.imbi.uni-freiburg.de
## ---------------------------------------------------------------------------------
## Description: Read and prepare a data set which can be passed to the function clos
## ---------------------------------------------------------------------------------
## Required Packages: -
## ---------------------------------------------------------------------------------
## Usage: read.los.data(file ,sep=";",header=TRUE,row.names=1,pos.columns)
##
## file:       the name of the file which the data are to be read from
## sep:        the field separator character
## header:     a logical value indicating whether the file contains the
##             names of the variables as its first line
## row.names  a vector of row names.  This can be a vector giving the
##            actual row names, or a single number giving the column of the
##            table which contains the row names, or character string
##            giving the name of the table column containing the row names.
##
##            If there is a header and the first row contains one fewer
##            field than the number of columns, the first column in the
##            input is used for the row names.  Otherwise if 'row.names' is
##            missing, the rows are numbered.
##
##            Using 'row.names = NULL' forces row numbering.
## pos.id     the position of the unique observation id (patient id, admision id)
## pos.columns the positions of the columns which are holding the observed times:
##             pos.columns[1]: jump form 0 to 1
##             pos.columns[2]: jump form 0 to 2
##             pos.columns[3]: jump form 0 to 3
##             pos.columns[4]: jump form 1 to 2
##             pos.columns[5]: jump form 1 to 3
##             pos.columns[6]: cens
## ---------------------------------------------------------------------------------
## Value: data.frame of the form data.frame( id, j.01, j.02, j.03, j.12, j.13, cens):
##
## id.      id (patient id, admision id, ...)
## j.01:    observed time for jump from 0 to 1
## j.02:    observed time for jump from 0 to 2
## j.03:    observed time for jump from 0 to 3
## j.12:    observed time for jump from 1 to 2
## j.13:    observed time for jump from 1 to 3
## cens:    observed time for censoring
## ---------------------------------------------------------------------------------
## Notes: -
## ---------------------------------------------------------------------------------
## Example: > los.data <- read.los.data(file="los.data.txt",pos.columns=c(2,3,4,5,6,7))
## ---------------------------------------------------------------------------------
## License: GPL 2
##----------------------------------------------------------------------------------
## History:    20.06.2004, Matthias Wangler
##                         first version
## ---------------------------------------------------------------------------------
  ## check the positions of the columns holding the observation id
  if( !(is.numeric(pos.id) && length( pos.id ) == 1) )
  {
    stop("The passed argument 'id' must hold one position")    
  }
  
  pos.columns <- unique( pos.columns )
  ## check the positions of the columns holding the observed times and censoring times
  if( length( pos.columns ) != 6 )
  {
    stop("The passed positions must hold 6 different positions")    
  }
  ## read the file
  los.data <- read.table(file=file,sep=sep,header=header,row.names=row.names)
  
  ## check the number of columns of the passed data.frame
  if( dim(los.data)[2] < 7 )
  {
    stop("The passed data.frame doesn't include 7 columns for holding the observation id, observed times and censoring times.")
  }
  
  ## check the observation id
  if( length(los.data[, pos.id]) !=  length(unique(los.data[, pos.id])) ) {
    stop("The values of the observation id in the passed data.frame must be unique")    
  }
  
  ## cut the irrelevant columns, which doesn't hold informations for the computation
  los.data <- los.data[,c(pos.id,pos.columns)]

  ## set my own columnnames
  names(los.data) <- c("id", "j.01", "j.02", "j.03", "j.12", "j.13", "cens")
  
  return(los.data)
}

