\name{read.los.data}
\alias{read.los.data}
\title{Read the data for clos}
\description{Read data fom a textfile, which contains one row per
  individual and one time variable per
  possible transition.}
\usage{read.los.data(file, sep = ";", header = TRUE, row.names=NULL,
pos.id=1, pos.columns) }
\arguments{
  \item{file}{name of the file to be read. If it does not contain an
           absolute path, the file name is relative to the current working directory, \code{getwd()}. }
\item{sep}{the field separator character}
\item{header}{a logical value indicating whether the file contains the names of the variables as its first line}
\item{row.names}{a vector of row names.  This can be a vector giving the
          actual row names, or a single number giving the column of the
          table which contains the row names, or character string
          giving the name of the table column containing the row names.

          If there is a header and the first row contains one fewer
          field than the number of columns, the first column in the
          input is used for the row names.  Otherwise if 'row.names' is
          missing, the rows are numbered.

          Using 'row.names = NULL' forces row numbering.}
\item{pos.id}{the position of the unique id (patient id, admision id)}	
\item{pos.columns}{the positions of the columns which are holding the observed times:
   \describe{ 
     \item{pos.columns[1]:}{ transition from initial state to intermediate state} 
     \item{pos.columns[2]:}{ transition from initial state to absorbing
       state (discharge)} 
     \item{pos.columns[3]:}{ transition from initial state to competing
       absorbing state (death)}
     \item{pos.columns[4]:}{ transition from intermediate state to absorbing state}
     \item{pos.columns[5]:}{ transition from intermediate state to competing absorbing state}
     \item{pos.columns[6]:}{ censoring time (either in initial or intermediate state)}
    }
 } 
}
\details{
  The data textfile to be read must contain one time variable per
  possible transition and one row per individual. An additional variable
  contains censoring times. If a transition was \emph{not} observed for
  an individual, the respective row of the textfile has entry \code{Inf}
  (`infinite') for that variable.

  This data structure is well suited for bootstrapping with respect to
  the number of individuals. It is still concise for the
  four-state-model for change in LOS; note that this model does not
  allow for backward transitions. However, for a larger number of possible
  transitions, this data structure will be less desirable.
}
\value{ 
  a data frame of the form data.frame(id,j.01,j.02,j.03,j.12,j.13,cens.0,cens.1)
  \item{id}{ id (patient id, admision id)}
  \item{j.01}{ observed time for jump from 0 (initial state) to 1
    (intermediate state)}
  \item{j.02}{ observed time for jump from 0 to 2 (discharge)}
  \item{j.03}{ observed time for jump from 0 to 3 (death)}
  \item{j.12}{ observed time for jump from 1 to 2}
  \item{j.13}{ observed time for jump from 1 to 3}
  \item{cens.0}{ observed time for censored in 0}
  \item{cens.1}{ observed time for censoerd in 1}
}
\author{ Matthias Wangler \email{mw@imbi.uni-freiburg.de}}
\seealso{ \code{\link{clos}},  \code{\link{prepare.los.data}} }
\examples{
## los.data <- read.los.data("los.data.txt",pos.id=1,pos.columns=c(2,3,4,5,6,7))
## Results in the same data frame as:
data(los.data)
}

\keyword{file}
\keyword{datasets}
