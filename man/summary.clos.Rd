\name{summary.clos}
\alias{summary.clos}
\alias{print.summary.clos}
\title{Summary method for objects of class 'clos'}
\description{
  Summary method for objects of class \code{clos}.
}
\usage{
summary.clos(object, ...)

print.summary.clos(x, ...)
}
\arguments{
  \item{object}{An object of class \code{clos}.}
  \item{x}{An object of class \code{summary.clos}.}
  \item{...}{other arguments}
}
\value{
  A list is returned by the function \code{summary.clos} with the
  following elements:
  \item{cLOS}{change in LOS}
  \item{patients}{total number of observed patients}
  \item{patients.discharge}{number of patients being discharged}
  \item{patients.death}{number of patients who die}
  \item{patients.cens}{number of patients being censored, i. e. for whom
  neither discharge or death was observed}
  \item{cases}{number of patients who experienced the intermediate event
    (IE)}
  \item{cases.discharge}{number of patients who experienced the IE being discharged}
  \item{cases.death}{number of patients who experienced the IE and died}
  \item{cases.cens}{number of patients who experienced the IE and were
    censored}  
}
\author{Matthias Wangler \email{mw@imbi.uni-freiburg.de}}
\seealso{  \code{\link{clos}}}
\examples{
data(los.data)
my.observ <- prepare.los.data(x=los.data)
trans <- matrix(FALSE,4,4)
diag(trans) <- TRUE
trans[1,] <- TRUE
trans[2,3:4] <- TRUE
my.model <- msmodel(c("0","1","2","3"),trans,cens.name="cens")
los <- clos(model=my.model,observ=my.observ)
summary(los)
}

\keyword{print}
