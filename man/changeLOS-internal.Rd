\name{changeLOS-internal}
\title{Internal changeLOS objects}
\alias{.First.lib}
\alias{.Last.lib}
\description{Internal changeLOS objects (functions) '.First.lib'and '.Last.lib'
  The package 'changeLOS' contains C/C++ code in its src subdirectory.
  'R CMD INSTALL' creates a shared object. This shared object is loaded
  via a call to 'library.dynam' in the '.First.lib' function.
  '.Last.lib.' is for unloading the shared object.  
}
\details{These are not to be called by the user.}
\keyword{internal}
