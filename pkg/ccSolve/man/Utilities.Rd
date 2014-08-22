\name{printCode}
\alias{printCode}
\alias{save.cc}
\alias{load.cc}
\title{
  Utilities for ccSolve.
}
\description{
  prints generated codes and code snippets; saves and loads compiled codes
}

\usage{
printCode(x)
save.cc(x, file)
load.cc(file)
}
\arguments{
  \item{x }{Code (snippit) or compiled object
  }
  \item{file }{The name of the file where the compiled code object will be saved
   or loaded from. For use with \code{load.cc}, the extension of the file name 
   should be `.cc'.
  } 

}
\details{
  As the compilation of a problem is the time-limiting factor, a compiled object
  can be saved with \code{save.cc} and then later loaded with \code{load.cc}.

  For \code{save.cc} two files will be saved, one with the dll or shared object 
  (extension ".dll" or ".so"), and one with the compiled code object 
  (without pointer), with extension ".cc".
  
  During \code{load.cc} the shared object or dll will be loaded and address pointers
  to each function set.
}
\author{
  Karline Soetaert <karline.soetaert@nioz.nl>
}
\examples{
# generate a code and print it
String <- " 
f(1) = 3*x(1) + x(2)*x(3) +2
f(2) = 2*x(2) + 3
f(3) = x(3)*(x(1)+1)"
rfun <- compile.multiroot(String)
printCode(rfun)
multiroot(rfun, 1:3)

\dontrun{

# save the compiled code for later use - will create two files
save.cc(rfun, file = "test.cc")

# load it, call it rr
rr <- load.cc("test.cc")
multiroot(rr, 1:3)
}
}

\references{
A lot of the code used to create the compiled codes is modified from the R-package
inline.
The modifications relate to the part of function \code{cfunction} 
to create inline code for fortran and fortran 95. 

Oleg Sklyar, Duncan Murdoch, Mike Smith, Dirk Eddelbuettel and Romain Francois (2013). 
inline: Inline C, C++, Fortran function calls from R. R package version
  0.3.13. http://CRAN.R-project.org/package=inline 
}



