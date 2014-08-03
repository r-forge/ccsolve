\name{print.code}
\alias{print.code}
\alias{save.cc}
\alias{load.cc}
\title{
  prints generated codes and code snippets; saves and loads compiled codes
}

\usage{
print.code(x)
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
# print a code snippit
print.code(tran1d)

# generate a code and print it
String <- " f(1) = 3*x(1) + x(2)*x(3)"
rfun <- compile.root(String)
print.code(rfun)

\dontrun{

# save the compiled code for later use - will create two files
save.cc(rfun, file = "test.cc")

# load it, call it rr
rr <- load.cc("test.cc")
multiroot(rr, 1:3, pos = TRUE)
}
}

\references{
}

\seealso{
}


