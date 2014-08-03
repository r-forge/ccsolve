\name{ccSolve-package}
\alias{ccSolve-package}
\alias{ccSolve}
\docType{package}
\title{
  creating compiled code problems for solver packages deSolve, bvpSolve, rootSolve.
}
\description{
  Functions to solve differential equation problems using inline 
   compiled code, written as the body of Fortran, F95, C (and C++) code.
   
  Can be used with solvers from packages deSolve, bvpSolve, deTestSet, 
  and rootSolve.
}
\details{
  \tabular{ll}{
    Package: \tab ccSolve\cr
    Type:    \tab Package\cr
  }

  The system of differential equations or function whose root is to be found,
  is written as the body of Fortran, F95 or C code. 
  
  These are then provided by the required (technical) interface,
  compiled and the function address stored in a list (class CFunc).   
  
  This list can then be passed to any of the solver routines.
  
  The idea of inline compilation is based on the R-package \code{inline}, maintained by Oleg Sklyvar.
  However, the current version of this package is not flexible enough 
  to generate Fortran inline code in a robust way. 
  Therefore the package \code{inlineSolve} for now uses functions very similar to functions 
  of package \code{inline}, but that extend them with the required functionality.
  
  Hopefully this code duplication will not be necessary in the future.
  
}
\author{
  Karline Soetaert (Maintainer)
}
\references{

}
\seealso{
  \code{\link{ode}} solver from the deSolve package.

  \code{\link{bvptwp}},  \code{\link{bvpcol}}, solvers from the bvpSolve package.
}

\examples{
\dontrun{
## show examples (see respective help pages for details)
example(compile.bvp)

example(compile.ode)

example(compile.dae)

example(compile.pde)

}
}

\keyword{ package }
