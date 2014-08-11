\name{compile.root}
\alias{compile.root}

\title{Creates compiled code for the rootSolve function multiroot and multiroot.1D }
\description{
  Uses inline compilation to generate compiled code from the body of 
  fortran or C function that specify aroot solving function.
}
\usage{
 compile.root (func, jacfunc = NULL, parms = NULL, 
      declaration = character(), includes = character(), language = "F95", ...)
}
\arguments{
  \item{func }{A character vector with F95, Fortran, C or C++ code, without 
    declarations, that specifies the root function. 
    The user gets as input a vector \code{x} with the x-values, 
    and must specify the function value whose root is sought, in a vector \code{f}.
  }
  \item{jacfunc }{A character vector with F95, Fortran, C or C++ code, without 
    declarations, that specifies the jacobian of the function whose root has to be found.
    The user gets as input a vector \code{x} with the x-values, 
    and must specify the jacobian as a matrix (language = Fortran, F95) or 
    vector (language = C) \code{df}.
  }
  \item{parms }{The parameter vector or list to be used in the equations. 
    Their names will be known in the functions, and their values passed 
    upon solving the model. Parameter values should not be changed in the compiled code.
    This vector or list should thus be compatible with the \code{parms} vector used in the call to the solver.
  }
  \item{declaration }{Text that eneters the declaration section in each function.
  }
  \item{includes }{Code that comes before the functions.
  }
  \item{language }{A character vector that specifies the source code; one of c("F95", "Fortran", "C++", "C") defaults to "F95".
  }
  \item{...}{optional arguments to the generic function (not used).
  }
}
\value{
  An object of class \code{CFunc} or \code{CFuncList}, as from the package  \code{inline}. 
  
  Each object of class \code{CFunc} can be called with the appropriate arguments
}  
\details{

The compiled functions are compatible with the C-functions as used in the solvers.

They are defined as:

 \code{func(n, t, x, f, rpar, ipar) }

 \code{jacfunc(n, t, x, df, rpar, ipar) }
 
In case Fortran or F95 is used, \code{n} is an integer value, 
\code{t} is a double value (not used),

\code{ipar} and \code{rpar} are an integer vector, and double vector that are not used).  

\code{x}, and \code{f} are double vectors of length \code{n}.
\code{df} is a double precision matrix of dimension \code{n, n} that contains the jacobian.

In case C or C++ are used all are pointers. 

the values of \code{df} are set to 0 at the beginning of the subroutine.

The variables n, t, ipar and rpar should be ignored.
The user needs to specify \code{f, df}, based on \code{x}

}

\note{


}
\examples{
\dontrun{

## =======================================================================
## example 1  
## 2 simultaneous equations
## =======================================================================

# R-function
model <- function(x) c(F1 = x[1]^2+ x[2]^2 -1, 
                       F2 = x[1]^2- x[2]^2 +0.5)

(ss <- multiroot(f = model, start = c(1, 1)))

# fortran
model.f <- "
  f(1) = x(1)**2. + x(2)**2. - 1
  f(2) = x(1)**2. - x(2)**2. + 0.5
"
cmodel <- compile.root(model.f)
(SS <- multiroot(f = cmodel, start = c(1,1)))

## =======================================================================
## example 2
## 3 equations, two solutions
## =======================================================================

model <- "
f(1) = x(1) + x(2) + x(3)**2. - 12.
f(2) = x(1)**2 - x(2) + x(3) - 2.
f(3) = 2.d0 * x(1) - x(2)**2 + x(3) - 1.d0
"
cmod <- compile.root(model)

# two solutions
(ss <- multiroot(f = cmod, start = c(1, 1, 1)))
(ss <- multiroot(cmod, c(0, 0, 0)))
####  cmod(ss$root) karline: should be fixed so that this can be queried

## =======================================================================
## example 3: same, with 3 parameters
## =======================================================================

model2 <- "
f(1) = x(1) + x(2) + x(3)**2. - a
f(2) = x(1)**2 - x(2) + x(3) - b
f(3) = 2.d0 * x(1) - x(2)**2 + x(3) - c
"
parms <- c(a = 12, b = 2, c = 1)
cmod2 <- compile.root(model2, parms = parms)

# two solutions
(ss <- multiroot(f = cmod2, start = c(1, 1, 1), parms = parms))
(ss <- multiroot(cmod2, c(10, 10, 10), parms = parms*20))
####  cmod(ss$root) karline: should be fixed so that this can be queried
}
}
\author{
  Karline Soetaert <karline.soetaert@nioz.nl>
}

\keyword{ utilities }