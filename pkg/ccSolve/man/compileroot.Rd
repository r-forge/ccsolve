\name{compile.multiroot}
\alias{compile.multiroot}

\title{Creates compiled code to solve for the root of multiple nonlinear
equations}
\description{
  Generates compiled code from the body of fortran subroutines 
  or C functions that specify multiple (nonlinear) functions whose roots have to 
  be solved for. These codes can then be used in the functions \code{multiroot} 
  and \code{multiroot.1D} from the R-package \code{rootSolve}. 
}
\usage{
 compile.multiroot (func, jacfunc = NULL, parms = NULL, x = NULL,
      declaration = character(), includes = character(), language = "F95", ...)
}
\arguments{
  \item{func }{A character vector with F95, Fortran, or C code, without 
    declarations, that specifies the functions whose root have to be found. 
    The function is called with arguments \code{n, t, x, f, rpar, ipar} and the 
    The user gets as input a vector \code{x} with the x-values, 
    and must specify the function value whose root is sought, in a vector \code{f}.
    The variables \code{n}, \code{t}, and \code{x} should not be changed in the code.
    \code{rpar} and \code{ipar} are a double and integer vector, passed upon runtime to the solver.     
  }
  \item{jacfunc }{A character vector with F95, Fortran, or C code, without 
    declarations, that specifies the jacobian of the functions whose root have to be found.
    The function is called with arguments \code{n, t, y, df, rpar, ipar} and the 
    jacobian should be written to matrix (language = Fortran, F95) or 
    vector (language = C) \code{df}.
    The variables \code{n}, \code{t}, and \code{x} should not be changed in the code.
  }
  \item{parms }{The parameter vector or list, or a vector with parameter names, 
    declaring the parameters to be used in the differential equations. 
    Only the parameter names will be used and these names declared in the compiled code. 
    The parameter values will be set at run time.
    Parameter values should not be changed in the code.
    The parameter vector or list should be compatible with the \code{parms} vector used in the call to the solver.
  }
  \item{x }{A vector that contains the names of the variables used in the equations. 
    These names will be known in the compiled code. 
    This vector should be compatible with the \code{start} vector used in the call to the solver 
    (see e.g. \link{multiroot}).
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

# F95
model.f <- "
  f(1) = x(1)**2. + x(2)**2. - 1
  f(2) = x(1)**2. - x(2)**2. + 0.5
"
cmodel <- compile.multiroot(model.f)
(SS <- multiroot(f = cmodel, start = c(1,1)))

# retrieve function value
ccfunc(cmodel, SS$root)

# same with x and parameter 'names' and in C

model.C <- "
  f[0] = x1*x1 + x2*x2 - p1;
  f[1] = x1*x1 - x2*x2 + p2;
"
parms <- c(p1 = 1, p2 = 0.5)
cmodel.C <- compile.multiroot(model.C, language = "C", 
  x = c("x1", "x2"), parms = parms)
(SS <- multiroot(f = cmodel.C, start = c(1,1), parms = parms))

## =======================================================================
## example 2
## 3 equations, two solutions
## =======================================================================

model <- "
f(1) = x(1) + x(2) + x(3)**2. - 12.
f(2) = x(1)**2 - x(2) + x(3) - 2.
f(3) = 2.d0 * x(1) - x(2)**2 + x(3) - 1.d0
"
cmod <- compile.multiroot(model)

# two solutions
(ss <- multiroot(f = cmod, start = c(1, 1, 1)))
(ss <- multiroot(cmod, c(0, 0, 0)))

## =======================================================================
## example 3: same, with 3 parameters
## =======================================================================

model2 <- "
f(1) = x(1) + x(2) + x(3)**2. - a
f(2) = x(1)**2 - x(2) + x(3) - b
f(3) = 2.d0 * x(1) - x(2)**2 + x(3) - c
"
parms <- c(a = 12, b = 2, c = 1)
cmod2 <- compile.multiroot(model2, parms = parms)

# two solutions
(ss <- multiroot(f = cmod2, start = c(1, 1, 1), parms = parms))
(ss <- multiroot(cmod2, c(10, 10, 10), parms = parms*20))

}
}
\author{
  Karline Soetaert <karline.soetaert@nioz.nl>
}
\references{
Soetaert K. (2014). rootSolve: Nonlinear root finding, equilibrium and 
steady-state analysis of ordinary differential equations.  R-package
  version 1.7
}
\keyword{ utilities }