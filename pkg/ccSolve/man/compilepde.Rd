\name{compile.pde}
\alias{compile.pde}
\title{
  compiles partial differential equation problems
}
\description{
  Uses inline compilation to generate compiled code from the body of 
  fortran or C function that specify an partial differential equation problem 
  to be solved with functions \code{ode.1D}, \code{ode.2D} or \code{ode.3D}
  from the R-package \pkg{deSolve} or  \code{steady.1D}, \code{steady.2D} or 
  \code{steady.3D} from \pkg{rootSolve}.
}

\usage{
 compile.pde (func, dim, names, 
      parms = NULL, forcings = NULL, outnames = NULL, 
      declaration = character(), includes = character(), language = "F95", ...)
}
\arguments{
  \item{func }{A character vector with F95, Fortran, C or C++ code, without 
    declarations, that specifies the derivative function. 
    The function is called with arguments \code{n, t, y, f, rpar, ipar} and the 
    derivative should be written to vector \code{f}.
    The variables n (integer), t, and y (double) should not be changed in the code.
  }
  \item{dim }{The dimension of the problem, one value (1-D), or a 2-valued (2.D) 
  or 3-valued (3.D) vector with the number of boxes in x, y, and z direction.
    This vector should be compatible with the \code{dim} argument used in the call to the solver.
  }  
  \item{names }{The names of the state variables used in the equations. 
    These names will be known in the functions; 
    also the derivatives will be declared in the function; if a state varialbe is called sv, 
    then its corresponding derivative will be called dsv and the user must specify the derivative in the code.
    This vector should be compatible with the \code{names} argument used in the call to the solver.
  }
  \item{parms }{The parameter vector or list to be used in the equations. 
    Their names will be known in the functions, and their values passed 
    upon solving the model. Parameter values should not be changed in the compiled code.
    This vector or list should thus be compatible with the \code{parms} vector used in the call to the solver.
  }
  \item{forcings }{The names of the forcings to be used in the equations. 
    These names will be known in the functions, and their values updated for each \code{t} when solving the model.
    This vector should thus be compatible with the \code{forcingss} vector used in the call to the solver.
  }
  \item{outnames }{The names of the output variables used in the equations. 
    These names will be known in the functions, and their values used to update the ouput.
    This vector should thus be compatible with the \code{outnames} vector used in the call to the solver.
  }
  \item{declaration }{Text that enters the declaration section in each function.
  }
  \item{includes }{Code that comes before the functions.
  }
  \item{language }{A character vector that specifies the source code; one of c("F95", "Fortran", "C++", "C") defaults to "F95".
  }
  \item{...}{optional arguments to the generic function (not used).
  }
}
\value{
  An object of class \code{CFunc} or \code{CFuncList}, as defined in the package \code{inline}. 
  
  Each object of class \code{CFunc} can be called with the appropriate arguments
}  
\details{

The compiled function is defined as:

 \code{func(n, t, y, f, rpar, ipar) }

In case Fortran or F95 is used, \code{n} are an integer value, \code{t} is a double value (typically time), 

\code{ipar} and \code{rpar} are an integer vector, and double vector, whose length 
is defined by arguments \code{ipar} and \code{rpar} when calling the solver).  

\code{y}, and \code{f} are double vectors of length \code{n}.
\code{df} is a double precision matrix of dimension \code{n, n} that contains the jacobian.
\code{root} is a double precision vector of length \code{n}.

In case C or C++ are used all are pointers. 

the values of \code{df} are set to 0 at the beginning of the subroutine.

The variables n, t, and -except for eventfunc- y, should not be changed within the function; they are input only.
The user needs to specify \code{f, df, root}, and \code{y} for \code{func, jacfunc, rootfunc} 
and \code{eventfunc} respectively.

}

\author{
  Karline Soetaert <karline.soetaert@nioz.nl>
}
\examples{
\dontrun{
printCode(tran1d)

## ===========================================================================
## a simple 1D transport model
## ===========================================================================

declaration <- "   double precision :: Fin, Fout"

transport <- "

    ! tran1d    (nx,C, Cup, Cdown, Fin, Fout,BcUp,BcDown,K,v,dx,dC) 
    CALL tran1d (nx,C, Cup, Cdown, Fin, Fout, 2, 2, Kx, v, dx, dC) 

    ! reaction
    dC = dC - k*C

"

parms <- c(k = 0.05, Cup = 20, Cdown = 30, Kx = 1, v = 1, dx = 1)
ctran <- compile.pde(func = transport, parms = parms, 
   includes = tran1d, declaration = declaration, 
   names = "C", dim = 100)

out <- ode.1D(y = rep(0,100), parms = parms, func = ctran, times = 0:365,
   names = "C", dim = 100)
std <- steady.1D(y = rep(0,100), parms = parms, func = ctran,
    names = "C", dim = 100)
image(out)

## ===========================================================================
## a simple 2D transport model
## ===========================================================================

declar2D <- "   double precision :: Fxup(ny), Fxdown(ny), Fyup(nx), Fydown(nx)"

trans2D <- "
    ! tran2d (nx, ny, C, C_xup, C_xdown, C_yup, C_ydown, 
    !    F_xup, F_xdown, F_yup, F_ydown, BcxUp,BcxDown,BcyUp, Bcydown,
    !    Kx,Ky,vx,vy,dx,dy,dC) 

    CALL tran2d (nx, ny, C, Cxup, Cxdown, Cyup, Cydown,       &
      Fxup, Fxdown, Fyup, Fydown, 2, 2, 2, 2,                 &
      Kx, Ky, vx, vy, dx, dy, dC) 

    ! reaction
    dC = dC - k*C

"

parms <- c(k = 0.05, Cxup = 20, Cxdown = 30, Cyup = 20, Cydown = 30,
  Kx = 1, Ky = 1, vx = 1, vy = 1, dx = 1, dy = 1)
nx <- 30
ny <- 30
ctran2D <- compile.pde(func = trans2D, parms = parms, 
   includes = tran2d, declaration = declar2D, 
   names = "C", dim = c(nx, ny))

out <- ode.2D(y = rep(0,nx*ny), parms = parms, func = ctran2D, times = 0:365,
   names = "C", dim = c(nx,ny), lrw = 50000)
std <- steady.2D(y = rep(0,nx*ny), parms = parms, func = ctran2D,
    names = "C", dim = c(nx,ny), lrw = 50000)
image(std)

## ===========================================================================
## a two-species transport model in 1-D
## ===========================================================================

declaration = " 
   DOUBLE PRECISION :: Ingestion(nx), GrowthPrey(nx), MortPredator(nx)
   DOUBLE PRECISION :: Fin, Fout"

LV1D <- "

 ! transport
   CALL tran1d (nx, prey, preyup, preydown, Fin, Fout, 2, 2,      &
                 D, v, dx, dprey) 
   
   CALL tran1d (nx, Predator, predup, preddown, Fin, Fout, 2, 2,  &
                 D, v, dx, dPredator) 

 ! biology
   Ingestion    = rIng  * Prey * Predator
   GrowthPrey   = rGrow * Prey * (1 - Prey/Kprey)
   MortPredator = rMort * Predator
 
   dPrey        = dprey + GrowthPrey - Ingestion
   dPredator    = dpredator + Ingestion * assEff - MortPredator

"

nb <- 2000
parms <- c(preyup = 2, preydown = 2, predup = 1, preddown = 1,
           D = 10, v = 10, dx = 0.5, 
           rIng = 0.2, rGrow = 1.0, rMort = 0.2, assEff = 0.5, Kprey = 10)

cpredprey1D <- compile.pde(func = LV1D, parms = parms, includes = tran1d, 
  names = c("prey", "predator"), dim = nb, declaration = declaration)

yini  <- c(rep(1, nb), rep(2, nb))
times <- seq(0, 200, by = 1)

print(system.time(
out   <- ode.1D(y = yini, times = times, func = cpredprey1D,
  names = c("prey", "predator"), dim = nb,  
  parms = parms, nspec = 2, method = lsodes)
))

image(out)

}
}
\references{
}

\seealso{
}

\keyword{datasets}

