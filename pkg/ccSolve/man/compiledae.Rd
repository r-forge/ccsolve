\name{compile.dae}
\alias{compile.dae}

\title{Creates compiled code for solving differential algebraic equations written in implicit form }
\description{
  Generates compiled code from the body of 
  fortran subroutines or C functions that specify an initial value problem of
  differential algebraic equations (DAE) written in implicit form.
  These problems can then be solved with function \link{daspk} from the R-package \code{deSolve} or \code{mebdfi} from the R-package \code{deTestSet}.
}
\usage{
 compile.dae (res, rootfunc = NULL, eventfunc = NULL,       
      parms = NULL, y = NULL, forcings = NULL, outnames = NULL, 
      declaration = character(), includes = character(), 
      language = "F95", ...)
}
\arguments{
  \item{res }{A character vector with F95, Fortran, or C code, without 
    declarations, that specifies the residual function. 
    The function is called with arguments \code{t, y, dy, cj, r, ires, rpar, ipar} 
    and the calculated residuals should be written to vector \code{r}.
    The variables n, t, y, dy and cj should not be changed in the code.
    \code{rpar} and \code{ipar} are a double and integer vector, passed upon runtime
    to the solver. See vignette \code{"compiledCode"} from R-package \code{deSolve}.
  }
  \item{rootfunc }{A character vector with F95, Fortran, C or C++ code, without 
    declarations, that specifies the root function. Only works with \code{daspk}.
    The function is called with arguments \code{n, t, y, nroot, root, rpar, ipar} and the 
    constraint function whose root is sought should be written to vector \code{root}.
    The variables n, t, y and nroot should not be changed in the code.
  }
  \item{eventfunc }{A character vector with F95, Fortran, C or C++ code, without 
    declarations, that specifies the event function. Only works with \code{daspk}.
    The function is called with arguments \code{n, t, y} and the code needs to 
    redefine (elements) of vector \code{y}.
    The variables n, and t should not be changed in the code.
  }
  \item{parms }{The parameter vector or list, or a vector with parameter names, 
    declaring the parameters to be used in the differential equations. 
    Only the parameter names will be used and these names declared in the compiled code. 
    The parameter values will be set at run time.
    Parameter values should not be changed in the code.
    The parameter vector or list should be compatible with the \code{parms} vector used in the call to the solver
    (see e.g. \link{daspk}).
  }
  \item{y }{The initial condition vector or a vector that contains the names of 
    the state variables used in the equations. 
    These names will be known in the compiled code. 
    In addition the derivatives will be declared in the function, as follows: for a state varialbe called sv, 
    then its corresponding derivative will be called dsv. 
    The user must specify these derivatives in the code.
    This vector should be compatible with the \code{y} vector used in the call to the solver 
    (see e.g. \link{daspk}).
  }   
  \item{forcings }{A vector with the names of the forcings to be used in the equations. 
    These names will be declared in the compiled code, 
    and their values updated for each \code{t} when solving the model.
    This vector should be compatible with the \code{forcings} vector used in the call to the solver 
    (see e.g. \link{daspk}).
  }
  \item{outnames }{The names of the output variables used in the equations. 
    These names will be known in the compiled code, and the user needs to specify their values.
    These values will then be used to update the ouput.
    This vector should be compatible with the \code{outnames} vector used in the call to the solver
    (see e.g. \link{daspk}).
  }
  \item{declaration }{Text that enters the declaration section in each function.
  }
  \item{includes }{Code that comes before the functions.
  }
  \item{language }{A character vector that specifies the source code; one of c("F95", "Fortran", "C"); defaults to "F95".
  }
  \item{...}{optional arguments to the generic function (not used).
  }
}
\value{
  An object of class \code{CFunc} or \code{CFuncList}, as defined in the package \code{inline}. 
}  
\details{

The compiled functions are compatible with the C-functions as used in the solvers.

They are defined as:

 \code{res(t, y, dy, cj, r, ires, rpar, ipar) }

 \code{rootfunc(n, t, y, nroot, root, rpar, ipar) }

 \code{eventfunc(n, t, y)}

In case Fortran or F95 is used, \code{n} is an integer value, \code{t} is a double value (typically time), 

\code{ipar} and \code{rpar} are an integer vector, and double vector, whose length 
is defined by arguments \code{ipar} and \code{rpar} when calling the solver).  

\code{y}, and \code{dy} are double vectors of length \code{n}.
\code{r} is a double precision vector of length \code{n} that contains the residual function.
\code{root} is a double precision vector of length \code{n}.

In case C or C++ are used all are pointers. 

\code{ires} and \code{cj} are single values (integer and double), that are not used.

The variables n, t, and -except for eventfunc- y, should not be changed within the function; they are input only.

The user needs to specify \code{r, root}, and \code{y} for \code{res, rootfunc} 
and \code{eventfunc} respectively.
}

\note{

Specifying the jacobian of the residual function cannot be done in compiled code (this is because the number of equations is not known) and this value is not passed to the subroutines. 
See however \code{vignette("compiledCode")} how to write this using regular text files.

}
\examples{
\dontrun{
## =======================================================================
## ex. 1
## A simple DAE
## y1' = y2
## y1 = cos(x)
## (analytical solution = y1 = cos(x); y2 = -sin(x)
## =======================================================================

# the R-version
resDAE <- function(t, y, dy, parms) {

    f1     <- dy[1] - y[2]
    f2     <- y[1] - cos(t)

    list(c(f1, f2))
}

yini  <- c(y1 = cos(0), y2 = -sin(0))
dyini <- c(-sin(0), -cos(0))
times <- seq(from = 0, to = 10, by = 0.1)
index <- c(1, 1, 0)
print(system.time(out   <- daspk(times = times, res = resDAE, y = yini,
  dy = dyini, parms = NULL, nind = index)))

# The same, in F95
fresDAE <- " 
   r(1) = dy(1) - y(2)
   r(2) = y(1) - cos(t)
"
cDAE <- compile.dae(res = fresDAE) 
print(system.time(out2 <- daspk(times = times, res = cDAE, y = yini,
  dy = dyini, parms = NULL, nind = index)))

plot(out, out2)

# function value
ccfunc(cDAE, y = yini, dy = dyini, times = 1, parms = NULL)

## =======================================================================
##
## The chemical model example of daspk, implemented as a DLL
##
## =======================================================================

# Dissociation constant
K <- 1

# parameters
pars <- c(K     = K  ,
          ka    = 1e6,     # forward rate
          rr    = 1  ,
          prod  = 0.1)

#---------------------------------------------------------
# problem formulation as R-function
#---------------------------------------------------------

Chemres <- function (t, y, dy, pars){
  with (as.list(c(y, dy, pars)), {
    ra  <- ka * D        # forward rate
    rb  <- ka/K * A * B  # backward rate
    # residuals of rates of changes
    res1 <- -dD  - ra + rb + prod
    res2 <- -dA  + ra - rb
    res3 <- -dB  + ra - rb - rr*B
    return(list(res = c(res1, res2, res3),
                CONC = A + B + D))
  })
}
times <- seq(0, 100, by = 2)

# Initial conc and rate of change; D is in equilibrium with A,B
y  <- c(A = 2, B = 3, D = 2*3/K)
dy <- c(dA = 0, dB = 0, dD = 0)

PD <- matrix(nr = 3, nc = 3, 0)

print(system.time(
  outR <- daspk(y = y, dy = dy, times = times, res = Chemres,
                parms = pars, atol = 1e-10, rtol = 1e-10)
))

# -------------------
# F95
# -------------------
DAE.f <- "
      ra  = ka* y(3)            ! forward rate
      rb  = ka/K *y(1) * y(2)   ! backward rate

      ! residuals of rates of changes
      r(3) = -dy(3)  -  ra + rb + prod
      r(1) = -dy(1)  +  ra - rb
      r(2) = -dy(2)  +  ra - rb - rr*y(2)
      CONC   = y(1) + y(2) + y(3)
"      
cDAE <- compile.dae(res = DAE.f,  parms = pars, outnames = "CONC",
  declaration = "DOUBLE PRECISION :: ra, rb")

print(system.time(
out <- daspk(y = y, dy = dy, times = times, res = cDAE, outnames = "CONC",
               nout = 1, parms = pars, atol = 1e-10, rtol = 1e-10)
))

## =======================================================================
##  
## The pendulum equation - with parameter L
## x' = u
## y' = v
## u' = -lam*x
## v' = -lam*y - g
## 0 = x^2 + y^2 - L^2
## =======================================================================

# implementation in R - using daspk (!)
pendulum <- function (t, y, dy, L) {
 list(c(-dy[1] + y[3]          ,
        -dy[2] + y[4]          ,
        -dy[3] -y[5]*y[1]      ,
        -dy[4] -y[5]*y[2] - 9.8,
         y[1]^2 + y[2]^2 - L
     ))
}

yini  <- c(x = 1,  y = 0, u = 0,  v = 1   , lam = 1)
dyini <- c(dx = 0,dy = 1,du = -1,dv = -9.8,dlam = 3*9.8)
times <- seq(from = 0, to = 10, by = 0.01)

parms <- c(L = 1.)

index3 <- c(2, 2, 1)
out <- daspk (y = yini, dy = dyini, res = pendulum, 
              times = times, parms = parms, nind = index3)

# in fortran
fpendulum <- "
  r(1) = -dy(1) + y(3)
  r(2) = -dy(2) + y(4)
  r(3) = -dy(3) -y(5)*y(1)
  r(4) = -dy(4) -y(5)*y(2) - 9.8d0
  r(5) = y(1)**2. + y(2)**2. - L
"
cpendulum <- compile.dae (res = fpendulum, parms = parms)
outb <- daspk (y = yini, dy = dyini, res = cpendulum, 
              times = times, parms = parms, nind = index3)

plot(out, outb, which = 1:4)
}
}                                                                                          
\author{
  Karline Soetaert <karline.soetaert@nioz.nl>
}
\references{
Karline Soetaert, Thomas Petzoldt, R. Woodrow Setzer (2010).  
  Solving Differential Equations in R: Package deSolve. Journal of Statistical
  Software, 33(9), 1--25.  URL http://www.jstatsoft.org/v33/i09/.

Karline Soetaert, Thomas Petzoldt, R. Woodrow Setzer.
  R Package deSolve, Writing Code in Compiled Languages.
  package vignette.
  
Karline Soetaert, Jeff Cash and Francesca Mazzia (2014). 
deTestSet: Testset for differential equations. R package version 1.1.1.
  
}
\keyword{ utilities }
