\name{compile.bvp}
\alias{compile.bvp}

\title{Creates compiled code for the bvpSolve functions }
\description{
  Uses inline compilation to generate compiled code from the body of 
  fortran or C function that specify a boundary value problem.
}
\usage{
 compile.bvp    (func, jacfunc = NULL, bound = NULL, jacbound = NULL, 
      parms = NULL, forcings = NULL, outnames = NULL, yini = NULL,
      declaration = character(), includes = character(), language = "F95", ...)
}
\arguments{
  \item{func }{A character vector with F95, Fortran, C or C++ code, without 
    declarations, that specifies the derivative function. 
    The function is called with arguments \code{n, x, y, f, rpar, ipar} and the 
    derivative should be written to vector \code{f}. 
    The variables n (integer), x, and y (double) should not be changed in the code.
  }
  \item{jacfunc }{A character vector with F95, Fortran, C or C++ code, without 
    declarations, that specifies the jacobian of the derivative function.
    The function is called with arguments \code{n, x, y, df, rpar, ipar} and the 
    jacobian should be written to the matrix (language = Fortran) or 
    vector (language C) \code{df}.
  }
  \item{bound }{A character vector with F95, Fortran, C or C++ code, without 
    declarations, that specifies the boundary function.
    The function is called with arguments \code{i, n, y, g, rpar, ipar} and the 
    result should be written to value \code{g}.
  }
  \item{jacbound }{A character vector with F95, Fortran, C or C++ code, without 
    declarations, that specifies the jacobian of the boundary function.
    The function is called with arguments \code{i, n, y, dg, rpar, ipar} and the 
    result should be written to vector \code{dg}.
  }
  \item{parms }{The parameter vector or list to be used in the equations. 
    Their names will be known in the functions, and their values passed 
    upon solving the model. Parameter values should not be changed in the compiled code.
    This vector or list should thus be compatible with the \code{parms} vector used in the call to the solver.
  }
  \item{forcings }{The names of the forcings to be used in the equations. 
    These names will be known in the functions, and their values updated for each x when solving the model.
    This vector should thus be compatible with the \code{forcingss} vector used in the call to the solver.
  }
  \item{outnames }{The names of the output variables used in the equations. 
    These names will be known in the functions, and their values used to update the ouput.
    This vector should thus be compatible with the \code{outnames} vector used in the call to the solver.
  }
  \item{yini }{The initial condition vector or a vector that contains the names of the state variables used in the equations. 
    These names will be known in the functions; also the derivatives will be declared in the function.
    If a state varialbe is called sv, 
    then its corresponding derivative will be called dsv; the user must specify the derivative in the code.    
    This vector should be compatible with the \code{yini} vector used in the call to the solver.
  }
  \item{declaration }{Text that enters the declaration section in each function.
  }
  \item{includes }{Code that comes before the functions.
  }
  \item{language }{A character vector that specifies the source code; 
    one of c("F95", "Fortran", "C++", "C") defaults to "F95".
  }
  \item{...}{optional arguments to the generic function (not used).
  }
}
\value{
  An object of class \code{CFunc} or \code{CFuncList}, as from the package  \code{inline}. 
  
  Each object of class \code{CFunc} can be called with the appropriate arguments.
}  
\details{

The  compiled functions are compatible with the C-functions as used in the solvers.

They are defined as:

 \code{func(n, x, y, f, rpar, ipar) }

 \code{jacfunc(n, x, y, df, rpar, ipar) }

 \code{bound(i, n, x, y, g, rpar, ipar) }

 \code{jacbound(i, n, x, y, dg, rpar, ipar) }

In case Fortran or F95 is used, \code{n, i} are an integer value, \code{x} is a double value, 

\code{ipar} and \code{rpar} are an integer vector, and double vector, whose length 
is defined by arguments \code{ipar} and \code{rpar} when calling the solver).  

\code{y}, and \code{f} and \code{dg} are double vectors of length \code{n}.
\code{df} is a double precision matrix of dimension \code{n, n}.

In case C or C++ are used all are pointers. 

The variables n, x, y, should not be changed within the function; they are input only.
The user needs to specify \code{f, df, g}, and \code{dg} respectively.

the values of \code{dg} and \code{df} are set to 0 at the beginning of the subroutine.
}

\note{


}
\examples{
\dontrun{
## =============================================================================
## Standard linear problem with boundary layer at the origin
##
##
## d2y/dt^2=-3py/(p+t^2)^2
## y(t= -0.1)=-0.1/sqrt(p+0.01)
## y(t=  0.1)= 0.1/sqrt(p+0.01)
## where p = 1e-5
##
## analytical solution y(t) = t/sqrt(p + t^2).
##
## =============================================================================

## ------------------------
## Inline Fortran, rewritten as first-order equations
## ------------------------
x <- seq(-1, 1, length.out = 100)
parms <- c(p = 1e-5)

fbnd <- "
  F(1) = y(2)
  F(2) = - 3.d0 * p * y(1)/(p + x*x)**2.
"
fjac <- "
  df(1, 2) = 1.d0
  df(2, 1) = -3.d0*p/(p +x*x)**2.
"
fbound <- "
  if (i == 1) then
     g = (y(1) + 0.1/sqrt(p + 0.01d0))
  else  
     g = (y(1) - 0.1/sqrt(p + 0.01d0))
  end if
"
fjacbound <- "
  dg(1) = 1.d0
"

# compile it - passing parms defines common block
cbnd <- compile.bvp(func = fbnd, jacfunc = fjac, bound = fbound, 
  jacbound = fjacbound, parms = parms)

# show the functions
#print.CFunc(cbnd[[1]])

## ------------------------
## Solution
## ------------------------

print(system.time(Sol <- bvptwp(x = x, leftbc = 1, func = cbnd, 
        n = 2, verbose = TRUE, parms = parms)))
plot(Sol, which = 1)

## ============================================================================
## Same using C and second-order equation
## ============================================================================

## ------------------------
## Inline C++,
## ------------------------
x <- seq(-1, 1, length.out = 100)
parms <- c(p = 1e-5)

fbnd.C <- "
  f[0] = - 3.0 * p * y[0]/((p + x[0]*x[0])*(p + x[0]*x[0]));
"

fbound.C <- "
  if (*i == 1) 
     g[0] = (y[0] + 0.1/sqrt(p + 0.01));
  else  
     g[0] = (y[0] - 0.1/sqrt(p + 0.01));
"


# compile it - passing parms defines common block
cbnd.C <- compile.bvp(func = fbnd.C, bound = fbound.C, parms = parms,
   language = "C")

# show the functions
#print.CFunc(cbnd.C[[1]])

## ------------------------
## Solution
## ------------------------

print(system.time(Sol <- bvpcol(x = x, leftbc = 1, func = cbnd.C, order = 2,
        n = 2, verbose = TRUE, parms = parms)))
plot(Sol, which = 1)


## ==================  tubular reactor with axial dispersion ===================
## y''=Pe(y'+Ry^n) Pe=1,R=2,n=2
## y'(0) = Pe (y(0)-1),
## y'(1)=0
##
## dy=y2
## dy2=Pe(dy-Ry^n)
##
## The initial condition y'(0) is a function of y(0)
## =============================================================================

# 2nd order equation solved and it uses the parameter vector 
freac2 <- "
        F(1) = Pe * (y(2) + R*(Y(1)**2.))
"

greac <- "
        if (i == 1) g = (y(2)-Pe*(y(1)-1.d0))
        if (i == 2) g = (y(2))
"

parms <- c(Pe = 1, R = 2)

cfreac2 <- compile.bvp(func = freac2, bound = greac, parms = parms, 
  language = "Fortran")

## ------------------------
## Solution
## ------------------------

sol <- bvpcol(func = cfreac2, x = seq(0, 1, by = 0.01), order = 2,
              leftbc = 1, ynames = c("y", "dy"), parms = parms)
# second run, with change parameters
sol2 <- bvpcol(func = cfreac2, x = seq(0, 1, by = 0.01), order = 2,
              leftbc = 1, ynames = c("y", "dy"), parms = parms*2)

plot(sol, sol2)

## ==================  tubular reactor with axial dispersion ===================
## but now R is a forcing function, linearly interpolated based on data
## ==================  tubular reactor with axial dispersion ===================

parms <- c(Pe = 1)
forcs <- "R"

cfreac.F <- compile.bvp(func = freac2, bound = greac, parms = parms, 
  forcs = forcs, language = "Fortran")

## ------------------------
## Solution 1
## ------------------------
forcdata <- cbind(x = seq(0, 1, length.out = 4), R = rep(2, 4)) 

sola <- bvpcol(func = cfreac.F, x = seq(0, 1, by = 0.01), order = 2,
              leftbc = 1, ynames = c("y", "dy"), parms = parms, 
              forcings = forcdata)
fdata2 <- cbind(x = seq(0, 1, length.out = 4), R = 1:4) 
              
# second run, with changed forcings
sol2 <- bvpcol(func = cfreac.F, x = seq(0, 1, by = 0.01), order = 2,
              leftbc = 1, ynames = c("y", "dy"), parms = parms,
              forcings = fdata2)

plot(sol, sol2)


## =============================================================================
## a multipoint bvp
## dy1 = (y2 - 1)/2
## dy2 = (y1*y2 - x)/mu
## over interval [0,1]
## y1(1) = 0; y2(0.5) = 1
## =============================================================================

# The R-version

multip <- function (x, y, p) {
  list(c((y[2] - 1)/2, 
         (y[1]*y[2] - x)/mu))
}

bound <- function (i, y, p) {
  if (i == 1) y[2] -1    # at x=0.5: y2=1
  else y[1]              # at x=  1: y1=0
}

mu  <- 0.1
sol <- bvpcol(func = multip, bound = bound, 
              x = seq(0, 1, 0.01), posbound = c(0.5, 1))

# The compiled version

fmultip <- "
 f(1) = (y(2) - 1.)/2.
 f(2) = (y(1)*y(2) - x)/mu
"
fbound <- "
  if (i == 1) then
    g = y(2) -1 
  else 
    g = y(1)
  endif  
"
parms <- "mu"
cmultip <- compile.bvp(func = fmultip, bound = fbound, parms = "mu")

solf <- bvpcol(func = cmultip, parms = 0.1,
              x = seq(0, 1, 0.01), posbound = c(0.5, 1))

plot(sol, solf)  # same

# check boundary value
sol [sol [,1] == 0.5,]
solf[solf[,1] == 0.5,]

## =============================================================================
## a bvp DAE 
## =============================================================================

# R-version
bvpdae <- function(t, x, ks, ...) {
  p1  <- p2 <- sin(t)
  dp1 <- dp2 <- cos(t)
  
  dx1 <- (ks + x[2] - p2)*x[4] + dp1
  dx2 <- dp2
  dx3 <- x[4]
  res <- (x[1] - p1)*(x[4] - exp(t))

  list(c(dx1, dx2, dx3, res))
}

x <- seq(0, 1, by = 0.01)

# solved using yini, yend
out1 <- bvpcol (func = bvpdae, x = x, parms = 1e-4, 
                yini = c(x1 = sin(0), x2 = NA, x3 = 1, res = NA), 
                yend = c(NA, sin(1), NA, NA),
                dae = list(index = 2,  nalg = 1)) 

# the analytic solution
ana <- cbind(x, "x1" = sin(x), "x2" = sin(x), "x3" = 1, "x4" = 0, res = 0)

# compiled version
bvpdaef <- "
  double precision :: ks = 1.d-4  
  f(1) = (ks + y(2) - sin(x))*y(4) + cos(x)
  f(2) = cos(x)
  f(3) = y(4)
  f(4) = (y(1) - sin(x)) * (y(4) - exp(x))
"
cdae <- compile.bvp(func = bvpdaef)

cout <- bvpcol (func = cdae, x = x, nmax = 1e4, 
                yini = c(x1 = sin(0), x2 = NA, x3 = 1, res = NA), 
                yend = c(NA, sin(1), NA, NA), 
                dae = list(index = 2,  nalg = 1)) 
plot(out1, cout, obs = ana)

}

}
\author{
  Karline Soetaert <karline.soetaert@nioz.nl>
}

\keyword{ utilities }