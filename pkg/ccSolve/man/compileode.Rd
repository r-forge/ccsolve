\name{compile.ode}
\alias{compile.ode}

\title{Creates compiled code for the deSolve functions }
\description{
  Uses inline compilation to generate compiled code from the body of 
  fortran or C function that specify an initial value problem or ordinary
  differential equations.
}
\usage{
 compile.ode (func, jacfunc = NULL, rootfunc = NULL, eventfunc = NULL,
      parms = NULL, forcings = NULL, outnames = NULL, y = NULL, 
      declaration = character(), includes = character(), language = "F95", ...)
}
\arguments{
  \item{func }{A character vector with F95, Fortran, C or C++ code, without 
    declarations, that specifies the derivative function. 
    The function is called with arguments \code{n, t, y, f, rpar, ipar} and the 
    derivative should be written to vector \code{f}.
    The variables n (integer), t, and y (double) should not be changed in the code.
    
  }
  \item{jacfunc }{A character vector with F95, Fortran, C or C++ code, without 
    declarations, that specifies the jacobian of the derivative function.
    The function is called with arguments \code{n, t, y, df, rpar, ipar} and the 
    jacobian should be written to matrix (language = Fortran, F95) or 
    vector (language = C) \code{df}.
  }
  \item{rootfunc }{A character vector with F95, Fortran, C or C++ code, without 
    declarations, that specifies the root function.
    The function is called with arguments \code{n, t, y, nroot, root, rpar, ipar} and the 
    constraint function whose root is sought should be written to vector \code{root}.
  }
  \item{eventfunc }{A character vector with F95, Fortran, C or C++ code, without 
    declarations, that specifies the event function.
    The function is called with arguments \code{n, t, y} and the code needs to 
    redefine (elements) of vector \code{y}.
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
  \item{y }{The initial condition vector or a vector that contains the names of the state variables used in the equations. 
    These names will be known in the functions; also the derivatives will be declared in the function; if a state varialbe is called sv, 
    then its corresponding derivative will be called dsv; the user must specify the derivative in the code.
    This vector should be compatible with the \code{y} vector used in the call to the solver.
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
  An object of class \code{CFunc} or \code{CFuncList}, as from the package  \code{inline}. 
  
  Each object of class \code{CFunc} can be called with the appropriate arguments
}  
\details{

The compiled functions are compatible with the C-functions as used in the solvers.

They are defined as:

 \code{func(n, t, y, f, rpar, ipar) }

 \code{jacfunc(n, t, y, df, rpar, ipar) }
 
 \code{rootfunc(n, t, y, nroot, root, rpar, ipar) }
 
 \code{eventfunc(n, t, y)}

In case Fortran or F95 is used, \code{n, i} are an integer value, \code{t} is a double value (typically time), 

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

\note{


}
\examples{
\dontrun{

## =======================================================================
##
## ex. 1
## The famous Lorenz equations: chaos in the earth's atmosphere
## Lorenz 1963. J. Atmos. Sci. 20, 130-141.
##
## =======================================================================

# ---------------------
# the R-version
# ---------------------
chaos <- function(t, state, parameters) {
  with(as.list(c(state)), {

    dx     <- -8/3 * x + y * z
    dy     <- -10 * (y - z)
    dz     <- -x * y + 28 * y - z

    list(c(dx, dy, dz))
  })
}

# ---------------------
# The fortran version
# ---------------------
Chaos.f <- " 
   f(1) = -8.d0/3 * y(1) + y(2)*y(3)
   f(2) = -10.d0*(y(2)-y(3))
   f(3) = -y(1)*y(2) + 28.d0 * y(2) - y(3)
"
cChaos <- compile.ode(Chaos.f) 

# ---------------------
# + parameters, using C
# ---------------------
Chaos.C <- " 
   f[0] = a * y[0] + y[1]*y[2];
   f[1] = b*(y[1]-y[2]);
   f[2] = -y[0]*y[1] + c * y[1] - y[2];
"
parms <- c(a = -8.0/3, b = -10.0, c = 28.0)
cChaos2 <- compile.ode(Chaos.C, language = "C", parms = parms) 

# ---------------------
# running all  
# ---------------------

state <- c(x = 1, y = 1, z = 1)
times <- seq(0, 100, 0.01)

print(system.time(out   <- vode(state, times, chaos, 0)))
print(system.time(cout   <- vode(state, times, func = cChaos, parms = 0)))
print(system.time(cout2  <- vode(state, times, func = cChaos2, parms = parms)))

plot(out, cout, cout2)    
plot(out[,"x"], out[,"y"], type = "l", main = "Lorenz butterfly",
  xlab = "x", ylab = "y")


## =======================================================================
##
## ex. 2 - a difference equation : f() has new value of y 
## The parasite model from vignette ("compiledCode") 
##
## =======================================================================

parms <- c(rH = 2.82, A = 100, ks = 15)

parasite.R <- function (t, y, parms) {
  with (as.list(parms), {
   P <- y[1]
   H <- y[2]
   f <- A * P / (ks +H)
   Pnew <- H* (1-exp(-f))
   Hnew <- H * exp(rH*(1.-H) - f)
   list (c(Pnew, Hnew))   
  })
}

# ---------------------
# fortran, using parms
# ---------------------

parasite.f <- "
	      double precision P, H, fr

    	  P = y(1)
    	  H = y(2)
        fr = A * P / (ks + H)

        f(1) = H * (1.d0 - exp(-fr))
        f(2) = H * exp (rH * (1.d0 - H) - fr)
"
cparasite <- compile.ode(func = parasite.f, parms = parms)

# ---------------------
# C, using parms and y
# the compiling function maps the variables
# ---------------------

parasite.C <- "
        double fr = A * P / (ks + H);

        dP = H * (1.0 - exp(-fr));
        dH = H * exp (rH * (1.0 - H) - fr);
"
yini <- c(P = 0.5, H = 0.5)
cparasite.c <- compile.ode(func = parasite.C, parms = parms, 
                           y = yini, language = "C")

# ---------------------
# solution
# ---------------------
system.time(out <- ode (func = parasite.R, y = yini, parms = parms, times = 0:500,
    method = "iteration"))
    
system.time(outc <- ode (func = cparasite, y = yini, parms = parms, times = 0:500,
     method = "iteration"))
plot(out, outc)

## ===========================================================================
##
## ode with event in a function
## derivative function: rate of change v1 = 0, v2 reduced at first-order rate
## events: add 1 to v1, multiply v2 with 0.5
##
## ===========================================================================

derivs.R <- function(t, var, parms) {
   list(c(0, -0.5 * var[2]))
}

event.R <- function(t, y, parms){
  with (as.list(y),{
    v1 <- v1 + 1
    v2 <- 5 * 0.5
    return(c(v1, v2))
  })
}

# ---------------------
# compiled code 
# ---------------------

derivs.f <- " 
 f(1) = 0.d0
 f(2) = -0.5d0 * y(2)
"
event.f <- "
 y(1) = y(1) + 1.d0
 y(2) = 5 * 0.5d0
"

cder <- compile.ode(func = derivs.f, eventfunc = event.f)

# ---------------------
# running
# ---------------------

yini <- c(v1 = 1, v2 = 2)
times <- seq(0, 10, by = 0.1)

out <- ode(func = derivs.R, y = yini, times = times, parms = NULL,
           events = list(func = event.R, time = 1:9) )
out2 <- ode(func = cder, y = yini, times = times, parms = NULL,
           events = list(time = 1:9) )

plot(out, out2, type = "l")

## ===========================================================================
##
## ODE with a root  + output variable
## Example from lsodar source code
##
## =======================================================================

Fun.R <- function (t, y, parms) {
  ydot <- vector(len = 3)
  ydot[1] <- -.04*y[1] + 1.e4*y[2]*y[3]
  ydot[3] <- 3.e7*y[2]*y[2]
  ydot[2] <- -ydot[1] - ydot[3]
  return(list(ydot, Sum = sum(y)))
}

rootFun.R <- function (t, y, parms) {
  yroot <- vector(len = 2)
  yroot[1] <- y[1] - 1.e-4
  yroot[2] <- y[3] - 1.e-2
  return(yroot)
}

# ---------------------
# Fortran
# ---------------------

Fun.f<- "
  f(1) = -.04d0*y(1) + 1.d4*y(2)*y(3)
  f(3) = 3.d7*y(2)*y(2)
  f(2) = -f(1) - f(3)
  Sum = y(1) + y(2) + y(3)
"

Root.f <- "
  root(1) = y(1) - 1.d-4
  root(2) = y(3) - 1.d-2
"
cRoot <- compile.ode(func = Fun.f, rootfunc = Root.f, outnames = "Sum")

# ---------------------
# application
# ---------------------

y     <- c(1, 0, 0)
times <- c(0, 0.4*10^(0:8))
Out   <- NULL
ny    <- length(y)

out   <- lsodar(y = y, times = times, func = Fun.R, rootfunc = rootFun.R,
                rtol = 1e-4, atol = c(1e-6, 1e-10, 1e-6), parms = NULL)
outc  <- lsodar(y = y, times = times, func = cRoot, nroot = 2, 
                outnames = "Sum", nout = 1,
                rtol = 1e-4, atol = c(1e-6, 1e-10, 1e-6), parms = NULL)

print(paste("root is found at t", attributes(outc)$troot))
print(paste("root is found at t", attributes(out)$troot))


## ===========================================================================
##
## ode with a root function and event
## the bouncing ball
## event triggered when the ball hits the ground (height = 0)
## bouncing means a change in spped direction and loss in speed
##
## =============================================================================

#-----------------------------
# the R implementation
#-----------------------------
ballode.R <- function(t, y, parms) {
  dy1 <- y[2]
  dy2 <- -9.8
  list(c(dy1, dy2))
}

root.R <- function(t, y, parms) 
  y[1]

event.R <- function(t, y, parms) {
  y[1] <- 0
  y[2] <- -0.9 * y[2]
 return(y)
}

# ---------------------
# in fortran
# ---------------------
ballode.f <- "
 f(1) = y(2)
 f(2) = -9.8
"

root.f  <-  " 
 root(1) = y(1) 
" 

event.f <- "
 y(1) = 0.d0
 y(2) = -0.9 * y(2)
"
cballode <- compile.ode(func = ballode.f, rootfunc = root.f, eventfunc = event.f)

# ---------------------
# application
# ---------------------

yini  <- c(height = 0, v = 20)
times <- seq(0, 40, 0.01)

out   <- lsodar(times = times, y = yini, func = ballode.R, parms = NULL,
  events = list(func = event.R, root = TRUE), rootfun = root.R)
out2  <- lsodar(times = times, y = yini, func = cballode, parms = NULL,
  events = list(root = TRUE), nroot = 1) 

attributes(out)$troot
attributes(out2)$troot

#-----------------------------
# display, plot results
#-----------------------------
plot(out, out2, which = "height", type = "l", lwd = 2, 
  main = "bouncing ball", ylab = "height")

## ===========================================================================
##
## An ode with a forcing function
## the SCOC model
##
## =============================================================================

# ---------------------
# Forcing function data
# ---------------------
Flux <- matrix(ncol = 2, byrow = TRUE, data = c(
  1, 0.654, 11, 0.167,   21, 0.060, 41, 0.070, 73,0.277, 83,0.186,
  93,0.140,103, 0.255,  113, 0.231,123, 0.309,133,1.127,143,1.923,
  153,1.091,163,1.001,  173, 1.691,183, 1.404,194,1.226,204,0.767,
  214, 0.893,224,0.737, 234,0.772,244, 0.726,254,0.624,264,0.439,
  274,0.168,284 ,0.280, 294,0.202,304, 0.193,315,0.286,325,0.599,
  335, 1.889,345, 0.996,355,0.681,365,1.135))

parms <- c(k = 0.01)

# ---------------------
#  model in R-code
# ---------------------
sediment.R <- function( t, state, parms) {
 with (as.list(c(parms, state)), {
  min <- k * O2
  depo <- Depo(t) 
  list (depo - min, Flux = depo, Min = min)
 }) 
}

# the forcing functions approx function
Depo <- approxfun(x = Flux[,1], y = Flux[,2])

# ---------------------
# model in FORTRAN
# ---------------------
sediment.f <- "
  f(1) = Depo - k*y(1)
  Flux = Depo
  Min = k*y(1)
"
outnames <- c("Flux", "Min")
cSed <- compile.ode(func = sediment.f, forcings = "Depo", 
  parms = parms, outnames = outnames)

# ---------------------
# and in FC
# ---------------------
sediment.c <- "
  f[0] = Depo - k*y[0];
  Flux = Depo ;
  Min = k*y[0];
"
outnames <- c("Flux", "Min")
# 
cSed2 <- compile.ode(func = sediment.c, forcings = "Depo",  
  parms = parms, outnames = outnames,language = "C")

# ---------------------
# application
# ---------------------

times <- 1:365
Out  <- ode(times = times, func = sediment.R, y = c(O2 = 63), parms = parms)
Out2 <- ode(times = times, func = cSed, y = c(O2 = 63), 
  outnames = outnames, nout = 2, parms = parms, forcings = Flux)
Out3 <- ode(times = times, func = cSed2, y = c(O2 = 63), 
  outnames = outnames, nout = 2, parms = parms, forcings = Flux)

plot(Out, Out2)
}
}
\author{
  Karline Soetaert <karline.soetaert@nioz.nl>
}

\keyword{ utilities }