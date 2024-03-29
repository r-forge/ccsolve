\name{compile.ode}
\alias{compile.ode}
\alias{compile.steady}

\title{Creates compiled code for ordinary differential equation solvers }
\description{
  Generates compiled code from the body of fortran subroutines or C functions 
  that specify an initial value problem or ordinary differential equations (ODE).  
  These problems can then be solved with differential equation functions from 
  the R-packages \code{deSolve}, or \code{deTestSet} or to find the 
  steady-state solution of ODEs as from the R-package \code{rootSolve}.
}  
\usage{
 compile.ode (func, jacfunc = NULL, rootfunc = NULL, eventfunc = NULL,
      parms = NULL, y = NULL, forcings = NULL, outnames = NULL, 
      declaration = character(), includes = character(), language = "F95", ...)
 compile.steady (func, jacfunc = NULL, rootfunc = NULL, eventfunc = NULL,
      parms = NULL, y = NULL, forcings = NULL, outnames = NULL, 
      declaration = character(), includes = character(), language = "F95", ...)
}
\arguments{
  \item{func }{A character vector with F95, Fortran, or C code, without 
    declarations, that specifies the derivative function. 
    The function is called with arguments \code{n, t, y, f, rpar, ipar} and the 
    derivatives should be written to vector \code{f}. 
    (see note for difference equations).
    The variables \code{n}, \code{t}, and \code{y} should not be changed in the code.
    \code{rpar} and \code{ipar} are a double and integer vector, passed upon runtime
    to the solver. See vignette \code{"compiledCode"} from R-package \code{deSolve}.
  }
  \item{jacfunc }{A character vector with F95, Fortran, or C code, without 
    declarations, that specifies the jacobian of the derivative function.
    The function is called with arguments \code{n, t, y, df, rpar, ipar} and the 
    jacobian should be written to matrix (language = Fortran, F95) or 
    vector (language = C) \code{df}.
    The variables \code{n}, \code{t}, and \code{y} should not be changed in the code.
  }
  \item{rootfunc }{A character vector with F95, Fortran, or C code, without 
    declarations, that specifies the root function.
    The function is called with arguments \code{n, t, y, nroot, root, rpar, ipar} and the 
    constraint function whose root is sought should be written to vector \code{root}.
    The variables n, t, y and nroot should not be changed in the code.
  }
  \item{eventfunc }{A character vector with F95, Fortran, or C code, without 
    declarations, that specifies the event function.
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
    (see e.g. \link{ode}).
  }
  \item{y }{The initial condition vector or a vector that contains the names of 
    the state variables used in the equations. 
    These names will be known in the compiled code. 
    In addition the derivatives will be declared in the function, as follows: for a state varialbe called sv, 
    then its corresponding derivative will be called dsv. 
    The user must specify these derivatives in the code.
    This vector should be compatible with the \code{y} vector used in the call to the solver 
    (see e.g. \link{ode}).
  }  
  \item{forcings }{A vector with the names of the forcings to be used in the equations. 
    These names will be declared in the compiled code, 
    and their values updated for each \code{t} when solving the model.
    This vector should be compatible with the \code{forcings} vector used in the call to the solver 
    (see e.g. \link{ode}).
  }
  \item{outnames }{The names of the output variables used in the equations. 
    These names will be known in the compiled code, and the user needs to specify their values.
    These values will then be used to update the ouput.
    This vector should be compatible with the \code{outnames} vector used in the call to the solver
    (see e.g. \link{ode}).
  }
  \item{declaration }{Text that enters the declaration section in each function.
  }
  \item{includes }{Code that comes before the functions (e.g. extra function or subroutine used, global variables (C) etc...
  }
  \item{language }{A character vector that specifies the source code; one of c("F95", "Fortran", "C"); defaults to "F95".
  }
  \item{...}{optional arguments to the function (not used).
  }
}
\value{
  An object of class \code{CFunc} or \code{CFuncList}, as defined in the package  \code{inline}. 
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

In case C is used all are pointers. 

the values of \code{df} are set to 0 at the beginning of the subroutine.

The variables \code{n}, \code{t}, and -except for eventfunc- \code{y}, 
should not be changed within the function; they are input only.
The user needs to specify \code{f, df, root} for \code{func, jacfunc, rootfunc}, 
and \code{y} for \code{eventfunc} respectively.

}

\note{
Function \code{compile.steady} is a synonym for \code{compile.ode}.

In case the model to be solved is a difference equation 
(\code{ode(method = "iteration", ...)}), the new value of the state variables 
should be put in \code{f}, not the derivative - see second example.

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
# The F95 version 
# ---------------------                        
# note 8.d0 assures double precision is used
Chaos.f <- " 
   f(1) = -8.d0/3 * y(1) + y(2)*y(3)
   f(2) = -10.d0*(y(2)-y(3))
   f(3) = -y(1)*y(2) + 28.d0 * y(2) - y(3)
"
cChaos <- compile.ode(Chaos.f) 
# Show the complete code
printCode(cChaos)

# ---------------------
# + parameters, C code
# ---------------------
Chaos.C <- " 
   f[0] = a * y[0] + y[1]*y[2];
   f[1] = b*(y[1]-y[2]);
   f[2] = -y[0]*y[1] + c * y[1] - y[2];
"
parms <- c(a = -8.0/3, b = -10.0, c = 28.0)
cChaos2 <- compile.ode(Chaos.C, language = "C", parms = parms) 
printCode(cChaos2)

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
    	  P = y(1)
    	  H = y(2)
        fr = A * P / (ks + H)

        f(1) = H * (1.d0 - exp(-fr))
        f(2) = H * exp (rH * (1.d0 - H) - fr)
"
cparasite <- compile.ode(func = parasite.f, parms = parms, 
  declaration = "double precision P, H, fr")

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
system.time(outc2 <- ode (func = cparasite.c, y = yini, parms = parms, times = 0:500,
     method = "iteration"))

plot(out, outc, outc2)

## ===========================================================================
##
## ex. 3. ode with event in a function
## derivative function: dv1 = 0, dv2 = -0.5*v2 
## events: add 1 to v1, multiply v2 with 0.5
##
## ===========================================================================

derivs.R <- function(t, var, parms) {
   list(c(0, -0.5 * var[2]))
}

event.R <- function(t, y, parms){
  with (as.list(y),{
    v1 <- v1 + 1
    v2 <- v2 * 0.5
    return(c(v1, v2))
  })
}

# ---------------------
# f95 code 
# ---------------------

derivs.f <- " 
 f(1) = 0.d0
 f(2) = -0.5d0 * y(2)
"
event.f <- "
 y(1) = y(1) + 1.d0
 y(2) = y(2) * 0.5d0
"

cder <- compile.ode(func = derivs.f, eventfunc = event.f)

# ---------------------
# C code 
# ---------------------

derivs.c <- " f[0] = 0.; f[1] = -0.5 * y[1];"
event.c <- " y[0] = y[0] + 1.; y[1] = y[1] * 0.5;"

cder.c <- compile.ode(func = derivs.c, eventfunc = event.c, 
  language = "C")

# ---------------------
# running
# ---------------------

yini <- c(v1 = 1, v2 = 2)
times <- seq(0, 10, by = 0.1)

out <- ode(func = derivs.R, y = yini, times = times, parms = NULL,
           events = list(func = event.R, time = 1:9) )
out2 <- ode(func = cder, y = yini, times = times, parms = NULL,
           events = list(time = 1:9) )
out3 <- ode(func = cder.c, y = yini, times = times, parms = NULL,
           events = list(time = 1:9) )

plot(out, out2, out3, type = "l")

## ===========================================================================
##
## ex. 4. ODE with a root  + output variable
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
# F95
# ---------------------

Fun.f <- "
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
## ex 5. ode with a root function and an event
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
# in F95
# ---------------------
ballode.f <- "
 f(1) = y(2)
 f(2) = -9.8d0
"

root.f  <-  " 
 root(1) = y(1) 
" 

event.f <- "
 y(1) = 0.d0
 y(2) = -0.9d0 * y(2)
"
cballode <- compile.ode(func = ballode.f, rootfunc = root.f, eventfunc = event.f)

# ---------------------
# application
# ---------------------

yini  <- c(height = 0, v = 20)
times <- seq(0, 40, 0.01)

print(system.time(
  out   <- lsodar(times = times, y = yini, func = ballode.R, parms = NULL,
    events = list(func = event.R, root = TRUE), rootfun = root.R)
))
    
print(system.time(
  out2  <- lsodar(times = times, y = yini, func = cballode, parms = NULL,
    events = list(root = TRUE), nroot = 1) 
))

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
## =========================================================================== 

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
# and in C
# ---------------------
sediment.c <- "
  f[0] = Depo - k*y[0];
  Flux = Depo;
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
print(system.time(
  Out  <- ode(times = times, func = sediment.R, y = c(O2 = 63), 
    parms = parms)
))
print(system.time(
  Out2 <- ode(times = times, func = cSed, y = c(O2 = 63), 
   outnames = outnames, nout = 2, parms = parms, forcings = Flux)
))
print(system.time(
  Out3 <- ode(times = times, func = cSed2, y = c(O2 = 63), 
    outnames = outnames, nout = 2, parms = parms, forcings = Flux)
))
plot(Out, Out2, Out3)

## =======================================================================
##
##  Steady-state application (from example(steady))
##  Bacteria (Bac) growing on a substrate (Sub)
##  Example from Soetaert and Herman (2009)
##
## =======================================================================

# R-function  
model <- function(t, state, pars) {
  with (as.list(c(state,pars)), {
  #       substrate uptake             death  respiration
  dBact = gmax*eff*Sub/(Sub+ks)*Bact - dB*Bact - rB*Bact
  dSub  =-gmax    *Sub/(Sub+ks)*Bact + dB*Bact          +input
  
  return(list(c(dBact, dSub)))
 })
}
 
pars <- list(gmax = 0.5,eff = 0.5,
             ks = 0.5, rB = 0.01, dB = 0.01, input = 0.1)
yini <- c(Bact = 0.1, Sub = 0)

# C version
model.C <- "
  dBact = gmax*eff*Sub/(Sub+ks)*Bact - dB*Bact - rB*Bact;
  dSub  =-gmax    *Sub/(Sub+ks)*Bact + dB*Bact          +input;
"
cmodel <- compile.steady( func = model.C, parms = pars, 
  y = yini, language = "C")

# Newton-Raphson
steady(y = yini, time = 0,
       func = model, parms = pars, pos = TRUE)  
steady(y = yini, time = 0,
       func = cmodel, parms = pars, pos = TRUE)  

# Dynamic run to steady-state
steady(y = c(Bact = 0.1, Sub = 0), time = c(0, 1e5),
             func = cmodel, parms = pars, method = "runsteady")
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

Karline Soetaert and Peter.M.J. Herman (2009).  
A Practical Guide to Ecological Modelling. Using R as a Simulation Platform.  
Springer, 372 pp.

Karline Soetaert (2014). rootSolve: Nonlinear root finding, equilibrium and 
steady-state analysis of ordinary differential equations.  R-package
  version 1.7

Karline Soetaert, Jeff Cash and Francesca Mazzia (2014). 
deTestSet: Testset for differential equations. R package version 1.1.1.
}

\keyword{ utilities }