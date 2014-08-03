\name{compile.pde}
\alias{compile.pde}
\title{
  compiles partial differential equation problems
}

\usage{
print.code(x)
}
\arguments{
  \item{x }{Code snippit to be printed.
  }

}

\author{
  Karline Soetaert <karline.soetaert@nioz.nl>
}
\examples{
print.code(tran1d)

## ===========================================================================
## a simple 1D transport model
## ===========================================================================

declaration <- "   double precision :: Fin, Fout"

transport <- "

    ! tran1d    (nx,C, Cup, Cdown, Fin, Fout,BcUp,BcDown,K,v,dx,dC) 
    CALL tran1d (nx,C, Cup, Cdown, Fin, Fout, 2, 2, D, v, dx, dC) 

    ! reaction
    dC = dC - k*C

"

parms <- c(k = 0.05, Cup = 20, Cdown = 30, D = 1, v = 1, dx = 1)
ctran <- compile.pde(func = transport, parms = parms, 
   includes = tran1d, declaration = declaration, 
   names = "C", dim = 100)

out <- ode.1D(y = rep(0,100), parms = parms, func = ctran, times = 0:365,
   names = "C", dim = 100)
std <- steady.1D(y = rep(0,100), parms = parms, func = ctran,
    names = "C", dim = 100)
image(out)

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
  names = c("prey", "predator"), dim = nx,  
  parms = parms, nspec = 2, method = lsodes)
))

image(out)


## ===========================================================================
## a two-species transport model in 2-D
## ===========================================================================

declaration = " 
           DOUBLE PRECISION :: Ingestion(nx, ny), GrowthPrey(nx, ny), MortPredator(nx, ny)
           DOUBLE PRECISION :: Fxin(ny), Fxout(ny), Fyin(nx), Fyout(nx)"

LV2D <- "

 ! transport
   CALL tran2d (nx, ny, prey, 2.d0, 2.d0, 2.d0, 2.d0,      &
                 Fxin, Fxout, Fyin, Fyout, 2, 2, 2, 2,     &
                 Kx, Ky, vx, vy, dx, dy, dprey) 
   
   CALL tran2d (nx, ny, Predator, 1.d0, 1.d0, 1.d0, 1.d0,  &
                 Fxin, Fxout, Fyin, Fyout, 2, 2, 2, 2,     &
                 Kx, Ky, vx, vy, dx, dy, dpredator)  

 ! biology
   Ingestion    = rIng  * Prey * Predator
   GrowthPrey   = rGrow * Prey * (1 - Prey/Kprey)
   MortPredator = rMort * Predator
 
   dPrey        = dprey + GrowthPrey - Ingestion
   dPredator    = dpredator + Ingestion * assEff - MortPredator

"

nx <- 50
ny <- 50
parms <- c(Kx = 10, Ky = 10, vx = 1, vy = 1, dx = 0.5, dy = 0.5, 
           rIng = 0.2, rGrow = 1.0, rMort = 0.2, assEff = 0.5, Kprey = 10)

cpredprey1D <- compile.pde(func = LV2D, parms = parms, includes = tran2d, 
  names = c("prey", "predator"), dim = c(nx, ny), 
  declaration = declaration)

yini  <- c(rep(1, nx*ny), rep(2, nx*ny))

print(system.time(
out   <- steady.2D(y = yini, pos = TRUE, func = cpredprey1D,
  names = c("prey", "predator"), dim = c(nx, ny), 
  parms = parms, nspec = 2, method = lsodes)
))

image(out)


}

\references{
}

\seealso{
}

\keyword{datasets}

