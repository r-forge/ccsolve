\name{compile.optim}
\alias{compile.optim}
\alias{ccoptim}

\title{Compiled code optimization}
\description{
  \code{compile.optim} generates compiled code 
  for multi-variable optimization. It should be solved with \code{ccoptim}.

  \code{ccoptim} provides a compiled code 
  extension to R's function \link{optim}.
 
}
\usage{
 compile.optim (func, jacfunc = NULL, data = NULL, par = NULL,  
      declaration = character(), includes = character(), language = "F95", ...)

 ccoptim(par, fn, gr = NULL, ..., method = c("Nelder-Mead",        
     "BFGS", "CG", "L-BFGS-B", "SANN"), lower = -Inf, upper = Inf,   
     control = list(), hessian = FALSE, data = NULL, 
     dllname = NULL, initfunc = NULL, rpar = NULL, 
     ipar = NULL) 
}
\arguments{
  \item{par }{Initial values for the parameters to be optimized over.
  }
  \item{fn }{A compiled subroutine that defines the function to be minimized 
    (or maximized), e.g. as generated by \code{compile.optim} or 
    a character string giving the name of a compiled subroutine in a 
    dynamically loaded shared library. 
  }
  \item{func }{A character vector with F95, Fortran, or C code, without 
    declarations, that specifies the optimization function. 
 
    The subroutine will be defined as \code{func(n, x, f, rpar, ipar)},
     where \code{x} is the vector of parameters over which minimization is 
     to take place, \code{n} is the number of parameters, \code{rpar} and 
     \code{ipar} are vectors of double and integer values, as passed with the 
     arguments of the same name to \code{ccoptim}. 
     The result, a scalar should be put in \code{f}.
  }
  \item{gr }{A compiled subroutine that specifies the gradient for the 
    "BFGS", "CG" and "L-BFGS-B" methods. 
    If it is NULL, a finite-difference approximation will be used.
    For the "SANN" method it specifies a function to generate a new candidate point. 
    If it is NULL a default Gaussian Markov kernel is used.
   }
  \item{jacfunc }{A character vector with F95, Fortran, or C code, without 
    declarations, that specifies the jacobian of the function whose optimum 
    has to be found. 
    The subroutine is defined as \code{jacfunc(n, x, df, rpar, ipar)},
    where \code{x} is the vector of parameters over which minimization is 
     to take place, \code{n} is the number of parameters, \code{rpar} and 
     \code{ipar} are vectors of double and integer values, as passed with the
     arguments of the same name. The result, a vector of length \code{n} 
     should be put in \code{df}.
  }
  \item{dllname }{A string giving the name of the shared library (without extension) 
   that contains all the compiled function or subroutine definitions refered to in 
   \code{fn} and \code{gr}.
  }   
  \item{initfunc }{if not NULL, the name of the initialisation function 
    (which initialises values of \code{data}), 
    as provided in "dllname". Will be generated with \code{compile.optim}, if 
    \code{data} is present. 
  }  
   \item{data }{A dataset with the same structure as will be used during run-time. 
     It should have named columns, and these names will be used to define global 
     variables (vectors) that will be given the actual data values upon run time.
     The names of these variables are available in the subroutines func and jacfunc.
     Also, the number of data values, passed when solving, will be known by
     a global variable called \code{ndat}.
     This is implemented for consistency with \code{ccnls}, and to facilitate using
     \code{ccoptim} for least-squares fitting - see last example.
   }
  \item{method }{The method to be used. See "Details" of \link{optim}.
   }
   \item{lower, upper }{Bounds on the variables for the "L-BFGS-B" method.
   }
   \item{control }{A list of control parameters. See "Details" of \link{optim}.
   }
   \item{hessian }{Logical. Should a numerically differentiated Hessian matrix 
     be returned?
   }  
   \item{rpar, ipar }{double and integer vector to be passed upon running the model.
   } 

  \item{declaration }{Text that eneters the declaration section in each function.
  }
  \item{includes }{Code that comes before the functions.
  }
  \item{language }{A character vector that specifies the source code; one of c("F95", "Fortran", "C") defaults to "F95".
  }
  \item{...}{optional arguments to the generic function (not used).
  }
}
\value{
  \code{compile.optim} returns an object of class \code{CFunc} or 
  \code{CFuncList}, as defined in the package  \code{inline}. 
}  

\details{

The compiled function that is generated by \code{compile.optim} is defined as:

 \code{func(n, x, f, rpar, ipar) }

 \code{jacfunc(n, x, df, rpar, ipar) }
 
In case Fortran or F95 is used, \code{n} is an integer value,

\code{ipar} and \code{rpar} are an integer vector, and double vector, that can be 
used to pass values at runtime.  

\code{x}, and \code{df} are a double vector of length \code{n}.

\code{f} is one value, the function to optimize.

In case C is used all are pointers. 

The values of \code{df} are set to 0 at the beginning of the subroutine.

The user needs to specify \code{f, df}, based on \code{x}

}

\note{

it is not possible to specify \code{parms} here (this was deemed to be too
confusing with the \code{par} argument. In case constants need to be passed
upon runtime, \code{rpar} should be used.
}
\examples{
\dontrun{

## =======================================================================
## example 1  
## small Rosenbrock Banana function (as from example page of optim)
## =======================================================================

fr <- function(x) {   ## Rosenbrock Banana function
    x1 <- x[1]
    x2 <- x[2]
    100 * (x2 - x1 * x1)^2 + (1 - x1)^2
}
grr <- function(x) { ## Gradient of 'fr'
    x1 <- x[1]
    x2 <- x[2]
    c(-400 * x1 * (x2 - x1 * x1) - 2 * (1 - x1),
       200 *      (x2 - x1 * x1))
}

## ---------------------------------------------------------------------------
## The F95 version 
## ---------------------------------------------------------------------------

fr.f95 <- "
  f = 100.d0 * (x(2) - x(1) * x(1))**2 + (1. - x(1))**2"

grr.f95 <- "
  df(1) = -400.d0 * x(1) *(x(2) - x(1) * x(1)) -2.d0* (1. - x(1))
  df(2) = 200d0          *(x(2) - x(1) * x(1))"

ccfr  <- compile.optim(func = fr.f95, jacfunc = grr.f95)

optim(c(-1.2,1), fr)
ccoptim(c(-1.2, 1), ccfr)
optim(c(-1.2,1), fr, grr, method = "BFGS", hessian = TRUE)
ccoptim(c(-1.2,1), ccfr, method = "BFGS", hessian = TRUE)

## These do not converge in the default number of steps
optim(c(-1.2,1), fr, grr, method = "CG")
ccoptim(c(-1.2,1), ccfr, method = "CG")
optim(c(-1.2,1), fr, grr, method = "CG", control = list(type = 2))
ccoptim(c(-1.2,1), ccfr, method = "CG", control = list(type = 2))
optim(c(-1.2,1), fr, grr, method = "L-BFGS-B")
ccoptim(c(-1.2,1), ccfr, method = "L-BFGS-B")

# retrieve function value for initial estimates and 'best'
ccfunc(ccfr, c(-1.2, 1))
ccfunc(ccfr, c(1.00026, 1.000506))

## =======================================================================
## example 2 -  a larger Rosenbrock function  (as from optimx package)
## =======================================================================

genrose.f <- function(x, gs=NULL){ # objective function
	n <- length(x)
  ii <- seq(1, n-1, by = 2)
 	fval <- sum (gs*(x[ii+1]^2 - x[ii])^2 + (x[ii] - 1)^2)
        return(fval)
}

startx <- 4*seq(1:1000)/3.

print(system.time(ans <- optim(startx, fn = genrose.f, gs=100)))

genrose.f95 <- "
  integer i
  f = 0.d0
  do i = 1, n-1, 2
    f = f + rpar(1) * (x(i+1)**2 - x(i))**2 + (x(i) - 1)**2
  enddo
  f = 1.d0 + f
"
cGenrose <- compile.optim(genrose.f95)
print(system.time(anscc <- ccoptim(startx, fn = cGenrose, rpar = 100)))

## =======================================================================
## example 3: from nls - see also help file of ccnls
## =======================================================================

DNase1 <- subset(DNase, Run == 1)

fm2DNase1 <- nls(density ~ 1/(1 + exp((xmid - log(conc))/scal)),
                 data = DNase1,
                 start = list(xmid = 0, scal = 1))
summary(fm2DNase1)

# compile it such that the compiler knows the names of the par to be solved for 
# and the names of the data; the values of par and data can differ during
# the actual application; the ordering of parameters and datacolumns should
# stay the same.

# F95 works with vectors
head (DNase1 [, -1])      # names conc, density

fo = "f = sum((density - 1.0/(1.d0 + dexp((xmid - dlog(conc))/scal)))**2)"
ccDNasef <- compile.optim(func = fo, par = c(xmid = 0, scal = 1), data = DNase1[,-1])  
printCode(ccDNasef)

# for C you need to write a loop ; length of data = ndat :
ccDNase.C <- compile.optim(func = '
 int i;
 double ff; 
 *f = 0.;
 for (i = 0; i < ndata; i++){
   ff = density[i] - 1.0/(1.0 + exp((xmid - log(conc[i]))/scal));
   *f = *f + ff*ff;
 }  
   ',   
 par = c(xmid = 0, scal = 1), 
 data = DNase1[,-1], language = "C")  

fmDNase2 <- ccoptim(f = ccDNasef, data = DNase1[,-1], 
  par = c(xmid = 0, scal = 1))

fmDNase3 <- ccoptim(f = ccDNase.C, data = DNase1[,-1], 
  par = c(xmid = 0, scal = 1))

fmDNase2
}
}
\author{
  Karline Soetaert <karline.soetaert@nioz.nl>

  Function \code{ccoptim} has copied significant parts of the R-function \link{optim},
  and the C-code in file \code{optim.c} from base R (stats).
}


\keyword{ utilities }
