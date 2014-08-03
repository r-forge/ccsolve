\name{snippets}
\docType{data}
\alias{tran1d}
\alias{tran2d}
\alias{tran2dmat}
\alias{print.code}
\title{
  Fortran code snippets that can be used when creating an inline model.
}

\description{
The definition of function \code{tran1d} is:

tran1d (N, C, Cup, Cdown, Fup, Fdown, BcUp, BcDown, K, v, dx ,dC) 
\itemize{
 \item INTEGER N                   ! length of C
 \item DOUBLE PRECISION C(N)       ! input concentration
 \item DOUBLE PRECISION Cup, Cdown ! Boundary concentrations (used if Bc = 2)
 \item DOUBLE PRECISION Fup, Fdown ! Boundary fluxes  (input if Bc = 1), else output
 \item DOUBLE PRECISION K, dx      ! Diffusion, grid size
 \item DOUBLE PRECISION v          ! velocity
 \item INTEGER BcUp, BcDown        ! boundary concitions (1= flux, 2=conc, 3 = 0-gradient) 
 \item DOUBLE PRECISION dC(N)      ! output: rate of change  (mass/sec)

}
The definition of function \code{tran2d} is:

tran2d (nx, ny, C, C_xup, C_xdown, C_yup, C_ydown, 
  F_xup, F_xdown, F_yup, F_ydown, BcxUp, BcxDown, Bcyup, Bcydown, 
  K_x, K_y, v_x, v_y, dx, dy, dC) 
\itemize{
 \item INTEGER nx, ny                          ! cells in x- and y
 \item DOUBLE PRECISION C(nx, ny)              ! input concentration
 \item DOUBLE PRECISION C_xup, C_xdown         ! Boundary concentrations in x (used if Bcx.. = 2)
 \item DOUBLE PRECISION C_yup, C_ydown         ! Boundary concentrations in y (used if Bcy.. = 2)
 \item DOUBLE PRECISION F_xup(ny), F_xdown(ny) ! Boundary fluxes in x (input if Bcx.. = 1), else output
 \item DOUBLE PRECISION F_yup(nx), F_ydown(nx) ! Boundary fluxes in y (input if Bcy.. = 1), else output
 \item DOUBLE PRECISION K_x, K_y               ! Diffusion
 \item DOUBLE PRECISION dx, dy                 ! Grid size
 \item DOUBLE PRECISION v_x, v_y               ! Velocity
 \item INTEGER BcxUp, BcxDown, BCyup, BCydown  ! boundary concitions (1= flux, 2=conc, 3 = 0-gradient) 
 \item DOUBLE PRECISION dC(nx, ny)             ! output: rate of change  (mass/sec)
}


tran2dmat (nx, ny, C, C_xup, C_xdown, C_yup, C_ydown, 
  F_xup, F_xdown, F_yup, F_ydown, BcxUp, BcxDown, Bcyup, Bcydown, 
  K_x, K_y, v_x, v_y, dx, dy, dxaux, dyaux, dC) 
\itemize{
 \item INTEGER nx, ny                          ! cells in x- and y
 \item DOUBLE PRECISION C(nx, ny)              ! input concentration
 \item DOUBLE PRECISION C_xup(ny), C_xdown(ny) ! Boundary concentrations in x (used if Bcx.. = 2)
 \item DOUBLE PRECISION C_yup(nx), C_ydown(nx) ! Boundary concentrations in y (used if Bcy.. = 2)
 \item DOUBLE PRECISION F_xup(ny), F_xdown(ny) ! Boundary fluxes in x (input if Bcx.. = 1), else output
 \item DOUBLE PRECISION F_yup(nx), F_ydown(nx) ! Boundary fluxes in y (input if Bcy.. = 1), else output
 \item DOUBLE PRECISION K_x(nx+1, ny)          !
 \item DOUBLE PRECISION K_y(nx, ny+1)          ! Diffusion on upstream interface
 \item DOUBLE PRECISION dx(nx), dxaux(nx+1)    !
 \item DOUBLE PRECISION dy(ny), dyaux(ny+1)    ! Grid size and distance from mid to mid
 \item DOUBLE PRECISION v_x(nx+1, ny)          ! Velocity in x, on upstream interface   
 \item DOUBLE PRECISION v_y(nx, ny+1)          ! Velocity in y, on upstream interface
 \item INTEGER BcxUp, BcxDown, BCyup, BCydown  ! boundary concitions (1= flux, 2=conc, 3 = 0-gradient) 
 \item DOUBLE PRECISION dC(nx, ny)             ! output: rate of change  (mass/sec)
}

}

\usage{
print.code(x)
data(tran1d)
}
\format{
  subroutines e.g. for transport that can be used as \code{includes} in the 
  compile functions
}

\author{
  Karline Soetaert <karline.soetaert@nioz.nl>
}
\examples{
print.code(tran1d)
}

\references{
}

\seealso{
}

\keyword{datasets}

