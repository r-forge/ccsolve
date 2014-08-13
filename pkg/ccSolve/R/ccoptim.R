# ------------------------------------------------------------------------------
#
#  Adaptation to compiled code of the functions optim and optimize
#
#  Copyright (C) optim and optimize 2000-12 The R Core Team
#
#  Adaptation for use with compiled code by Karline Soetaert
#
# ------------------------------------------------------------------------------

isvalid <- function(f) {        # check if a function points to compiled code
  is.character(f) | class(f) == "CFunc"
}


ccoptim <-
    function(par, fn, gr = NULL, ...,
             method = c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN", "Brent"), 
             lower = -Inf, upper = Inf, 
             control = list(), hessian = FALSE, 
             dllname = NULL, rpar = NULL, ipar = NULL)
{

  if (is.list(fn)) {            ### IF a list
      if (!is.null(gr) & "jacfunc" %in% names(fn))
         stop("If 'fn' is a list that contains 'jacfunc', argument 'gr' should be NULL")
      if (!is.null(dllname) & "dllname" %in% names(fn))
         stop("If 'fn' is a list that contains dllname, argument 'dllname' should be NULL")

     gr <- fn$jacfunc  
     dllname <- fn$dllname
     fn <- fn$func
  }

    if (! isvalid(fn))
      stop("'fn' should be either a compiled function or a character vector")

    if (! is.null(gr))
      if (! isvalid(gr))
        stop("'gr' should be either a compiled function or a character vector")

    if (is.character(fn) & is.null(dllname))  
      stop("'dllname' should have a value if 'fn' is a character")

    if (sum(duplicated (c(fn, gr))) > 0)
      stop("fn, and gr cannot be the same")


      if (class (fn) == "CFunc")
        fn1 <- body(fn)[[2]]
      else if (is.loaded(fn, PACKAGE = dllname, type = "") ||
        is.loaded(fn, PACKAGE = dllname, type = "Fortran"))  {
        fn1 <- getNativeSymbolInfo(fn, PACKAGE = dllname)$address
      } else 
        stop(paste("'fn' not loaded ", fn))

    gr1 <- NULL
    if (!is.null(gr)) {
      if (class (gr) == "CFunc")
        gr1 <- body(gr)[[2]]
      else if (is.loaded(gr, PACKAGE = dllname, type = "") ||
        is.loaded(gr, PACKAGE = dllname, type = "Fortran"))  {
        gr1 <- getNativeSymbolInfo(gr, PACKAGE = dllname)$address
      } else 
        stop(paste("'gr' not loaded ", gr))
    }

    method <- match.arg(method)
    if((length(lower) > 1L || length(upper) > 1L ||
       lower[1L] != -Inf || upper[1L] != Inf)
       && !any(method == c("L-BFGS-B","Brent"))) {
	warning("bounds can only be used with method L-BFGS-B (or Brent)")
	method <- "L-BFGS-B"
    }
    npar <- length(par)
    
    if (is.null(rpar))
      rpar <- 0.
    if (is.null(ipar))
      ipar <- 0  
    
    ## Defaults :
    con <- list(trace = 0, fnscale = 1, parscale = rep.int(1, npar),
		ndeps = rep.int(1e-3, npar),
		maxit = 100L, abstol = -Inf, reltol = sqrt(.Machine$double.eps),
		alpha = 1.0, beta = 0.5, gamma = 2.0,
		REPORT = 10,
		type = 1,
		lmm = 5, factr = 1e7, pgtol = 0,
		tmax = 10, temp = 10.0)
    nmsC <- names(con)
    if (method == "Nelder-Mead") con$maxit <- 500
    if (method == "SANN") {
	con$maxit <- 10000
	con$REPORT <- 100
    }
    con[(namc <- names(control))] <- control
    if(length(noNms <- namc[!namc %in% nmsC]))
	warning("unknown names in control: ", paste(noNms,collapse=", "))
    if(con$trace < 0)
	warning("read the documentation for 'trace' more carefully")
    else if (method == "SANN" && con$trace && as.integer(con$REPORT) == 0)
	stop("'trace != 0' needs 'REPORT >= 1'")
    if (method == "L-BFGS-B" &&
	any(!is.na(match(c("reltol","abstol"), namc))))
	warning("method L-BFGS-B uses 'factr' (and 'pgtol') instead of 'reltol' and 'abstol'")
    if(npar == 1 && method == "Nelder-Mead")
        warning("one-dimensional optimization by Nelder-Mead is unreliable:\nuse \"Brent\" or optimize() directly")
    if(npar > 1 && method == "Brent")
	stop('method = "Brent" is only available for one-dimensional optimization')
    lower <- as.double(rep_len(lower, npar))
    upper <- as.double(rep_len(upper, npar))
    res <- if(method == "Brent") { ## 1-D
        if(any(!is.finite(c(upper, lower))))
           stop("'lower' and 'upper' must be finite values")
    # Karline: cannot work with fnscale
    if (con$fnscale != 1)
      stop("'Brent' method in compiled code only works with 'fnscale' = 1")
               
	res <- ccoptimize(fn, lower = lower, upper = upper, tol = con$reltol)
	names(res)[names(res) == c("minimum", "objective")] <- c("par", "value")
 #       res$value <- res$value * con$fnscale
	c(res, list(counts = c(`function` = NA, gradient = NA),
                    convergence = 0L, message = NULL))
    } else .Call("call_optim", par, fn1, gr1, method, con, lower, 
      upper, as.double(rpar), as.integer(ipar), package = "ccSolve")
    if (hessian)
        res$hessian <- .Call("call_optimhess", res$par, fn1, gr1, con, 
          as.double(rpar), as.integer(ipar), PACKAGE = "ccSolve")
    res
}

# ------------------------------------------------------------------------------

ccoptimize <- function(f, interval, ...,
		     lower=min(interval), upper=max(interval),
		     maximum=FALSE, tol=.Machine$double.eps^0.25, 
         dllname = NULL, rpar = NULL, ipar = NULL)
{
  if (is.list(f)) { 
    if (!is.null(dllname) & "dllname" %in% names(f))
      stop("If 'f' is a list that contains dllname, argument 'dllname' should be NULL")

    dllname <- f$dllname
    f <- f$func
  }

  if (! isvalid(f))
    stop("'f' should be either a compiled function or a character vector")
  
  if (class (f) == "CFunc")
    f1 <- body(f)[[2]]
  else if (is.loaded(f, PACKAGE = dllname, type = "") ||
    is.loaded(f, PACKAGE = dllname, type = "Fortran"))  {
      f1 <- getNativeSymbolInfo(f, PACKAGE = dllname)$address
  } else 
    stop(paste("'f' not loaded ", f))

  if (is.null(rpar))
    rpar <- 0.
  if (is.null(ipar))
    ipar <- 0  

  val <- .Call("cc_do_fmin", f1, as.double(lower), as.double(upper), as.double(tol), 
     as.integer(maximum), as.double(rpar), as.integer(ipar), PACKAGE = "ccSolve")

  if (maximum)
   	list(maximum = val[1], objective = val[2])
  else
  	list(minimum = val[1], objective = val[2])
}

##nice to the English (or rather the Scots)
ccoptimise <- ccoptimize
