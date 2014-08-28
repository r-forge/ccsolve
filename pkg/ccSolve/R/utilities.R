## ============================================================================
## function evaluation of compiled code objects
## ============================================================================

ccfunc <- function(fn, ...) {
  Call <- attributes (fn)$call

  if (Call == "compile.nls")
    out <- ccfunc.ccnls(fn, ...)
  else if (Call %in% c("compile.ode", "compile.steady"))
    out <- ccfunc.de(fn, ...)
  else if (Call %in% c("compile.bvp"))
    out <- ccfunc.bvp(fn, ...)
  else if (Call == "compile.multiroot")
    out <- ccfunc.multiroot(fn, ...)
  else if (Call == "compile.optim")
    out <- ccfunc.optim(fn, ...)
  else if (Call %in% c("compile.optimize", "compile.uniroot"))
    out <- ccfunc.optimize(fn, ...)
  else if (Call == "compile.dae")
    out <- ccfunc.dae(fn, ...)
  else  
    stop("No function evaluation possible for fn compiled with ", Call)

  out    
}

ccfunc.multiroot <- function(fn, start, parms = NULL, ...) {
  ccfunc.de(fn, times = 0, y = start, parms = parms, ...)$y  
}
ccfunc.bvp <- function(fn, yini, x, ...) {
  ccfunc.de(fn, times = x, y = yini, ...)$y  
}
  

getnames <- function (x) 
  if(is.matrix(x)) return(colnames(x)) else return(names(x))

isvalid <- function(f) {        # check if a function points to compiled code
  is.character(f) | class(f) == "CFunc"
}



## ============================================================================
## DLL functions
## ============================================================================

# ---------------------------------------------------------------------------
# Save a compiled code to a DLL and CFunc structure
# ---------------------------------------------------------------------------

save.cc <- function(x, file) {
  DLLname <- attributes(x)$DLL

  if (is.null(DLLname))
    stop ("'x' not of correct type: does not point to a DLL")
  if (!file.exists(DLLname))
    stop ("'x' does not point to an existing DLL")

  # correct extension of filename  (dll, so)  
  dname <- dirname(file)
  bname <- unlist(strsplit(basename(file), ".", fixed = TRUE))[1]
  extension <- unlist(strsplit(basename(DLLname), ".", fixed = TRUE))[2]
  file <- paste(bname, extension, sep = ".")

  try(dyn.unload(file), silent = TRUE)

  file.copy(from = DLLname, to = file, overwrite = TRUE)

  # accessory file with compiled code information (DLL name has changed)
  fileCF <- paste(bname, ".cc", sep = "")
  attributes(x)$DLL <- paste(getwd(), "/", file, sep = "")
  save(file = fileCF, x)
}

load.cc <- function(file) {

# open all the required files
  extension <- unlist(strsplit(basename(file), ".", fixed = TRUE))[2]

  if (extension != "cc")
    stop ("'file' should point to a compiled code object, extension '.cc'")

  if (!file.exists(file))
    stop ("'file' does not exist")

  CF <- get(load(file = file))
  attrs <- attributes(CF)
  DLLname <- attrs$DLL
    
  if (!file.exists(DLLname))
    stop ("'file' does not point to valid cc object: DLL ", DLLname, " does not exist")

#    cleanup <- function(env) {
#        unlink(DLLname)
#    }
#    reg.finalizer(environment(), cleanup, onexit = TRUE)


# load routines in DLL

  DLL <- dyn.load(DLLname)
  fn <- names(CF)
  for (i in 1:length(CF))  {
    CFi <- CF[[i]]
    code <- CF[[i]]@code
    body(CFi)[[2]] <- getNativeSymbolInfo(fn[i], DLL)$address
    CF[[i]]@.Data <- CFi
  }  
  attributes(CF) <- attrs  
  return(CF)  
}


## ============================================================================
## String functions
## ============================================================================

# ---------------------------------------------------------------------------
# Print code snippits, CFunc objects or lists of RFunc objects
# ---------------------------------------------------------------------------

printCode <- function(x) {
  if (is.list(x))
    x <- x[[1]]@code
  else if (! is.character(x))
    x <- x@code  
  cat("Program source:\n")
  lines <- strsplit(x, "\n")
  for (i in 1:length(lines[[1]])) cat(format(i, width = 3), 
    ": ", lines[[1]][i], "\n", sep = "")
}

# ---------------------------------------------------------------------------
# if FORTRAN text longer than 72 characters - split it
# ---------------------------------------------------------------------------

toFortran <- function (text) {
   if ((cl <- nchar(text)) >= 72) {
     fstring <- substr(text, 72, cl)
     text <- substr(text, 1, 71)
     while ((cf <- nchar(fstring)) > 66) {
        text <- paste(text, "\n     &", substr(fstring, 1, 66), sep = "")
        fstring <- substr(fstring, 67, cf)
     } 
     if (cf > 0)   
        text <- paste(text, "\n     &", fstring, sep = "")
        text <- paste(text, "\n")   
     } 
   return (text)
} 

# ---------------------------------------------------------------------------
# get and print dimensionality of object
# A <- list(B = 1, C= 1:2, D = diag(2))  => "B,C(2),D(2,2)"
# ---------------------------------------------------------------------------

getdims <- function (pars, language = "F95") {
  dimfun.f <- function(x) {
    nn <- ""
    if (!is.vector(x)) 
      paste(nn, "(",paste(dim(x), collapse = ","),")", sep = "") 
    else if (length(x) == 1) 
      nn 
    else 
      paste(nn,"(",length(x),")",sep="")
  }                             

  dimfun.C <- function(x) {
    nn <- ""
    if (!is.vector(x)) 
      paste(nn, "[",paste(dim(x), collapse = ","),"]", sep = "") 
    else if (length(x) == 1) 
      nn 
    else 
      paste(nn,"[",length(x),"]",sep="")
  }                             
  if (language == "C")
    return(paste(names(pars),lapply(pars, FUN = dimfun.C),sep=""))
  else
    return(paste(names(pars),lapply(pars, FUN = dimfun.f),sep=""))
}

# ---------------------------------------------------------------------------

create.ynamesc <- function (ynames) {
    npar <- length(ynames)
    head <- paste("#define ", ynames, " y[", 0:(npar-1), "]\n", sep = "", collapse = "")
    head
}

# ---------------------------------------------------------------------------

declare.ynames <- function(y, language, header, label = "y", dy = TRUE) {

  if (language == "Fortran")
    lead <- "        "
  else
    lead <- ""  

  tail <- character() 
  header2 <- character()
  if (! is.null(y))
    if (is.character(y))
      ynames <- y
    else
      ynames <- names(y)
  else
    ynames <- NULL  

  if (! is.null(ynames) & language %in% c("F95", "Fortran")){
     dynames <- paste("d", ynames, sep="")
     header <- paste(header, 
         "\n", lead, " double precision ", paste(ynames, collapse = ", "))
     if (dy) 
       header <- paste(header, 
         "\n", lead, " double precision ", paste(dynames, collapse = ", "),"\n")

     dol <- "\n"
     for (i in 1:length(ynames))
       dol <- paste(dol,    
           lead , " ", ynames[i]," = ", label, "(", i,")\n", collapse = "", sep = "")
     header2 <- dol
  
     if (dy) {
      tail <- "\n"
      for (i in 1:length(dynames))
        tail <- paste(tail,    
            lead, " f(" , i,") = ",dynames[i],"\n", collapse = "", sep = "")
     }       
  } else if (! is.null(ynames)){ # C
     dynames <- paste("d", ynames, sep="")
     header <- paste(header, "\n double ", paste(ynames, collapse = ", "),";\n")
     if (dy)
       header <- paste(header, "\n double ", paste(dynames, collapse = ", "),";\n")

     dol <- "\n"
     for (i in 1:length(ynames))
       dol <- paste(dol, 
           lead, " " , ynames[i]," = ",label,"[", i-1,"];\n", collapse = "", sep = "")
     header2 <- dol
  
     if (dy) {
      tail <- "\n"
      for (i in 1:length(dynames))
        tail <- paste(tail,    
           lead, " f[" , i-1,"] = ",dynames[i],";\n", collapse = "", sep = "")
     }
  }
  list(header = header, header2 = header2, tail = tail)
}

# ---------------------------------------------------------------------------
# Check input 
# ---------------------------------------------------------------------------

checkinput <- function (body, header, fun) {
  if (! is.character(body))
    stop ("'body' should be a character string in ", fun)
    
  if (! is.null(header))
    if (! is.character(header))
    stop ("'header' should be a character string or NULL, in ", fun)
} 
  
## ==========================================================================
ccfunc.de <- function (func, times, y, parms, dllname = NULL, initfunc = dllname, 
    rpar = NULL, ipar = NULL, nout = 0, outnames = NULL, forcings = NULL, 
    initforc = NULL, fcontrol = NULL) 
{
    cl <- attributes(func)$call
    if (! is.null(cl))
      if (!cl %in% c("compile.ode", "compile.steady", "compile.bvp", "compile.multiroot"))
        stop ("problem is not compiled with 'compile.ode' or 'compile.steady'")

    out <- DLLfunc(func, times, y, parms, dllname, initfunc, 
      rpar, ipar, nout, outnames, forcings, initforc, fcontrol)
    return(out)
}

## ==========================================================================
ccfunc.dae <- function (res, times, y, dy, parms, dllname = NULL, initfunc = dllname, 
    rpar = NULL, ipar = NULL, nout = 0, outnames = NULL, forcings = NULL, 
    initforc = NULL, fcontrol = NULL) 
{
    cl <- attributes(res)$call
    if (! is.null(cl))
      if (!cl == "compile.dae")
        stop ("problem is not compiled with 'compile.dae'")

    out <- DLLres(res, times, y, dy, parms, dllname, initfunc, 
      rpar, ipar, nout, outnames, forcings, initforc, fcontrol)
    return(out)
}
