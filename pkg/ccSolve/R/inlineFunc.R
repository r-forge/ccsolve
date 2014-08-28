
## ============================================================================
## inline basics - slight rewrites from package 'inline
## ============================================================================
    
setClass("CFunc",
  representation(
    code="character"
  ),
  contains="function"
)
setClass( "CFuncList", contains = "list" )

## =============================================================================
## these functions are similar, but extensions to similar functions
## with lowercase first letter from package inline
## it was not possible to use the inline functions in a way
## that is compatible with what I have in mind.
## =============================================================================

Cfunction <- function (sig = character(), body = character(), includes = character(), 
    otherdefs = character(), language = c("C++", "C", "Fortran", 
        "F95", "ObjectiveC", "ObjectiveC++"), verbose = FALSE, 
        dim = NULL, implicit = NULL, module = NULL, ## Karline: added
    convention = c(".Call", ".C", ".Fortran"), #Rcpp = FALSE, 
    cppargs = character(), cxxargs = character(), libargs = character()) 
{
    Rcpp <- FALSE   # Will not use this

# karline comment 4: to avoid using '.Call with cpp arguments if language is Fortran or F95 (or C???)
    if (missing (convention) & !missing(language))
      convention <- switch (EXPR = language, "Fortran" = ".Fortran", "F95" = ".Fortran", 
        ".C" = ".C", ObjectiveC = ".Call", "ObjectiveC++" = ".Call", "C++" = ".Call")

    convention <- match.arg(convention)
    if (missing(language)) 
        language <- ifelse(convention == ".Fortran", "Fortran", 
            "C++")
    else language <- match.arg(language)
    language <- switch(EXPR = tolower(language), cpp = "C++", 
        f = "Fortran", f95 = "F95", objc = "ObjectiveC", objcpp = , 
        `objc++` = "ObjectiveC++", language)
    f <- basename(tempfile())
    if (!is.list(sig)) {
        sig <- list(sig)
        names(sig) <- f
        names(body) <- f
    }
    if (length(sig) != length(body)) 
        stop("mismatch between the number of functions declared in 'sig' and the number of function bodies provided in 'body'")

# karline: dim should be a list, similar to sig         
    
    if (is.null(dim)) # this assumes fortran style
      dim <- as.list(rep("(*)", length(sig)))
    else {
      if (!is.list(dim)) 
        dim <- list(dim)
      if (length(dim) != length(sig)) 
        stop("mismatch between the number of functions declared in 'sig' and the number of dimensions declared in 'dim'")
    }
## end karline

#    if (Rcpp) {
#        if (!require(Rcpp)) 
#            stop("Rcpp cannot be loaded, install it or use the default Rcpp=FALSE")
#        cxxargs <- c(Rcpp:::RcppCxxFlags(), cxxargs)
#        libargs <- c(Rcpp:::RcppLdFlags(), libargs)
#    }
    if (length(cppargs) != 0) {
        args <- paste(cppargs, collapse = " ")
        if (verbose) 
            cat("Setting PKG_CPPFLAGS to", args, "\n")
        Sys.setenv(PKG_CPPFLAGS = args)
    }
    if (length(cxxargs) != 0) {
        args <- paste(cxxargs, collapse = " ")
        if (verbose) 
            cat("Setting PKG_CXXFLAGS to", args, "\n")
        Sys.setenv(PKG_CXXFLAGS = args)
    }
    if (length(libargs) != 0) {
        args <- paste(libargs, collapse = " ")
        if (verbose) 
            cat("Setting PKG_LIBS to", args, "\n")
        Sys.setenv(PKG_LIBS = args)
    }
    types <- vector(mode = "list", length = length(sig))
    for (i in seq_along(sig)) {
        if (convention == ".Call") {
            if (i == 1) {
                code <- ifelse(Rcpp, "#include <Rcpp.h>\n", paste("#include <R.h>\n#include <Rdefines.h>\n", 
                  "#include <R_ext/Error.h>\n", sep = ""))
                code <- paste(c(code, includes, ""), collapse = "\n")
                code <- paste(c(code, otherdefs, ""), collapse = "\n")
            }
            if (length(sig[[i]]) > 0) {
                funCsig <- paste("SEXP", names(sig[[i]]), collapse = ", ")
            }
            else funCsig <- ""
            funCsig <- paste("SEXP", names(sig)[i], "(", funCsig, 
                ")", sep = " ")
            if (language == "C++" || language == "ObjectiveC++") 
                code <- paste(code, "extern \"C\" {\n  ", funCsig, 
                  ";\n}\n\n", sep = "")
            code <- paste(code, funCsig, " {\n", sep = "")
            code <- paste(code, paste(body[[i]], collapse = "\n"), 
                sep = "")
            code <- paste(code, "\n  ", ifelse(Rcpp, "Rf_warning", 
                "warning"), "(\"your C program does not return anything!\");\n  return R_NilValue;\n}\n", 
                sep = "")
        }
        else if (convention == ".C") {
            if (i == 1) {
                code <- ifelse(Rcpp, "#include <Rcpp.h>\n", "#include <R.h>\n")
                code <- paste(c(code, includes, ""), collapse = "\n")
                code <- paste(c(code, otherdefs, ""), collapse = "\n")
            }
            if (length(sig[[i]]) > 0) {
                types[[i]] <- pmatch(sig[[i]], c("logical", "integer", 
                  "double", "complex", "character", "raw", "numeric"), 
                  duplicates.ok = TRUE)
                if (any(is.na(types[[i]]))) 
                  stop(paste("Unrecognized type", sig[[i]][is.na(types[[i]])]))
                decls <- c("int *", "int *", "double *", "Rcomplex *", 
                  "char **", "unsigned char *", "double *")[types[[i]]]
                funCsig <- paste(decls, names(sig[[i]]), collapse = ", ")
            }
            else funCsig <- ""
            funCsig <- paste("void", names(sig)[i], "(", funCsig, 
                ")", sep = " ")
            if (language == "C++" || language == "ObjectiveC++") 
                code <- paste(code, "extern \"C\" {\n  ", funCsig, 
                  ";\n}\n\n", sep = "")
            code <- paste(code, funCsig, " {\n", sep = "")
            code <- paste(code, paste(body[[i]], collapse = "\n"), 
                sep = "")
            code <- paste(code, "\n}\n", sep = "")
        }
        else {
            lead <- ifelse (language == "Fortran", "      ","") # Karline: first 7 columns blank only required in oldstyle fortran
            if (i == 1) {
                code <- paste(includes, collapse = "\n")
                code <- paste(c(code, otherdefs, ""), collapse = "\n")
            }                          
            if (length(sig[[i]]) > 0) {
                types[[i]] <- pmatch(sig[[i]], c("logical", "integer", 
                  "double", "complex", "character", "raw", "numeric"), 
                  duplicates.ok = TRUE)
                if (any(is.na(types[[i]]))) 
                  stop(paste("Unrecognized type", sig[[i]][is.na(types[[i]])]))
                if (6 %in% types[[i]]) 
                  stop("raw type unsupported by .Fortran()")
                decls <- c("INTEGER", "INTEGER", "DOUBLE PRECISION", 
                  "DOUBLE COMPLEX", "CHARACTER*255", "Unsupported", 
                  "DOUBLE PRECISION")[types[[i]]]
#                if (is.null(dim[[i]]))
#                  dim[[i]] <- rep("(*)", length(sig[[i]]))  
# KARLINE - changed that !!!!!!!!                  
                decls <- paste(lead, decls, " ", names(sig[[i]]), 
 #                 "(*)", sep = "", collapse = "\n")
                 dim[[i]], sep = "", collapse = "\n")               # KARLINE
                funCsig <- paste(names(sig[[i]]), collapse = ", ")
            }
            else {
                decls <- ""
                funCsig <- ""
            }
                
# Karline: if FORTRAN and text longer than 72 characters - split it
            funCsig <- paste(lead, "SUBROUTINE", names(sig)[i], 
                "(", funCsig, ")\n", sep = " ")
            if (language == "Fortran") {
              
              if ((cl <- nchar(funCsig)) >= 72) {
                fstring <- substr(funCsig, 72, cl)
                funCsig <- substr(funCsig, 1, 71)
                while ((cf <- nchar(fstring)) > 66) {
                  funCsig <- paste(funCsig, "\n     &", substr(fstring, 1, 66), sep = "")
                  fstring <- substr(fstring, 67, cf)
                } 
              if (cf > 0)   
                 funCsig <- paste(funCsig, "\n     &", fstring, sep = "")
              funCsig <- paste(funCsig, "\n")   
              } 
            } 
             
            if (is.character(module))         # implicit declaration overruled
              funCsig <- paste(funCsig, lead, "USE ", module, "\n", sep = "")    ### KARLINE            if (is.character(implicit))         # implicit declaration overruled
              funCsig <- paste(funCsig, lead, "IMPLICIT ", implicit, "\n", sep = "")    ### KARLINE ADDED
            code <- paste(code, funCsig, decls, collapse = "\n", sep = "")
            code <- paste(code, paste(body[[i]], collapse = "\n", sep = ""), ### Karline added sep=""
                sep = "")
            code <- paste(code, "\n", lead, "RETURN\n", lead, "END\n\n", 
                sep = "")
        }
    }
#    libLFile <- compileCode(f, code, language, verbose)  # changed name to CoompileCode
    libLFile <- CompileCode(f, code, language, verbose)  # to use independent from inline
    cleanup <- function(env) {
        if (f %in% names(getLoadedDLLs())) 
            dyn.unload(libLFile)
        unlink(libLFile)
    }
    reg.finalizer(environment(), cleanup, onexit = TRUE)
    res <- vector("list", length(sig))
    names(res) <- names(sig)
    for (i in seq_along(sig)) {
        res[[i]] <- new("CFunc", code = code)
        fn <- function(arg) {
            NULL
        }
        DLL <- dyn.load(libLFile)
        args <- formals(fn)[rep(1, length(sig[[i]]))]
        names(args) <- names(sig[[i]])
        formals(fn) <- args
        if (convention == ".Call") {
            body <- quote(CONVENTION("EXTERNALNAME", ARG))[c(1:2, 
                rep(3, length(sig[[i]])))]
            for (j in seq(along = sig[[i]])) body[[j + 2]] <- as.name(names(sig[[i]])[j])
        }
        else {
            body <- quote(CONVENTION("EXTERNALNAME", as.logical(ARG), 
                as.integer(ARG), as.double(ARG), as.complex(ARG), 
                as.character(ARG), as.raw(ARG), as.double(ARG)))[c(1:2, 
                types[[i]] + 2)]
            names(body) <- c(NA, "", names(sig[[i]]))
            for (j in seq(along = sig[[i]])) body[[j + 2]][[2]] <- as.name(names(sig[[i]])[j])
        }
        body[[1]] <- get(convention)
        body[[2]] <- getNativeSymbolInfo(names(sig)[i], DLL)$address
        body(fn) <- body
        res[[i]]@.Data <- fn
    }
    if (verbose) {
        cat("Program source:\n")
        lines <- strsplit(code, "\n")
        for (i in 1:length(lines[[1]])) cat(format(i, width = 3), 
            ": ", lines[[1]][i], "\n", sep = "")
    }
    remove(list = c("args", "body", "fn", "funCsig", "i", "includes", 
        "j"))
        
    
    if (length(res) == 1 && names(res) == f) 
        R <- res[[1]]
        
    else R <- new("CFuncList", res)
    attr(R, "DLL") <- libLFile  
    return(R)
}


CompileCode <- function (f, code, language, verbose) 
{
    wd = getwd()
    on.exit(setwd(wd))
    if (.Platform$OS.type == "windows") {
        dir <- gsub("\\\\", "/", tempdir())
        libCFile <- paste(dir, "/", f, ".EXT", sep = "")
        libLFile <- paste(dir, "/", f, ".dll", sep = "")
        libLFile2 <- paste(dir, "/", f, ".dll", sep = "")
    }
    else {
        libCFile <- paste(tempdir(), "/", f, ".EXT", sep = "")
        libLFile <- paste(tempdir(), "/", f, .Platform$dynlib.ext, 
            sep = "")
        libLFile2 <- paste(tempdir(), "/", f, ".sl", sep = "")
    }
    extension <- switch(language, `C++` = ".cpp", C = ".c", Fortran = ".f", 
        F95 = ".f95", ObjectiveC = ".m", `ObjectiveC++` = ".mm")
    libCFile <- sub(".EXT$", extension, libCFile)
    write(code, libCFile)
    if (file.exists(libLFile)) 
        file.remove(libLFile)
    if (file.exists(libLFile2)) 
        file.remove(libLFile2)
    setwd(dirname(libCFile))
    errfile <- paste(basename(libCFile), ".err.txt", sep = "")
    cmd <- paste(R.home(component = "bin"), "/R CMD SHLIB ", 
        basename(libCFile), " 2> ", errfile, sep = "")
    if (verbose) 
        cat("Compilation argument:\n", cmd, "\n")
    compiled <- system(cmd, intern = !verbose)
    errmsg <- readLines(errfile)
    unlink(errfile)
    writeLines(errmsg)
    setwd(wd)
    if (!file.exists(libLFile) && file.exists(libLFile2)) 
        libLFile <- libLFile2
    if (!file.exists(libLFile)) {
        cat("\nERROR(s) during compilation: source code errors or compiler configuration errors!\n")
        cat("\nProgram source:\n")
        code <- strsplit(code, "\n")
        for (i in 1:length(code[[1]])) cat(format(i, width = 3), 
            ": ", code[[1]][i], "\n", sep = "")
        stop(paste("Compilation ERROR, function(s)/method(s) not created!", 
            paste(errmsg, collapse = "\n")))
    }

    return(libLFile)
}
