/*
 This is an adaptation of the file "optim.c" that is part of the R-core 
 software.
 It has been adapted to work with compiled code solutions 
 bu Karline Soetaert.

 .. removed underscore in error(_(
 .. changed .External2 calls to .Call
*/

#include <R.h>
#include <Rdefines.h>
#include "ccSolve.h"


/* global variables */
C_func_type *fcall = NULL;
C_jac_type *gcall = NULL;

double* parcopy;
double* grads;
int hasgn; 
double* rpar;
int* ipar; 

SEXP getListElement(SEXP list, char *str)
{
    SEXP elmt = R_NilValue, names = getAttrib(list, R_NamesSymbol);
    int i;

    for (i = 0; i < length(list); i++)
	if (strcmp(CHAR(STRING_ELT(names, i)), str) == 0) {
	    elmt = VECTOR_ELT(list, i);
	    break;
	}
    return elmt;
}

static double * vect(int n)
{
    return (double *)R_alloc(n, sizeof(double));
}


typedef struct opt_struct
{
    SEXP R_fcall;    /* function */
    SEXP R_gcall;    /* gradient */
    SEXP R_env;      /* where to evaluate the calls */
    double* ndeps;   /* tolerances for numerical derivatives */
    double fnscale;  /* scaling for objective */
    double* parscale;/* scaling for parameters */
    int usebounds;
    double* lower, *upper;
    SEXP names;	     /* names for par */
} opt_struct, *OptStruct;

/* compiled code version ; OS->R_fcall (misnomer) points to function address */
static double fminfn_cc(int n, double *p, void *ex)
{
    int i;
    double val;

    OptStruct OS = (OptStruct) ex;

    for (i = 0; i < n; i++) {
  	  if (!R_FINITE(p[i])) error("non-finite value supplied by optim");
	    parcopy[i] = p[i] * (OS->parscale[i]);
    }

    fcall(&n, parcopy, &val, rpar, ipar);    


    val = val/(OS->fnscale);
//  	error("till here val %g %g %i", val, rpar[0], ipar[0]);	
    return val;
}

static void fmingr_cc(int n, double *p, double *df, void *ex)
{
  int i;
  double val1, val2, eps, epsused, tmp;
  OptStruct OS = (OptStruct) ex;

  if (hasgn == 1) { /* analytical derivatives */

	  for (i = 0; i < n; i++) {
	    if (!R_FINITE(p[i])) error("non-finite value supplied by optim");
	    parcopy[i] = p[i] * (OS->parscale[i]);
	  }

    gcall(&n, parcopy, grads, rpar, ipar);    

	  for (i = 0; i < n; i++)
	    df[i] = grads[i] * (OS->parscale[i])/(OS->fnscale);

  } else { /* numerical derivatives */
	  for (i = 0; i < n; i++) {
	    if (!R_FINITE(p[i])) error("non-finite value supplied by optim");
	    parcopy[i] = p[i] * (OS->parscale[i]);
	  }

  
  	if(OS->usebounds == 0) {
	    for (i = 0; i < n; i++) {
  		eps = OS->ndeps[i];

		  parcopy[i] = (p[i] + eps) * (OS->parscale[i]);
	    fcall(&n, parcopy, &val1, rpar, ipar);    
   		val1 = val1/(OS->fnscale);

		  parcopy[i] = (p[i] - eps) * (OS->parscale[i]);
	    fcall(&n, parcopy, &val2, rpar, ipar);    
   		val2 = val2/(OS->fnscale);

  		df[i] = (val1 - val2)/(2 * eps);
		  if(!R_FINITE(df[i]))
		    error(("non-finite finite-difference value [%d]"), i+1);
		  parcopy[i] = p[i] * (OS->parscale[i]);
	    }
	  } else { /* usebounds */
	    for (i = 0; i < n; i++) {
		    epsused = eps = OS->ndeps[i];
		    tmp = p[i] + eps;
		    if (tmp > OS->upper[i]) {
		      tmp = OS->upper[i];
		      epsused = tmp - p[i];
		    }
		   parcopy[i] = tmp * (OS->parscale[i]);
  	   fcall(&n, parcopy, &val1, rpar, ipar);    
   		 val1 = val1/(OS->fnscale);

  		 tmp = p[i] - eps;
		   if (tmp < OS->lower[i]) {
		    tmp = OS->lower[i];
		    eps = p[i] - tmp;
	   	 }
		   parcopy[i] = tmp * (OS->parscale[i]);
  	   fcall(&n, parcopy, &val2, rpar, ipar);    
   		 val2 = val2/(OS->fnscale);


   		 df[i] = (val1 - val2)/(epsused + eps);
		   if(!R_FINITE(df[i]))
		    error(("non-finite finite-difference value [%d]"), i+1);
		  parcopy[i] = p[i] * (OS->parscale[i]);
	    }
	  }
  }
}


/* par fn gr method options upper, lower*/                              
SEXP call_optim(SEXP par, SEXP fn, SEXP gr, SEXP method, SEXP options, 
                SEXP slower, SEXP supper, SEXP Rpar, SEXP Ipar)
{
    SEXP tmp;
    SEXP res, value, counts, conv;
    int i, np, npar=0, *mask, trace, maxit, fncount = 0, grcount = 0, nREPORT, tmax;
    int ifail = 0;
    double *dpar, *opar, val = 0.0, abstol, reltol, temp;
    const char *tn;

    OptStruct OS;

    OS = (OptStruct) R_alloc(1, sizeof(opt_struct));
    OS->usebounds = 0;

    OS->names = getAttrib(par, R_NamesSymbol);

    np = LENGTH(Ipar);
    ipar = (int *) R_alloc(np, sizeof(int));
    for (i = 0; i < np; i++)
      ipar[i] = INTEGER(Ipar)[i];

    np = LENGTH(Rpar);
    rpar = (double *) R_alloc(np, sizeof(double));
    for (i = 0; i < np; i++)
      rpar[i] = REAL(Rpar)[i];

    if (!inherits(fn, "NativeSymbol")) 
       error("'fn' is not a compiled function");
    fcall = (C_func_type *) R_ExternalPtrAddr(fn);  

    if (!isString(method)|| LENGTH(method) != 1)
  	  error("invalid '%s' argument", "method");
    tn = CHAR(STRING_ELT(method, 0));

    npar = LENGTH(par);
    dpar = vect(npar);
    opar = vect(npar);
    
    /* karline : global pars */
    parcopy = vect(npar);
    grads = vect(npar);
    
    trace = asInteger(getListElement(options, "trace"));
    OS->fnscale = asReal(getListElement(options, "fnscale"));
    tmp = getListElement(options, "parscale");
    if (LENGTH(tmp) != npar)
  	  error("'parscale' is of the wrong length");
    
    PROTECT(tmp = coerceVector(tmp, REALSXP));
    OS->parscale = vect(npar);
    for (i = 0; i < npar; i++) OS->parscale[i] = REAL(tmp)[i];
    UNPROTECT(1);
    
    for (i = 0; i < npar; i++)
	    dpar[i] = REAL(par)[i] / (OS->parscale[i]);
    PROTECT(res = allocVector(VECSXP, 5));
    SEXP names;
    PROTECT(names = allocVector(STRSXP, 5));
    SET_STRING_ELT(names, 0, mkChar("par"));
    SET_STRING_ELT(names, 1, mkChar("value"));
    SET_STRING_ELT(names, 2, mkChar("counts"));
    SET_STRING_ELT(names, 3, mkChar("convergence"));
    SET_STRING_ELT(names, 4, mkChar("message"));
    setAttrib(res, R_NamesSymbol, names);
    UNPROTECT(1);
    
    PROTECT(value = allocVector(REALSXP, 1));
    PROTECT(counts = allocVector(INTSXP, 2));
    SEXP countnames;
    PROTECT(countnames = allocVector(STRSXP, 2));
    SET_STRING_ELT(countnames, 0, mkChar("function"));
    SET_STRING_ELT(countnames, 1, mkChar("gradient"));
    setAttrib(counts, R_NamesSymbol, countnames);
    UNPROTECT(1);
    PROTECT(conv = allocVector(INTSXP, 1));
    abstol = asReal(getListElement(options, "abstol"));
    reltol = asReal(getListElement(options, "reltol"));
    maxit = asInteger(getListElement(options, "maxit"));
    if (maxit == NA_INTEGER) error("'maxit' is not an integer");

    if (strcmp(tn, "Nelder-Mead") != 0) {
    
     if (!isNull(gr)) {
       hasgn = 1; 
       if (!inherits(gr, "NativeSymbol")) 
        error("'gr' is not a compiled function");
       gcall = (C_jac_type *) R_ExternalPtrAddr(gr);  
	   } else {
       hasgn = 0;
 	   }
   }
    if (strcmp(tn, "Nelder-Mead") == 0) {
	double alpha, beta, gamm;

	alpha = asReal(getListElement(options, "alpha"));
	beta = asReal(getListElement(options, "beta"));
	gamm = asReal(getListElement(options, "gamma"));
	nmmin(npar, dpar, opar, &val, fminfn_cc, &ifail, abstol, reltol,
	      (void *)OS, alpha, beta, gamm, trace, &fncount, maxit);
	for (i = 0; i < npar; i++)
	    REAL(par)[i] = opar[i] * (OS->parscale[i]);
	grcount = NA_INTEGER;

    }
    else if (strcmp(tn, "SANN") == 0) {
	tmax = asInteger(getListElement(options, "tmax"));
	temp = asReal(getListElement(options, "temp"));
	if (trace) trace = asInteger(getListElement(options, "REPORT"));
	if (tmax == NA_INTEGER || tmax < 1) // PR#15194
	    error("'tmax' is not a positive integer");

  PROTECT(OS->R_gcall = R_NilValue); /* to generate new values... */

	samin (npar, dpar, &val, fminfn_cc, maxit, tmax, temp, trace, (void *)OS);
  UNPROTECT(1);
  for (i = 0; i < npar; i++)
	    REAL(par)[i] = dpar[i] * (OS->parscale[i]);
	fncount = npar > 0 ? maxit : 1;
	grcount = NA_INTEGER;

    } else if (strcmp(tn, "BFGS") == 0) {
	SEXP ndeps;

	nREPORT = asInteger(getListElement(options, "REPORT"));
	if (isNull(gr)) {
	    ndeps = getListElement(options, "ndeps");
	    if (LENGTH(ndeps) != npar)
		error("'ndeps' is of the wrong length");
	    OS->ndeps = vect(npar);
	    PROTECT(ndeps = coerceVector(ndeps, REALSXP));
	    for (i = 0; i < npar; i++) OS->ndeps[i] = REAL(ndeps)[i];
	    UNPROTECT(1);
	}
	mask = (int *) R_alloc(npar, sizeof(int));
	for (i = 0; i < npar; i++) mask[i] = 1;
	vmmin(npar, dpar, &val, fminfn_cc, fmingr_cc, maxit, trace, mask, abstol,
	      reltol, nREPORT, (void *)OS, &fncount, &grcount, &ifail);
	for (i = 0; i < npar; i++)
	    REAL(par)[i] = dpar[i] * (OS->parscale[i]);
    } else if (strcmp(tn, "CG") == 0) {
	int type;
	SEXP ndeps;

	type = asInteger(getListElement(options, "type"));
	if (isNull(gr)) {
	    ndeps = getListElement(options, "ndeps");
	    if (LENGTH(ndeps) != npar)
		error("'ndeps' is of the wrong length");
	    OS->ndeps = vect(npar);
	    PROTECT(ndeps = coerceVector(ndeps, REALSXP));
	    for (i = 0; i < npar; i++) OS->ndeps[i] = REAL(ndeps)[i];
	    UNPROTECT(1);
	}
	cgmin(npar, dpar, opar, &val, fminfn_cc, fmingr_cc, &ifail, abstol,
	      reltol, (void *)OS, type, trace, &fncount, &grcount, maxit);
	for (i = 0; i < npar; i++)
	    REAL(par)[i] = opar[i] * (OS->parscale[i]);

    } else if (strcmp(tn, "L-BFGS-B") == 0) {
	SEXP ndeps, smsg;
	double *lower = vect(npar), *upper = vect(npar);
	int lmm, *nbd = (int *) R_alloc(npar, sizeof(int));
	double factr, pgtol;
	char msg[60];

	nREPORT = asInteger(getListElement(options, "REPORT"));
	factr = asReal(getListElement(options, "factr"));
	pgtol = asReal(getListElement(options, "pgtol"));
	lmm = asInteger(getListElement(options, "lmm"));
	if (isNull(gr)) {
	    ndeps = getListElement(options, "ndeps");
	    if (LENGTH(ndeps) != npar)
		error("'ndeps' is of the wrong length");
	    OS->ndeps = vect(npar);
	    PROTECT(ndeps = coerceVector(ndeps, REALSXP));
	    for (i = 0; i < npar; i++) OS->ndeps[i] = REAL(ndeps)[i];
	    UNPROTECT(1);
	}

	for (i = 0; i < npar; i++) {
	    lower[i] = REAL(slower)[i] / (OS->parscale[i]);
	    upper[i] = REAL(supper)[i] / (OS->parscale[i]);
	    if (!R_FINITE(lower[i])) {
		if (!R_FINITE(upper[i])) nbd[i] = 0; else nbd[i] = 3;
	    } else {
		if (!R_FINITE(upper[i])) nbd[i] = 1; else nbd[i] = 2;
	    }
	}
	OS->usebounds = 1;
	OS->lower = lower;
	OS->upper = upper;
	lbfgsb(npar, lmm, dpar, lower, upper, nbd, &val, fminfn_cc, fmingr_cc,
	       &ifail, (void *)OS, factr, pgtol, &fncount, &grcount,
	       maxit, msg, trace, nREPORT);
	for (i = 0; i < npar; i++)
	    REAL(par)[i] = dpar[i] * (OS->parscale[i]);
	PROTECT(smsg = mkString(msg));
	SET_VECTOR_ELT(res, 4, smsg);
	UNPROTECT(1);
    } else
	error("unknown 'method'");

    if(!isNull(OS->names)) setAttrib(par, R_NamesSymbol, OS->names);
    REAL(value)[0] = val * (OS->fnscale);
    SET_VECTOR_ELT(res, 0, par); SET_VECTOR_ELT(res, 1, value);
    INTEGER(counts)[0] = fncount; INTEGER(counts)[1] = grcount;
    SET_VECTOR_ELT(res, 2, counts);
    INTEGER(conv)[0] = ifail;
    SET_VECTOR_ELT(res, 3, conv);
    UNPROTECT(4);
    return res;
}

/* par fn gr options */ 
SEXP call_optimhess(SEXP par, SEXP fn, SEXP gr, SEXP options, SEXP Rpar, SEXP Ipar)
{
    SEXP tmp, ndeps, ans;
    OptStruct OS;
    int npar, np, i , j;
    double *dpar, *df1, *df2, eps;

    OS = (OptStruct) R_alloc(1, sizeof(opt_struct));
    OS->usebounds = 0;

    npar = LENGTH(par);
    OS->names = getAttrib(par, R_NamesSymbol);

    np = LENGTH(Ipar);
    ipar = (int *) R_alloc(np, sizeof(int));
    for (i = 0; i < np; i++)
      ipar[i] = INTEGER(Ipar)[i];

    np = LENGTH(Rpar);
    rpar = (double *) R_alloc(np, sizeof(double));
    for (i = 0; i < np; i++)
      rpar[i] = REAL(Rpar)[i];

    if (!inherits(fn, "NativeSymbol")) 
       error("'fn' is not a compiled function");
    fcall = (C_func_type *) R_ExternalPtrAddr(fn);  

    OS->fnscale = asReal(getListElement(options, "fnscale"));
    tmp = getListElement(options, "parscale");
    
    if (LENGTH(tmp) != npar)
 	    error("'parscale' is of the wrong length");

    PROTECT(tmp = coerceVector(tmp, REALSXP));
    OS->parscale = vect(npar);
    for (i = 0; i < npar; i++) OS->parscale[i] = REAL(tmp)[i];
    UNPROTECT(1);

    /* karline : global pars */
    parcopy = vect(npar);
    grads = vect(npar);

    if (!isNull(gr)) {
      hasgn = 1; 
      if (!inherits(gr, "NativeSymbol")) 
        error("'gr' is not a compiled function");
      gcall = (C_jac_type *) R_ExternalPtrAddr(gr);  
	  } else {
      hasgn = 0;
 	  }

    ndeps = getListElement(options, "ndeps");
    if (LENGTH(ndeps) != npar) error("'ndeps' is of the wrong length");
    OS->ndeps = vect(npar);
    PROTECT(ndeps = coerceVector(ndeps, REALSXP));
    for (i = 0; i < npar; i++) OS->ndeps[i] = REAL(ndeps)[i];
    UNPROTECT(1);
    PROTECT(ans = allocMatrix(REALSXP, npar, npar));
    dpar = vect(npar);
    for (i = 0; i < npar; i++)
	   dpar[i] = REAL(par)[i] / (OS->parscale[i]);
    df1 = vect(npar);
    df2 = vect(npar);
    for (i = 0; i < npar; i++) {
	   eps = OS->ndeps[i]/(OS->parscale[i]);
     dpar[i] = dpar[i] + eps;
	   fmingr_cc(npar, dpar, df1, (void *)OS);
	   dpar[i] = dpar[i] - 2 * eps;
	   fmingr_cc(npar, dpar, df2, (void *)OS);
	   for (j = 0; j < npar; j++)
	    REAL(ans)[i * npar + j] = (OS->fnscale) * (df1[j] - df2[j])/
		(2 * eps * (OS->parscale[i]) * (OS->parscale[j]));
	   dpar[i] = dpar[i] + eps;
    }
    // now symmetrize
    for (i = 0; i < npar; i++) 
	   for (j = 0; j < i; j++) {
	    double tmp =
		0.5 * (REAL(ans)[i * npar + j] + REAL(ans)[j * npar + i]);
	    REAL(ans)[i * npar + j] = REAL(ans)[j * npar + i] = tmp;
	}
    SEXP nm = getAttrib(par, R_NamesSymbol);
    if(!isNull(nm)) {
	SEXP dm;
	PROTECT(dm = allocVector(VECSXP, 2));
	SET_VECTOR_ELT(dm, 0, duplicate(nm));
	SET_VECTOR_ELT(dm, 1, duplicate(nm));
	setAttrib(ans, R_DimNamesSymbol, dm);
	UNPROTECT(1);
    }
    UNPROTECT(1);
    return ans;
}
                     
                     
                     
                     
                     