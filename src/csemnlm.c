/* ********************
 * Please note that this function "csenlm" is modified from the R core program main/optimize.c
 * by Zhenghua Nie.
 * The major goal of this modification is that we want to call "nlm" in C when
 * the objetive function and gradients, or hessians are computed in C. This is for 
 * speed-up of the package "sem". In future,  if we find more efficient solver of 
 * non-linear optimization problems,  we may replace the solver.
 * ********************/
/*
 *  R : A Computer Language for Statistical Data Analysis
 *  Copyright (C) 1995, 1996  Robert Gentleman and Ross Ihaka
 *  Copyright (C) 1998--2011  The R Development Core Team
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, a copy is available at
 *  http://www.r-project.org/Licenses/
 */
#include "csemnlm.h"

/* General Nonlinear Optimization */
/* Initialize the storage in the table of computed function values */
static int csem_isnan(double x)
{
		return(x != x);
}

static void FT_init(int n, int FT_size, function_info *state)
{
    int i, j;
    int have_gradient, have_hessian;
    ftable *Ftable;

		int modeln = state->model->n;
		int m = state->model->m;

    have_gradient = state->have_gradient;
    have_hessian = state->have_hessian;

    Ftable = (ftable *)R_alloc(FT_size, sizeof(ftable));

		for (i = 0; i < FT_size; i++) {
				Ftable[i].x = (double *)R_alloc(n, sizeof(double));
				Ftable[i].A = (double *)R_alloc(m*m, sizeof(double)); //A: m-by-m 
				Ftable[i].P = (double *)R_alloc(m*m, sizeof(double)); //P: m-by-m
				Ftable[i].C = (double *)R_alloc(modeln*modeln, sizeof(double)); //After the compuation,  C is modeln-by-modeln.

				/* initialize to unlikely parameter values */
				for (j = 0; j < n; j++) {
						Ftable[i].x[j] = DBL_MAX;
				}
				if (have_gradient) {
						Ftable[i].grad = (double *)R_alloc(n, sizeof(double));
						if (have_hessian) {
								Ftable[i].hess = (double *)R_alloc(n * n, sizeof(double));
						}
				}
		}
		state->Ftable = Ftable;
		state->FT_size = FT_size;
		state->FT_last = -1;
}

static void msem_FT_init(int n, int FT_size, msem_function_info *state)
{
    int i, j;
    int have_gradient, have_hessian;
    msem_ftable *Ftable;
		int G = state->model->G;   //number of groups

		int *modeln = (int *)R_alloc(G, sizeof(int));
		int *m = (int *)R_alloc(G, sizeof(int));

		Memcpy(modeln, INTEGER(AS_INTEGER(state->model->n)), G);
		Memcpy(m, INTEGER(AS_INTEGER(state->model->m)), G);

		int totalm=0, totaln=0;
		for(i = 0 ; i < G; ++i) {
				totalm += m[i]*m[i];
				totaln += modeln[i]*modeln[i];
		}
		state->sizeC = totaln;
		state->sizeAP = totalm;

    have_gradient = state->have_gradient;
    have_hessian = state->have_hessian;

    Ftable = (msem_ftable *)R_alloc(FT_size, sizeof(msem_ftable));

		for (i = 0; i < FT_size; i++) {
				Ftable[i].x = (double *)R_alloc(n, sizeof(double));
				Ftable[i].A = (double *)R_alloc(totalm, sizeof(double)); //A: m-by-m 
				Ftable[i].P = (double *)R_alloc(totalm, sizeof(double)); //P: m-by-m
				Ftable[i].C = (double *)R_alloc(totaln, sizeof(double)); //After the compuation,  C is modeln-by-modeln.
				Ftable[i].ff = (double *)R_alloc(G, sizeof(double));

				/* initialize to unlikely parameter values */
				for (j = 0; j < n; j++) {
						Ftable[i].x[j] = DBL_MAX;
				}
				if (have_gradient) {
						Ftable[i].grad = (double *)R_alloc(n, sizeof(double));
						if (have_hessian) {
								Ftable[i].hess = (double *)R_alloc(n * n, sizeof(double));
						}
				}
		}
		state->Ftable = Ftable;
		state->FT_size = FT_size;
		state->FT_last = -1;

		return;
}
/* Store an entry in the table of computed function values */

static void FT_store(int n, const double f, const double *x, const double *grad,
		     const double *hess, const double *A, const double *P, const double *C, function_info *state)
{
    int ind;

		ind = (++(state->FT_last)) % (state->FT_size);
		state->Ftable[ind].fval = f;
		Memcpy(state->Ftable[ind].x, x, n);
		Memcpy(state->Ftable[ind].C, C, state->model->n*state->model->n);
		Memcpy(state->Ftable[ind].A, A, state->model->m*state->model->m);
		Memcpy(state->Ftable[ind].P, P, state->model->m*state->model->m);
		if (grad) {
				Memcpy(state->Ftable[ind].grad, grad, n);
				if (hess) {
						Memcpy(state->Ftable[ind].hess, hess, n * n);
				}
		}
}

static void msem_FT_store(int n, const double f, const double *x, const double *grad,
		     const double *hess, const double *A, const double *P, const double *C, const double *ff, msem_function_info *state)
{
    int ind;

		ind = (++(state->FT_last)) % (state->FT_size);
		state->Ftable[ind].fval = f;
		Memcpy(state->Ftable[ind].x, x, n);
		Memcpy(state->Ftable[ind].C, C, state->sizeC);
		Memcpy(state->Ftable[ind].A, A, state->sizeAP);
		Memcpy(state->Ftable[ind].P, P, state->sizeAP);
		Memcpy(state->Ftable[ind].ff, ff, state->model->G);
		if (grad) {
				Memcpy(state->Ftable[ind].grad, grad, n);
				if (hess) {
						Memcpy(state->Ftable[ind].hess, hess, n * n);
				}
		}
}

/* Check for stored values in the table of computed function values.
   Returns the index in the table or -1 for failure */

static int FT_lookup(int n, const double *x, function_info *state)
{
    double *ftx;
    int i, j, ind, matched;
    int FT_size, FT_last;
    ftable *Ftable;

    FT_last = state->FT_last;
    FT_size = state->FT_size;
    Ftable = state->Ftable;

		for (i = 0; i < FT_size; i++) {
				ind = (FT_last - i) % FT_size;
				/* why can't they define modulus correctly */
				if (ind < 0) ind += FT_size;
				ftx = Ftable[ind].x;
				if (ftx) {
						matched = 1;
						for (j = 0; j < n; j++) {
								if (x[j] != ftx[j]) {
										matched = 0;
										break;
								}
						}
						if (matched) return ind;
				}
		}
		return -1;
}

static int msem_FT_lookup(int n, const double *x, msem_function_info *state)
{
    double *ftx;
    int i, j, ind, matched;
    int FT_size, FT_last;
    msem_ftable *Ftable;

    FT_last = state->FT_last;
    FT_size = state->FT_size;
    Ftable = state->Ftable;

		for (i = 0; i < FT_size; i++) {
				ind = (FT_last - i) % FT_size;
				/* why can't they define modulus correctly */
				if (ind < 0) ind += FT_size;
				ftx = Ftable[ind].x;
				if (ftx) {
						matched = 1;
						for (j = 0; j < n; j++) {
								if (x[j] != ftx[j]) {
										matched = 0;
										break;
								}
						}
						if (matched) return ind;
				}
		}
		return -1;
}

/* This how the optimizer sees them */

static void fcn(int n, const double x[], double *f, function_info *state)
{

		ftable *Ftable;
		double *g=NULL;
		double *h=NULL;
		double *C=NULL;
		double *A=NULL;
		double *P=NULL;
		int i;

		Ftable = state->Ftable;
		if ((i = FT_lookup(n, x, state)) >= 0) {
				*f = Ftable[i].fval;
				return;
		}

		for (i = 0; i < n; i++) {
				if (!R_FINITE(x[i])) error(("non-finite value supplied by 'nlm'"));
		}


		if(state->have_gradient) {
				g = (double *)R_alloc(n, sizeof(double));
				memset(g, 0, n*sizeof(double));
				if(state->have_hessian) {
						h = (double *)R_alloc(n*n, sizeof(double));
						memset(h, 0, n*n*sizeof(double));
				}
		}
	  int m = state->model->m;
		int modeln = state->model->n;
		int maxmn = (m > modeln ? m: modeln);
		C = (double *)R_alloc(maxmn*maxmn, sizeof(double)); //After the compuation,  C is n-by-n.
		A = (double *)R_alloc(m*m, sizeof(double));
		P = (double *)R_alloc(m*m, sizeof(double));

		myfcn_p myobjfun = (myfcn_p)state->myobjfun;

		(myobjfun)(n, x, f, g, h,  A, P, C, state);
		++state->n_eval;  //number of the evaluations.

		if((*f != *f) || !R_FINITE(*f)) {
				warning(("NA//Inf replaced by maximum positive value"));
				*f = DBL_MAX;
		}

    FT_store(n, *f, x, g, h, A, P, C, state);
		return;
}

static void msem_fcn(int n, const double x[], double *f, msem_function_info *state)
{

		msem_ftable *Ftable;
		double *g=NULL;
		double *h=NULL;
		double *C=NULL;
		double *A=NULL;
		double *P=NULL;
		double *ff=NULL;
		int i;

		Ftable = state->Ftable;
		if ((i = msem_FT_lookup(n, x, state)) >= 0) {
				*f = Ftable[i].fval;
				return;
		}

		for (i = 0; i < n; i++) {
				if (!R_FINITE(x[i])) error(("non-finite value supplied by 'nlm'"));
		}


		if(state->have_gradient) {
				g = (double *)R_alloc(n, sizeof(double));
				memset(g, 0, n*sizeof(double));
				if(state->have_hessian) {
						h = (double *)R_alloc(n*n, sizeof(double));
						memset(h, 0, n*n*sizeof(double));
				}
		}

		C = (double *)R_alloc(state->sizeC, sizeof(double)); //After the compuation,  C is n-by-n. I need to check the size.
		A = (double *)R_alloc(state->sizeAP, sizeof(double));
		P = (double *)R_alloc(state->sizeAP, sizeof(double));
		ff = (double *)R_alloc(state->model->G, sizeof(double));

		msem_fcn_p myobjfun = (msem_fcn_p)state->myobjfun;

		(myobjfun)(n, x, f, g, h,  A, P, C,ff,  state);
		++state->n_eval;  //number of the evaluations.

		if((*f != *f) || !R_FINITE(*f)) {
				warning(("NA/Inf replaced by maximum positive value"));
				*f = DBL_MAX;
		}

    msem_FT_store(n, *f, x, g, h, A, P, C, ff,  state);
		return;

}

/* gradient */
static void Cd1fcn(int n, const double x[], double *g, function_info *state)
{
    int ind;

		if ((ind = FT_lookup(n, x, state)) < 0) {	/* shouldn't happen */
				fcn(n, x, g, state);
				if ((ind = FT_lookup(n, x, state)) < 0) {
						error(("function value caching for optimization is seriously confused"));
				}
		}

		Memcpy(g, state->Ftable[ind].grad, n);

		return;
}

static void msem_Cd1fcn(int n, const double x[], double *g, msem_function_info *state)
{
    int ind;

		if ((ind = msem_FT_lookup(n, x, state)) < 0) {	/* shouldn't happen */
				msem_fcn(n, x, g, state);
				if ((ind = msem_FT_lookup(n, x, state)) < 0) {
						error(("function value caching for optimization is seriously confused"));
				}
		}

		Memcpy(g, state->Ftable[ind].grad, n);

		return;
}

/* hessian */
static void Cd2fcn(int nr, int n, const double x[], double *h, function_info *state)
{
		int j, ind;

		if ((ind = FT_lookup(n, x, state)) < 0) {	/* shouldn't happen */
				fcn(n, x, h, state);
				if ((ind = FT_lookup(n, x, state)) < 0) {
						error(("function value caching for optimization is seriously confused"));
				}
		}
		for (j = 0; j < n; j++) {  /* fill in lower triangle only */
				Memcpy( h + j*(n + 1), state->Ftable[ind].hess + j*(n + 1), n - j);
		}
		return;
}

static void msem_Cd2fcn(int nr, int n, const double x[], double *h, function_info *state)
{
		int j, ind;

		if ((ind = FT_lookup(n, x, state)) < 0) {	/* shouldn't happen */
				fcn(n, x, h, state);
				if ((ind = FT_lookup(n, x, state)) < 0) {
						error(("function value caching for optimization is seriously confused"));
				}
		}
		for (j = 0; j < n; j++) {  /* fill in lower triangle only */
				Memcpy( h + j*(n + 1), state->Ftable[ind].hess + j*(n + 1), n - j);
		}
		return;
}

/* A, P, C */
static void returnAPCfcn(int n, const double x[], double *A, double *P, double *C, function_info *state)
{
		int  ind;

		if ((ind = FT_lookup(n, x, state)) < 0) {	/* shouldn't happen */
				fcn(n, x, C, state);
				if ((ind = FT_lookup(n, x, state)) < 0) {
						error(("function value caching for optimization is seriously confused"));
				}
		}
		const int modeln = state->model->n;
		const int m = state->model->m;

		Memcpy(A, state->Ftable[ind].A, m*m);
		Memcpy(P, state->Ftable[ind].P, m*m);
		Memcpy(C, state->Ftable[ind].C, modeln*modeln);

		return;
}

static void msem_returnAPCfcn(int n, const double x[], double *A, double *P, double *C, double *ff, msem_function_info *state)
{
		int  ind;

		if ((ind = msem_FT_lookup(n, x, state)) < 0) {	/* shouldn't happen */
				msem_fcn(n, x, C, state);
				if ((ind = msem_FT_lookup(n, x, state)) < 0) {
						error(("function value caching for optimization is seriously confused"));
				}
		}

		Memcpy(A, state->Ftable[ind].A, state->sizeAP);
		Memcpy(P, state->Ftable[ind].P, state->sizeAP);
		Memcpy(C, state->Ftable[ind].C, state->sizeC);
		Memcpy(ff, state->Ftable[ind].ff, state->model->G);

		return;
}

	/* Fatal errors - we don't deliver an answer */

static void opterror(int nerr)
{
    switch(nerr) {
    case -1:
	error(("non-positive number of parameters in nlm"));
    case -2:
	error(("nlm is inefficient for 1-d problems"));
    case -3:
	error(("invalid gradient tolerance in nlm"));
    case -4:
	error(("invalid iteration limit in nlm"));
    case -5:
	error(("minimization function has no good digits in nlm"));
    case -6:
	error(("no analytic gradient to check in nlm!"));
    case -7:
	error(("no analytic Hessian to check in nlm!"));
    case -21:
	error(("probable coding error in analytic gradient"));
    case -22:
	error(("probable coding error in analytic Hessian"));
    default:
	error(("*** unknown error message (msg = %d) in nlm()\n*** should not happen!"), nerr);
    }
}


	/* Warnings - we return a value, but print a warning */

static void optcode(int code)
{
    switch(code) {
    case 1:
	Rprintf(("Relative gradient close to zero.\n"));
	Rprintf(("Current iterate is probably solution.\n"));
	break;
    case 2:
	Rprintf(("Successive iterates within tolerance.\n"));
	Rprintf(("Current iterate is probably solution.\n"));
	break;
    case 3:
	Rprintf(("Last global step failed to locate a point lower than x.\n"));
	Rprintf(("Either x is an approximate local minimum of the function,\n\
the function is too non-linear for this algorithm,\n\
or steptol is too large.\n"));
	break;
    case 4:
	Rprintf(("Iteration limit exceeded.  Algorithm failed.\n"));
	break;
    case 5:
	Rprintf(("Maximum step size exceeded 5 consecutive times.\n\
Either the function is unbounded below,\n\
becomes asymptotic to a finite value\n\
from above in some direction,\n"\
"or stepmx is too small.\n"));
	break;
    }
    Rprintf("\n");
}

/* NOTE: The actual Dennis-Schnabel algorithm `optif9' is in ../appl/uncmin.c */


    /* `x0' : inital parameter value */
    /* `want_hessian' : H. required? */
    /* `typsize' : typical size of parameter elements */
    /* `fscale' : expected function size */
    /* `msg' (bit pattern) */
    /* `iterlim' (def. 100) */
    //iagflg = 0;			/* No analytic gradient */
   // iahflg = 0;			/* No analytic hessian */
    // n = 0;
    //x = fixparam(x0, &n);

    //typsiz = fixparam(typsize, &n);
SEXP csemnlm(double *x0, int n, int iagflg,  int iahflg, int want_hessian, 
				double *typsiz, double fscale, int msg, int ndigit, double gradtl, 
				double stepmx, double steptol, int itnlim, model_info *model, myfcn_p myobjfun, 
				int optimize)
{
    SEXP value, names;

		if(SEM_DEBUG) Rprintf("Optimize: [%d]\n", optimize);

    double *x,*xpls, *gpls, fpls, *a, *wrk, dlt;

    int code, i, j, k, method, iexp, omsg,  itncnt;

		x = (double *)R_alloc(n, sizeof(double));
		Memcpy(x, x0, n);

		//initial function_info,  this will be transfered into nlm.
    function_info *state;

    state = (function_info *) R_alloc(1, sizeof(function_info));
		state->model = model;
		state->have_gradient = iagflg;
		state->have_hessian = iahflg;
		state->n_eval = 0;
		state->myobjfun = (myfcn_p *) myobjfun;


/* .Internal(
 *	nlm(function(x) f(x, ...), p, hessian, typsize, fscale,
 *	    msg, ndigit, gradtol, stepmax, steptol, iterlim)
 */

    omsg = msg ;

    if (((msg/4) % 2) && !iahflg) { /* skip check of analytic Hessian */
      msg -= 4;
    }
    if (((msg/2) % 2) && !iagflg) { /* skip check of analytic gradient */
      msg -= 2;
    }

		FT_init(n, FT_SIZE, state);


    method = 1;	/* Line Search */
    iexp = iahflg ? 0 : 1; /* Function calls are expensive */
    dlt = 1.0;

    xpls = (double*)R_alloc(n, sizeof(double));
    gpls = (double*)R_alloc(n, sizeof(double));
    a = (double*)R_alloc(n*n, sizeof(double));
    wrk = (double*)R_alloc(8*n, sizeof(double));


		//we will not optimize,  only return the objective value,  gradients and Hessian.
		if(optimize !=1 ) {

				int m = state->model->m;
				int modeln = state->model->n;
				int maxmn = (m > modeln ? m: modeln);
				double *matrixA = (double *)R_alloc(m*m, sizeof(double)); //After the compuation,  C is n-by-n.
				double *P = (double *)R_alloc(m*m, sizeof(double)); //After the compuation,  C is n-by-n.
				double *C = (double *)R_alloc(maxmn*maxmn, sizeof(double)); //After the compuation,  C is n-by-n.

				memset(gpls, 0, n*sizeof(double));
				memset(a, 0, n*n*sizeof(double));

				(*myobjfun)(n, x0, &fpls, gpls, a , matrixA, P,  C,  state);

				int num_objs=2;  //x0, objective,
				//A, P, C
				if(!csem_isnan(*matrixA)) ++num_objs;
				if(!csem_isnan(*P)) ++num_objs;
				if(!csem_isnan(*C)) ++num_objs;

				if(iagflg) {
						++num_objs;
						if(!iahflg && want_hessian) {//we need to compute Hessian if it is not provided.)
								fdhess(n, x0, fpls, (fcn_p)fcn, state, a,  n,  &wrk[0], &wrk[n], ndigit, typsiz );
						for (i = 0; i < n; i++)
								for (j = 0; j < i; j++)
										a[i + j * n] = a[j + i * n];
						}
				} else if(want_hessian) {
						fdhess(n, x0, fpls, (fcn_p)fcn, state, a,  n,  &wrk[0], &wrk[n], ndigit, typsiz );
						for (i = 0; i < n; i++)
								for (j = 0; j < i; j++)
										a[i + j * n] = a[j + i * n];
				}

				if(want_hessian) ++num_objs;

				PROTECT(value = allocVector(VECSXP, num_objs));
				PROTECT(names = allocVector(STRSXP, num_objs));
				k = 0;

				SET_STRING_ELT(names, k, mkChar("minimum"));
				SET_VECTOR_ELT(value, k, ScalarReal(fpls));
				k++;

				SET_STRING_ELT(names, k, mkChar("estimate"));
				SET_VECTOR_ELT(value, k, allocVector(REALSXP, n));
				for (i = 0; i < n; i++)
						REAL(VECTOR_ELT(value, k))[i] = x0[i];
				k++;

				if(iagflg) {
						SET_STRING_ELT(names, k, mkChar("gradient"));
						SET_VECTOR_ELT(value, k, allocVector(REALSXP, n));
						for (i = 0; i < n; i++)
								REAL(VECTOR_ELT(value, k))[i] = gpls[i];
						k++;
				}

				if(want_hessian){
						SET_STRING_ELT(names, k, mkChar("hessian"));
						SET_VECTOR_ELT(value, k, allocMatrix(REALSXP, n, n));
						for (i = 0; i < n * n; i++)
								REAL(VECTOR_ELT(value, k))[i] = a[i];
						k++;
				}

				/* A */
				if(!csem_isnan(*matrixA)) {
						SET_STRING_ELT(names, k, mkChar("A"));
						SET_VECTOR_ELT(value, k, allocMatrix(REALSXP, m, m));
						for (i = 0; i < m * m; i++)
								REAL(VECTOR_ELT(value, k))[i] = matrixA[i];
						k++;
				}

				/* P */
				if(!csem_isnan(*P)) {
						SET_STRING_ELT(names, k, mkChar("P"));
						SET_VECTOR_ELT(value, k, allocMatrix(REALSXP, m, m));
						for (i = 0; i < m * m; i++)
								REAL(VECTOR_ELT(value, k))[i] = P[i];
						k++;
				}

				/* C */
				if(!csem_isnan(*C)) { 
						SET_STRING_ELT(names, k, mkChar("C"));
						SET_VECTOR_ELT(value, k, allocMatrix(REALSXP, modeln, modeln));
						for (i = 0; i < modeln * modeln; i++)
								REAL(VECTOR_ELT(value, k))[i] = C[i];
						k++;
				}

				setAttrib(value, R_NamesSymbol, names);
				UNPROTECT(3);
		}
		else {


				/*
				 *	 Dennis + Schnabel Minimizer
				 *
				 *	  SUBROUTINE OPTIF9(NR,N,X,FCN,D1FCN,D2FCN,TYPSIZ,FSCALE,
				 *	 +	   METHOD,IEXP,MSG,NDIGIT,ITNLIM,IAGFLG,IAHFLG,IPR,
				 *	 +	   DLT,GRADTL,STEPMX,STEPTOL,
				 *	 +	   XPLS,FPLS,GPLS,ITRMCD,A,WRK)
				 *
				 *
				 *	 Note: I have figured out what msg does.
				 *	 It is actually a sum of bit flags as follows
				 *	   1 = don't check/warn for 1-d problems
				 *	   2 = don't check analytic gradients
				 *	   4 = don't check analytic hessians
				 *	   8 = don't print start and end info
				 *	  16 = print at every iteration
				 *	 Using msg=9 is absolutely minimal
				 *	 I think we always check gradients and hessians
				 */


				optif9(n, n, x, (fcn_p) fcn, (fcn_p) Cd1fcn, (d2fcn_p) Cd2fcn,
								state, typsiz, fscale, method, iexp, &msg, ndigit, itnlim,
								iagflg, iahflg, dlt, gradtl, stepmx, steptol, xpls, &fpls,
								gpls, &code, a, wrk, &itncnt);

				if(SEM_DEBUG) {
						if(state->have_hessian) 
								Rprintf("Hessian is provided.\n");
						else 
								Rprintf("Hessian is not provided.\n");

						Rprintf("The number of function evaluations: [%d]\n", state->n_eval);

				}

				if (msg < 0)
						opterror(msg);
				if (code != 0 && (omsg&8) == 0)
						optcode(code);

				int num_objs = 5;

				if (want_hessian) {
						++num_objs;
						fdhess(n, xpls, fpls, (fcn_p) fcn, state, a, n, &wrk[0], &wrk[n],
										ndigit, typsiz);
						for (i = 0; i < n; i++)
								for (j = 0; j < i; j++)
										a[i + j * n] = a[j + i * n];
				}

				//A, P, C
				int modeln = state->model->n;
				int m = state->model->m;
				double *matrixA = (double *)R_alloc(m*m, sizeof(double)); //After the compuation,  C is n-by-n.
				double *P = (double *)R_alloc(m*m, sizeof(double)); //After the compuation,  C is n-by-n.
				double *C = (double *)R_alloc(modeln*modeln, sizeof(double)); //After the compuation,  C is n-by-n.
				returnAPCfcn(n, xpls, matrixA, P, C, state);
				if(!csem_isnan(*matrixA)) ++num_objs;
				if(!csem_isnan(*P)) ++num_objs;
				if(!csem_isnan(*C)) ++num_objs;

				PROTECT(value = allocVector(VECSXP, num_objs));
				PROTECT(names = allocVector(STRSXP, num_objs));
				k = 0;

				SET_STRING_ELT(names, k, mkChar("minimum"));
				SET_VECTOR_ELT(value, k, ScalarReal(fpls));
				k++;

				SET_STRING_ELT(names, k, mkChar("estimate"));
				SET_VECTOR_ELT(value, k, allocVector(REALSXP, n));
				for (i = 0; i < n; i++)
						REAL(VECTOR_ELT(value, k))[i] = xpls[i];
				k++;

				SET_STRING_ELT(names, k, mkChar("gradient"));
				SET_VECTOR_ELT(value, k, allocVector(REALSXP, n));
				for (i = 0; i < n; i++)
						REAL(VECTOR_ELT(value, k))[i] = gpls[i];
				k++;

				if (want_hessian) {
						SET_STRING_ELT(names, k, mkChar("hessian"));
						SET_VECTOR_ELT(value, k, allocMatrix(REALSXP, n, n));
						for (i = 0; i < n * n; i++)
								REAL(VECTOR_ELT(value, k))[i] = a[i];
						k++;
				}

				SET_STRING_ELT(names, k, mkChar("code"));
				SET_VECTOR_ELT(value, k, allocVector(INTSXP, 1));
				INTEGER(VECTOR_ELT(value, k))[0] = code;
				k++;

				SET_STRING_ELT(names, k, mkChar("iterations"));
				SET_VECTOR_ELT(value, k, allocVector(INTSXP, 1));
				INTEGER(VECTOR_ELT(value, k))[0] = itncnt;
				k++;

				/* A */
				if(!csem_isnan(*matrixA)) {
						SET_STRING_ELT(names, k, mkChar("A"));
						SET_VECTOR_ELT(value, k, allocMatrix(REALSXP, m, m));
						for (i = 0; i < m * m; i++)
								REAL(VECTOR_ELT(value, k))[i] = matrixA[i];
						k++;
				}

				/* P */
				if(!csem_isnan(*P)) {
						SET_STRING_ELT(names, k, mkChar("P"));
						SET_VECTOR_ELT(value, k, allocMatrix(REALSXP, m, m));
						for (i = 0; i < m * m; i++)
								REAL(VECTOR_ELT(value, k))[i] = P[i];
						k++;
				}

				/* C */
				if(!csem_isnan(*C)) {
						SET_STRING_ELT(names, k, mkChar("C"));
						SET_VECTOR_ELT(value, k, allocMatrix(REALSXP, modeln, modeln));
						for (i = 0; i < modeln * modeln; i++)
								REAL(VECTOR_ELT(value, k))[i] = C[i];
						k++;
				}

				setAttrib(value, R_NamesSymbol, names);
				UNPROTECT(3);
		}

    return value;
}

SEXP cmsemnlm(double *x0, int n, int iagflg,  int iahflg, int want_hessian, 
				double *typsiz, double fscale, int msg, int ndigit, double gradtl, 
				double stepmx, double steptol, int itnlim, msem_model_info *model, msem_fcn_p myobjfun, 
				int optimize)
{
    SEXP value, names;

		if(SEM_DEBUG) Rprintf("Optimize: [%d]\n", optimize);

    double *x,*xpls, *gpls, fpls, *a, *wrk, dlt;

    int code, i, j, k, method, iexp, omsg,  itncnt;

		x = (double *)R_alloc(n, sizeof(double));
		Memcpy(x, x0, n);

		//initial function_info,  this will be transfered into nlm.
    msem_function_info *state;

    state = (msem_function_info *) R_alloc(1, sizeof(msem_function_info));
		state->model = model;
		state->have_gradient = iagflg;
		state->have_hessian = iahflg;
		state->n_eval = 0;
		state->myobjfun = (msem_fcn_p *) myobjfun;

/* .Internal(
 *	nlm(function(x) f(x, ...), p, hessian, typsize, fscale,
 *	    msg, ndigit, gradtol, stepmax, steptol, iterlim)
 */

    omsg = msg ;

    if (((msg/4) % 2) && !iahflg) { /* skip check of analytic Hessian */
      msg -= 4;
    }
    if (((msg/2) % 2) && !iagflg) { /* skip check of analytic gradient */
      msg -= 2;
    }

		msem_FT_init(n, 2, state);  //FT_SIZE


    method = 1;	/* Line Search */
    iexp = iahflg ? 0 : 1; /* Function calls are expensive */
    dlt = 1.0;

    xpls = (double*)R_alloc(n, sizeof(double));
    gpls = (double*)R_alloc(n, sizeof(double));
    a = (double*)R_alloc(n*n, sizeof(double));
    wrk = (double*)R_alloc(8*n, sizeof(double));


		//we will not optimize,  only return the objective value,  gradients and Hessian.
		if(optimize !=1 ) {

				int sizeAP = state->sizeAP;
				int sizeC = state->sizeC;
				int maxmn = (sizeAP > sizeC ? sizeAP: sizeC);
				double *matrixA = (double *)R_alloc(sizeAP, sizeof(double)); //After the compuation,  C is n-by-n.
				double *P = (double *)R_alloc(sizeAP, sizeof(double)); //After the compuation,  C is n-by-n.
				double *C = (double *)R_alloc(maxmn, sizeof(double)); //After the compuation,  C is n-by-n.
				double *ff = (double *)R_alloc(state->model->G, sizeof(double));

				memset(gpls, 0, n*sizeof(double));
				memset(a, 0, n*n*sizeof(double));

				(*myobjfun)(n, x0, &fpls, gpls, a , matrixA, P,  C, ff,  state);

				int num_objs=3;  //x0, objective, *ff
				//A, P, C
				if(!csem_isnan(*matrixA)) ++num_objs;
				if(!csem_isnan(*P)) ++num_objs;
				if(!csem_isnan(*C)) ++num_objs;

				if(iagflg) {
						++num_objs;
						if(!iahflg && want_hessian) {//we need to compute Hessian if it is not provided.)
								fdhess(n, x0, fpls, (fcn_p)msem_fcn, state, a,  n,  &wrk[0], &wrk[n], ndigit, typsiz );
						for (i = 0; i < n; i++)
								for (j = 0; j < i; j++)
										a[i + j * n] = a[j + i * n];
						}
				} else if(want_hessian) {
						fdhess(n, x0, fpls, (fcn_p)msem_fcn, state, a,  n,  &wrk[0], &wrk[n], ndigit, typsiz );
						for (i = 0; i < n; i++)
								for (j = 0; j < i; j++)
										a[i + j * n] = a[j + i * n];
				}

				if(want_hessian) ++num_objs;

				PROTECT(value = allocVector(VECSXP, num_objs));
				PROTECT(names = allocVector(STRSXP, num_objs));
				k = 0;

				SET_STRING_ELT(names, k, mkChar("minimum"));
				SET_VECTOR_ELT(value, k, ScalarReal(fpls));
				k++;

				SET_STRING_ELT(names, k, mkChar("estimate"));
				SET_VECTOR_ELT(value, k, allocVector(REALSXP, n));
				for (i = 0; i < n; i++)
						REAL(VECTOR_ELT(value, k))[i] = x0[i];
				k++;

				if(iagflg) {
						SET_STRING_ELT(names, k, mkChar("gradient"));
						SET_VECTOR_ELT(value, k, allocVector(REALSXP, n));
						for (i = 0; i < n; i++)
								REAL(VECTOR_ELT(value, k))[i] = gpls[i];
						k++;
				}

				if(want_hessian){
						SET_STRING_ELT(names, k, mkChar("hessian"));
						SET_VECTOR_ELT(value, k, allocMatrix(REALSXP, n, n));
						for (i = 0; i < n * n; i++)
								REAL(VECTOR_ELT(value, k))[i] = a[i];
						k++;
				}

				/* A */
				if(!csem_isnan(*matrixA)) {
						SET_STRING_ELT(names, k, mkChar("A"));
						SET_VECTOR_ELT(value, k, allocMatrix(REALSXP, sizeAP, 1));
						for (i = 0; i < sizeAP; i++)
								REAL(VECTOR_ELT(value, k))[i] = matrixA[i];
						k++;
				}

				/* P */
				if(!csem_isnan(*P)) {
						SET_STRING_ELT(names, k, mkChar("P"));
						SET_VECTOR_ELT(value, k, allocMatrix(REALSXP, sizeAP, 1));
						for (i = 0; i < sizeAP; i++)
								REAL(VECTOR_ELT(value, k))[i] = P[i];
						k++;
				}

				/* C */
				if(!csem_isnan(*C)) { 
						SET_STRING_ELT(names, k, mkChar("C"));
						SET_VECTOR_ELT(value, k, allocMatrix(REALSXP, sizeC, 1));
						for (i = 0; i < sizeC; i++)
								REAL(VECTOR_ELT(value, k))[i] = C[i];
						k++;
				}

				SET_STRING_ELT(names, k, mkChar("f"));
				SET_VECTOR_ELT(value, k, allocMatrix(REALSXP, 1, state->model->G));
				for (i = 0; i < state->model->G; i++)
						REAL(VECTOR_ELT(value, k))[i] = ff[i];
				k++;

				setAttrib(value, R_NamesSymbol, names);
				UNPROTECT(3);
		}
		else {


				/*
				 *	 Dennis + Schnabel Minimizer
				 *
				 *	  SUBROUTINE OPTIF9(NR,N,X,FCN,D1FCN,D2FCN,TYPSIZ,FSCALE,
				 *	 +	   METHOD,IEXP,MSG,NDIGIT,ITNLIM,IAGFLG,IAHFLG,IPR,
				 *	 +	   DLT,GRADTL,STEPMX,STEPTOL,
				 *	 +	   XPLS,FPLS,GPLS,ITRMCD,A,WRK)
				 *
				 *
				 *	 Note: I have figured out what msg does.
				 *	 It is actually a sum of bit flags as follows
				 *	   1 = don't check/warn for 1-d problems
				 *	   2 = don't check analytic gradients
				 *	   4 = don't check analytic hessians
				 *	   8 = don't print start and end info
				 *	  16 = print at every iteration
				 *	 Using msg=9 is absolutely minimal
				 *	 I think we always check gradients and hessians
				 */


				optif9(n, n, x, (fcn_p) msem_fcn, (fcn_p) msem_Cd1fcn, (d2fcn_p) msem_Cd2fcn,
								state, typsiz, fscale, method, iexp, &msg, ndigit, itnlim,
								iagflg, iahflg, dlt, gradtl, stepmx, steptol, xpls, &fpls,
								gpls, &code, a, wrk, &itncnt);

				if(SEM_DEBUG) {
						if(state->have_hessian) 
								Rprintf("Hessian is provided.\n");
						else 
								Rprintf("Hessian is not provided.\n");

						Rprintf("The number of function evaluations: [%d]\n", state->n_eval);

				}

				if (msg < 0)
						opterror(msg);
				if (code != 0 && (omsg&8) == 0)
						optcode(code);

				int num_objs = 6;  //  ff 

				if (want_hessian) {
						++num_objs;
						fdhess(n, xpls, fpls, (fcn_p) msem_fcn, state, a, n, &wrk[0], &wrk[n],
										ndigit, typsiz);
						for (i = 0; i < n; i++)
								for (j = 0; j < i; j++)
										a[i + j * n] = a[j + i * n];
				}

				//A, P, C
				int sizeAP = state->sizeAP;
				int sizeC = state->sizeC;
				double *matrixA = (double *)R_alloc(sizeAP, sizeof(double)); //After the compuation,  C is n-by-n.
				double *P = (double *)R_alloc(sizeAP, sizeof(double)); //After the compuation,  C is n-by-n.
				double *C = (double *)R_alloc(sizeC, sizeof(double)); //After the compuation,  C is n-by-n.
				double *ff = (double *)R_alloc(state->model->G, sizeof(double));
				msem_returnAPCfcn(n, xpls, matrixA, P, C,ff,  state);
				if(!csem_isnan(*matrixA)) ++num_objs;
				if(!csem_isnan(*P)) ++num_objs;
				if(!csem_isnan(*C)) ++num_objs;

				PROTECT(value = allocVector(VECSXP, num_objs));
				PROTECT(names = allocVector(STRSXP, num_objs));
				k = 0;

				SET_STRING_ELT(names, k, mkChar("minimum"));
				SET_VECTOR_ELT(value, k, ScalarReal(fpls));
				k++;

				SET_STRING_ELT(names, k, mkChar("estimate"));
				SET_VECTOR_ELT(value, k, allocVector(REALSXP, n));
				for (i = 0; i < n; i++)
						REAL(VECTOR_ELT(value, k))[i] = xpls[i];
				k++;

				SET_STRING_ELT(names, k, mkChar("gradient"));
				SET_VECTOR_ELT(value, k, allocVector(REALSXP, n));
				for (i = 0; i < n; i++)
						REAL(VECTOR_ELT(value, k))[i] = gpls[i];
				k++;

				if (want_hessian) {
						SET_STRING_ELT(names, k, mkChar("hessian"));
						SET_VECTOR_ELT(value, k, allocMatrix(REALSXP, n, n));
						for (i = 0; i < n * n; i++)
								REAL(VECTOR_ELT(value, k))[i] = a[i];
						k++;
				}

				SET_STRING_ELT(names, k, mkChar("code"));
				SET_VECTOR_ELT(value, k, allocVector(INTSXP, 1));
				INTEGER(VECTOR_ELT(value, k))[0] = code;
				k++;

				SET_STRING_ELT(names, k, mkChar("iterations"));
				SET_VECTOR_ELT(value, k, allocVector(INTSXP, 1));
				INTEGER(VECTOR_ELT(value, k))[0] = itncnt;
				k++;

				/* A */
				if(!csem_isnan(*matrixA)) {
						SET_STRING_ELT(names, k, mkChar("A"));
						SET_VECTOR_ELT(value, k, allocMatrix(REALSXP, sizeAP, 1));
						for (i = 0; i < sizeAP; i++)
								REAL(VECTOR_ELT(value, k))[i] = matrixA[i];
						k++;
				}

				/* P */
				if(!csem_isnan(*P)) {
						SET_STRING_ELT(names, k, mkChar("P"));
						SET_VECTOR_ELT(value, k, allocMatrix(REALSXP, sizeAP, 1));
						for (i = 0; i < sizeAP; i++)
								REAL(VECTOR_ELT(value, k))[i] = P[i];
						k++;
				}

				/* C */
				if(!csem_isnan(*C)) {
						SET_STRING_ELT(names, k, mkChar("C"));
						SET_VECTOR_ELT(value, k, allocMatrix(REALSXP, sizeC, 1));
						for (i = 0; i < sizeC; i++)
								REAL(VECTOR_ELT(value, k))[i] = C[i];
						k++;
				}

				/* ff */
				SET_STRING_ELT(names, k, mkChar("f"));
				SET_VECTOR_ELT(value, k, allocMatrix(REALSXP, 1, state->model->G));
				for (i = 0; i < state->model->G; i++)
						REAL(VECTOR_ELT(value, k))[i] = ff[i];
				k++;

				setAttrib(value, R_NamesSymbol, names);
				UNPROTECT(3);
		}

    return value;
}
