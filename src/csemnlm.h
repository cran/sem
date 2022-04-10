/*
 * =====================================================================================
 *
 *       Filename:  csemnlm.h
 *
 *    Description:  csemnlm
 *
 *        Version:  1.0
 *        Created:  27/12/2011 04:31:12
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Zhenghua Nie (ZHN), zhenghua.nie@gmail.com
 *        Company:  McMaster University
 *
 *    Copyright (C) 2011 Zhenghua Nie. All Rights Reserved.
 *    This code is published under GNU GENERAL PUBLIC LICENSE.
 *
 *    This program is free software; you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation; either version 3 of the License,  or
 *    (at your option) any later version.
 *      
 *    This program is distributed WITHOUT ANY WARRANTY. See the
 *    GNU General Public License for more details.
 *           
 *    If you do not have a copy of the GNU General Public License,  
 *    write to the Free Software Foundation, Inc., 
 *    59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *              
 *
 *
 *        
 * =====================================================================================
 */

#ifndef __CSEMNLM_HPP__
#define __CSEMNLM_HPP__

//#ifndef USE_FC_LEN_T
//#define USE_FC_LEN_T
//#endif

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <R.h>
#include <Rdefines.h>
#include <R_ext/Utils.h>
#include <float.h>		/* for DBL_MAX */
#include <Rmath.h>
#include <R_ext/Applic.h>
#include <R_ext/RS.h>	       	/* for Memcpy */
#include <R_ext/BLAS.h>
#include <R_ext/Linpack.h>

#ifndef FCONE
//#define FCONE
#endif

#ifdef DEBUGSEM
#define SEM_DEBUG 1
#else
#define SEM_DEBUG 0
#endif 

// the structure for the solution of sem  (object returned by sem)
// currently,  we onlu use it for calculating logLik for missing data (solved by objectiveFIML)
typedef struct sem_solution_Info {
	SEXP data;
	SEXP pattern_number; 
	SEXP valid_data_patterns;
	SEXP tri;
	int posn_intercept;
	int t;
} sem_object;

/* 
    "data" = model$data, 
		"pattern.number" = model$model.number, 
		"valid.data.patterns"= model$valid.data.patterns, 
		"S" = model$S, 
		"logdetS" = model$logdetS, 
		"invS" = model$invS, 
		"N" = as.integer(model$N), 
		"m" = as.integer(model$m), 
		"n" = as.integer(model$n), 
		"t" = as.integer(model$t), 
		"fixed" = model$fixed, 
		"ram" = model$ram, 
		"sel.free" = model$sel.free, 
		"arrows.1" = model$arrows.1, 
		"arrows.1.free" = model$arrows.1.free, 
		"one.head" = model$one.head, 
		"arrows.2t" = model$arrows.2t, 
		"arrows.2" = model$arrows.2, 
		"arrows.2.free" = model$arrows.2.free, 
		"unique.free.1" = model$unique.free.1, 
		"unique.free.2" = model$unique.free.2, 
		"J" = model$J, 
		"correct" = model$correct, 
		"param.names" = model$param.names, 
		"var.names" = model$var.names, 
		"one.free" = model$one.free, 
		"two.free" = model$two.free, 
		"raw" = model$raw)
 * */
typedef struct model_Info {
	SEXP data;
	SEXP pattern_number; 
	SEXP valid_data_patterns;
  SEXP S; //n-by-n
	double logdetS;
	SEXP invS; //n-by-n
	int N;
	int m;
	int n;
	int t;
	SEXP fixed;  //vector, t+m-n
	SEXP ram;
	SEXP sel_free;
	SEXP arrows_1;
	SEXP arrows_1_free;
	SEXP one_head;
	SEXP arrows_2t;
	SEXP arrows_2;
	SEXP arrows_2_free;
	SEXP unique_free_1;
	SEXP unique_free_2;
	SEXP J;   //n-by-m
	SEXP correct; ///m-by-m
	SEXP param_names;
	SEXP var_names;
	SEXP one_free;
	SEXP two_free;
	int raw;
	int *arrows_1_seq;
	int *arrows_2_seq;
	sem_object *semObject;
} model_info;

typedef struct msem_model_Info {
	int G;  //number of groups
	SEXP data;
	SEXP pattern_number; 
	SEXP valid_data_patterns;
  SEXP S; //n-by-n
	SEXP logdetS;
	SEXP invS; //n-by-n
	SEXP N;
	SEXP m;
	SEXP n;
	int t;
	SEXP fixed;  //vector, t+m-n
	SEXP ram;
	SEXP sel_free;
	SEXP arrows_1;
	SEXP arrows_1_free;
	SEXP one_head;
	SEXP arrows_2t;
	SEXP arrows_2;
	SEXP arrows_2_free;
	SEXP unique_free_1;
	SEXP unique_free_2;
	SEXP J;   //n-by-m
	SEXP correct; ///m-by-m
	SEXP param_names;
	SEXP var_names;
	SEXP one_free;
	SEXP two_free;
	int raw;
	SEXP arrows_1_seq;
	SEXP arrows_2_seq;
	model_info *gmodel;  //a pointer for each group's model. 
} msem_model_info;

//this will define the protocol of our objective function
typedef void (*myfcn_p)(int, const double *,  double *, double *, double *,  double *, double *, double *, void *);
//typedef void (*myfcn_p)(int n,  double *x,  double *f, double *g, double *h,  double *A,  double *P, double *C, void *state);

typedef void (*msem_fcn_p)(int, const double *,  double *, double *, double *,  double *, double *, double *, double *, void *);
//typedef void (*msem_fcn_p)(int n,  double *x,  double *f, double *g, double *h,  double *A,  double *P, double *C, double *ff, void *state);

#define FT_SIZE 3		/* default size of table to store computed
				   function values */

typedef struct {
  double fval;
  double *x;
  double *grad;
  double *hess;
	double *C;
	double *A;
	double *P;
} ftable;

typedef struct {
  double fval;
  double *x;
  double *grad;
  double *hess;
	double *C;
	double *A;
	double *P;
	double *ff;
} msem_ftable;

typedef struct {
  int  n_eval;	      /* the number of evaluations of the objective function. */
  myfcn_p *myobjfun; 
  int have_gradient;
  int have_hessian;
/*  int n;	      -* length of the parameter (x) vector */
  int FT_size;	      /* size of table to store computed
			 function values */
  int FT_last;	      /* Newest entry in the table */
  ftable *Ftable;
	model_info *model;   /*including sem_object*/
} function_info;

typedef struct {
  int  n_eval;	      /* the number of evaluations of the objective function. */
  msem_fcn_p *myobjfun; 
  int have_gradient;
  int have_hessian;
/*  int n;	      -* length of the parameter (x) vector */
  int FT_size;	      /* size of table to store computed
			 function values */
  int FT_last;	      /* Newest entry in the table */
  msem_ftable *Ftable;
	msem_model_info *model;
	int sizeAP;
	int sizeC;
} msem_function_info;

#ifdef __cplusplus
 extern "C" {
#endif


void fdhess(int n, double *x, double fval, fcn_p fun, void *state,
	    double *h, int nfd, double *step, double *f,
	    int ndigit, double *typx);
void
optif9(int nr, int n, double *x, fcn_p fcn, fcn_p d1fcn, d2fcn_p d2fcn,
       void *state, double *typsiz, double fscale, int method,
       int iexp, int *msg, int ndigit, int itnlim, int iagflg, int iahflg,
       double dlt, double gradtl, double stepmx, double steptl,
       double *xpls, double *fpls, double *gpls, int *itrmcd, double *a,
       double *wrk, int *itncnt);
SEXP csemnlm(double *x0, int n, int iagflg,  int iahflg, int want_hessian, 
				double *typsize, double fscale, int msg, int ndigit, double gradtl, 
				double stepmx, double steptol, int itnlim, void *model, myfcn_p myobjfun,
				int optimize);

SEXP cmsemnlm(double *x0, int n, int iagflg,  int iahflg, int want_hessian, 
				double *typsize, double fscale, int msg, int ndigit, double gradtl, 
				double stepmx, double steptol, int itnlim, msem_model_info *model, msem_fcn_p myobjfun,
				int optimize);

 #ifdef __cplusplus
 }
 #endif

#endif
