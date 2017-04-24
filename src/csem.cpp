/*
 * =====================================================================================
 *
 *       Filename:  csem.cpp
 *
 *    Description:  csem
 *
 *        Version:  1.0
 *        Created:  27/12/2011 00:28:29
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


#include "csem.h"
#include "utils.h"

using namespace std;

// Environment for evaluation and object hunting
static SEXP theenv;
static SEXP thefun;   // call R functions. Currently,  we call "print" in R.

static double csem_NaN;
static const double log2Pi = 1.83787706640934533908193770912476;   //log(2*pi)

void printSEXP(SEXP sexp, const string msg);

void semprintRealVector(const double *x,  int n,  int index)
{
		SEXP rargs,Rcall,result;

		// Allocate memory for a vector of reals.
		// This vector will contain the elements of x,
		// x is the argument to the R function R_eval_f
		PROTECT(rargs = allocVector(REALSXP,n));
		for (int i=0;i<n;i++) {
				REAL(rargs)[i] = x[i];
		}

		// evaluate R function R_eval_f with the control x as an argument
		PROTECT(Rcall = lang2(thefun,rargs));
		PROTECT(result = eval(Rcall,theenv));

		UNPROTECT(3);
		return;

}
//
// Extracts element with name 'str' from R object 'list'
// and returns that element.
//
SEXP getListElement(SEXP list,   int ind)
{
		SEXP elmt = R_NilValue;

		if(ind >= 0 && ind < length(list))
		{
				elmt = VECTOR_ELT(list, ind);
		}
		else
				error(("The index is not in the range of the list."));

		return elmt;
}

SEXP getListElement(SEXP list, std::string str)
{
		SEXP elmt = R_NilValue, names = getAttrib(list, R_NamesSymbol);
		int i;
		for (i = 0; i < length(list); i++)
				if(str.compare(CHAR(STRING_ELT(names, i))) == 0) {
						elmt = VECTOR_ELT(list, i);
						break;
				}
		return elmt;
}

double  getVectorElement(SEXP vect,  int ind )
{
		double elmt = csem_NaN;
		if(ind >= 0 && ind < length(vect))
				elmt = REAL(AS_NUMERIC(vect))[ind];
		else
				error(("The index is not in the range of the vector."));

		return elmt;
}
// if ind==-1,  we will search the names,  and then return the object.
double  getVectorElement(SEXP vect,  std::string str )
{
		SEXP names = getAttrib(vect, R_NamesSymbol);
		double elmt = csem_NaN;
		int i;
		for (i = 0; i < length(vect); i++)
				if(str.compare(CHAR(STRING_ELT(names, i))) == 0) {
						elmt = REAL(AS_NUMERIC(vect))[i];
						break;
				}
		return elmt;
}

SEXP generateMatrix(double *A, int nrow,  int ncol)
{
		SEXP elmt;
		elmt = allocMatrix(REALSXP, nrow, ncol);
		for(int i=0; i < nrow*ncol; ++i) 
				REAL(elmt)[i] = A[i];
		return(elmt);
}
/*
 * We commented this function. This function is for debugging.
SEXP showArgs1(SEXP largs)
{
		int i, nargs = LENGTH(largs);
		Rcomplex cpl;
		SEXP el, names;
		PROTECT(names = getAttrib(largs, R_NamesSymbol));

		const char *name;

		for(i = 0; i < nargs; i++) {
				el = VECTOR_ELT(largs, i);
				name = isNull(names) ? "" : CHAR(STRING_ELT(names, i));
				switch(TYPEOF(el)) {
						case REALSXP:
								Rprintf("[%d] '%s' %f\n", i+1, name, REAL(el)[0]);
								break;
						case LGLSXP:
						case INTSXP:
								Rprintf("[%d] '%s' %d\n", i+1, name, INTEGER(el)[0]);
								break;
						case CPLXSXP:
								cpl = COMPLEX(el)[0];
								Rprintf("[%d] '%s' %f + %fi\n", i+1, name, cpl.r, cpl.i);
								break;
						case STRSXP:
								Rprintf("[%d] '%s' %s\n", i+1, name,
												CHAR(STRING_ELT(el, 0)));
								break;
						default:
								Rprintf("[%d] '%s' R type\n", i+1, name);
				}
		}
		UNPROTECT(1);

		return(R_NilValue);
}
*/
//
void setApplicationOptions(int &hessian, double &fscale, double &steptol, double &stepmax, int &iterlim, int &ndigit, 
				int &print_level, int &check_analyticals,  double &gradtol, SEXP opts ) {

		const int options_integer_num = 5;
		const string option_names_integer[options_integer_num]={
				"hessian", 
				"iterlim", 
				"ndigit", 
				"print.level", 
				"check.analyticals" 
		};
		int option_integer[options_integer_num];
		// Default value
		option_integer[0] = 0;  //hessian
		option_integer[1] = 100;  //iterlim
		option_integer[2] = 12; //ndigit
		option_integer[3] = 0;  //print_level
		option_integer[4] = 1;  //check.analytics

		const int options_numeric_num = 4;
		const string option_names_numeric[options_numeric_num]={
				"fscale", 
				"steptol", 
				"stepmax", 
				"gradtol"  
		};
		double option_double[options_numeric_num];
		// Default Values
		option_double[0] = 1.0; //fscale
		option_double[1] = 1.0e-6;  //steptol
		option_double[2] = stepmax;   //stepmax,  if not given,  max(1000*sqrt(sum(x0/typsiz)^2), 1000)
		option_double[3] = 1.0e-6; //gradtol

		// extract the sub-lists with options of the different types into separate lists
		SEXP opts_integer = getListElement(opts, "integer");
		SEXP opts_numeric = getListElement(opts, "numeric");
		SEXP opts_string = getListElement(opts, "string");

		// loop over the integer options and set them
		SEXP opts_integer_names;
		PROTECT(opts_integer_names = getAttrib(opts_integer, R_NamesSymbol));
		for (int list_cnt=0;list_cnt<length( opts_integer );list_cnt++) {

				SEXP opt_value;
				PROTECT(opt_value = AS_INTEGER(VECTOR_ELT(opts_integer, list_cnt)));

				string names = CHAR(STRING_ELT(opts_integer_names, list_cnt));
				for (int i=0;i<options_integer_num;++i) {
						if(names == option_names_integer[i]) option_integer[i] = INTEGER(opt_value)[0];
				}
				UNPROTECT(1);	
		}
		UNPROTECT(1);

		// loop over the numeric options and set them
		SEXP opts_numeric_names;
		PROTECT(opts_numeric_names = getAttrib(opts_numeric, R_NamesSymbol));
		for (int list_cnt=0;list_cnt<length( opts_numeric );list_cnt++) {

				SEXP opt_value;
				PROTECT(opt_value = VECTOR_ELT(opts_numeric, list_cnt));

				string names = CHAR(STRING_ELT(opts_numeric_names, list_cnt));
				for(int i=0; i< options_numeric_num; ++i) {
						if(names == option_names_numeric[i]) option_double[i] = REAL(opt_value)[0];
				}
				UNPROTECT(1);	
		}
		UNPROTECT(1);

		// loop over the string options and set them
		// Currently,  we don't have string options, so we commented the following two lines.
		// SEXP opts_string_names;
		// opts_string_names = getAttrib(opts_string, R_NamesSymbol);
		for (int list_cnt=0;list_cnt<length( opts_string );list_cnt++) {

				// opt_value will contain the first (should be the only one) element of the list
				SEXP opt_value;
				PROTECT(opt_value = STRING_ELT(VECTOR_ELT(opts_string, list_cnt),0));
				//we don't have string options.

				UNPROTECT(1);	
		}

		hessian = option_integer[0];
		iterlim = option_integer[1];
		ndigit = option_integer[2];
		print_level = option_integer[3];
		check_analyticals = option_integer[4];

		fscale = option_double[0];
		steptol = option_double[1];
		stepmax = option_double[2];
		gradtol = option_double[3];

		return;

}
//print an R object
void printSEXP(SEXP sexp, const string msg)
{
		Rprintf("%s\n", msg.c_str());
		PrintValue(sexp);

		return;
}
//print a mtrix (double)
void printMatrix(const double *x, int row, int col, const string msg, int column_wise) 
{
		Rprintf("Matrix: %s [%d*%d]\t", msg.c_str(), row, col);
		if(column_wise) {
				Rprintf("(column-wise)\n");
				for(int i=0;i<row;++i){
						for(int j=0;j<col; ++j)
								Rprintf("%f ", x[i+j*row]);
						Rprintf("\n");
				}
		}
		else {
				Rprintf("(row-wise)\n");
				for(int i=0;i<row;++i){
						for(int j=0;j<col; ++j)
								Rprintf("%f ", x[i*col+j]);
						Rprintf("\n");
				}
		}
}
//
//print a mtrix (double)
void printMatrix(const int *x, int row, int col, const string msg, int column_wise) 
{
		Rprintf("Matrix: %s [%d*%d]\t", msg.c_str(), row, col);
		if(column_wise) {
				Rprintf("(column-wise)\n");
				for(int i=0;i<row;++i){
						for(int j=0;j<col; ++j)
								Rprintf("%d ", x[i+j*row]);
						Rprintf("\n");
				}
		}
		else {
				Rprintf("(row-wise)\n");
				for(int i=0;i<row;++i){
						for(int j=0;j<col; ++j)
								Rprintf("%d ", x[i*col+j]);
						Rprintf("\n");
				}
		}
}

// A = A^T
static void MatrixTranspose(double *A, const int &row,  const int &col)
{
		double *copy_A = new double[row*col];
		Memcpy(copy_A, A, row*col);

		for(int i =0; i < row*col; ++i)
		{
				int r = i / col;
				int c = i % col;

				A[c * row + r] = copy_A[i];
		}

		delete[] copy_A;
		return;
}

static void MatrixTranspose(int *A, const int &row,  const int &col)
{
		int *copy_A = new int[row*col];
		Memcpy(copy_A, A, row*col);

		for(int i =0; i < row*col; ++i)
		{
				int r = i / col;
				int c = i % col;

				A[c * row + r] = copy_A[i];
		}

		delete[] copy_A;
		return;
}

extern "C" {
		void test_objective(int n, const double x[],  double *f,  double *g, double *h,double *A,  double *P,  double *C, function_info *state)
		{
				int ind;
				*A = csem_NaN; 
				*P = csem_NaN;
				*C = csem_NaN;

				*f = 0.0;
				for(ind=0;ind<n;++ind) *f += (x[ind]-ind-1)*(x[ind]-ind-1); 

				if(state->have_gradient) {
						for(ind = 0; ind < n; ++ind){
								g[ind] = 2.0*(x[ind]-ind-1);
						}
						//The interested thing is that we need more function evaluations when Hessian is provided (n_eval:33).
						//If not provided,  n_eval is 18. IS there anything wrong?
						if(state->have_hessian){
								for(ind = 0; ind < n; ++ind){
										h[ind*n+ind] = 2.0;
								}
						}
				}

				return;
		}

		void msem_test_objective(int n, const double x[],  double *f,  double *g, double *h,double *A,  double *P,  double *C, double *ff, msem_function_info *state)
		{
				int ind;
				*A = csem_NaN; 
				*P = csem_NaN;
				*C = csem_NaN;
				*ff = csem_NaN;

				*f = 0.0;
				for(ind=0;ind<n;++ind) *f += (x[ind]-ind-1)*(x[ind]-ind-1); 

				if(state->have_gradient) {
						for(ind = 0; ind < n; ++ind){
								g[ind] = 2.0*(x[ind]-ind-1);
						}
						//The interested thing is that we need more function evaluations when Hessian is provided (n_eval:33).
						//If not provided,  n_eval is 18. IS there anything wrong?
						if(state->have_hessian){
								for(ind = 0; ind < n; ++ind){
										h[ind*n+ind] = 2.0;
								}
						}
				}

				return;
		}


		static void MatrixInverse(double* A,  const int n)
		{
				int *IPIV = new int[n+1];
				int LWORK = n*n;
				double *WORK = new double[LWORK];
				int INFO;

				F77_CALL(dgetrf)(&n, &n, A, &n, IPIV, &INFO);
				if(INFO != 0) {
						error(( "The matrix is non-invertable."));
				}
				F77_CALL(dgetri)(&n, A, &n, IPIV, WORK, &LWORK, &INFO);

				delete[] IPIV;
				delete[] WORK;
				return;
		}

		static double MatrixDeterminant(const double *A,const int nrow,const int ncol)
		{
				double det;

				double *tA;  //temp matrix	
				int *IPIV;
				int INFO;

				if(nrow != ncol) {
						error(("We cannot comptue the determinant of a non-square matrix.\n"));
				}

				IPIV = new int[nrow+1];
				tA = new double[nrow*nrow];
				Memcpy(tA, A, nrow*ncol);


				F77_CALL(dgetrf)(&nrow, &ncol, tA, &ncol, IPIV, &INFO);

				if(INFO != 0) {
						error(("Nonsingular matrix." ));
				}

				det = 1;
				for(int i = 0; i < nrow; i++) {
						det *= tA[i*nrow+ i];  //diagonal
						if(IPIV[i] != (i+1)) 
								det *= -1;
				}

				delete[] tA;
				delete[] IPIV;

				return(det);
		}

		static double MatrixTrace(const double *A, int nrow, int ncol)
		{
				double trace;

				if(nrow != ncol) {
						error(("We fail to comptue a trace of a non-square matrix.\n"));
				}

				trace = 0.0;
				for(int j = 0; j < nrow; j++) {
						trace += A[j*nrow+j];
				}

				return(trace);
		}

		//C=A.B ,  no transpose.
		static void MatrixMulti(const double *A, int rowA, int colA, const double *B,  int rowB,  int colB,  double *C)
		{
				if(SEM_DEBUG) Rprintf("A[%4d-by-%4d]*B[%4d-by-%4d]=C[%4d-by-%4d]\n", rowA, colA, rowB, colB, rowA, colB);

				if(colA != rowB) error(("The matrices are not conformable."));
				///			  C := alpha*op( A )*op( B ) + beta*C, 
				// R_ext/Blas.h
				// F77_NAME(dgemm)(const char *transa,  const char *transb,  const int *m, 
				//     const int *n,  const int *k,  const double *alpha, 
				//     const double *a,  const int *lda, 
				//     const double *b,  const int *ldb, 
				//     const double *beta,  double *c,  const int *ldc);
				memset(C, 0, rowA*colB*sizeof(double));
				char Tran = 'n';
				const double alpha = 1.0;
				const double beta = 0.0;
				F77_CALL(dgemm)(&Tran, &Tran, &rowA, &colB, &colA, &alpha, A, &rowA, B, &colA, &beta, C, &rowA);

				return;
		}

		//C=A.Transpose(B).
		static void MatrixMultiTransB(const double *A, int rowA, int colA, const double *B,  int rowB,  int colB,  double *C)
		{
				if(SEM_DEBUG) Rprintf("A[%4d-by-%4d]*B'[%4d-by-%4d]=C[%4d-by-%4d]\n", rowA, colA, rowB, colB, rowA, rowB);

				if(colA != colB) error(("The matrices are not conformable."));
				///			  C := alpha*op( A )*op( B ) + beta*C, 
				// R_ext/Blas.h
				// F77_NAME(dgemm)(const char *transa,  const char *transb,  const int *m, 
				//     const int *n,  const int *k,  const double *alpha, 
				//     const double *a,  const int *lda, 
				//     const double *b,  const int *ldb, 
				//     const double *beta,  double *c,  const int *ldc);
				memset(C, 0, rowA*rowB*sizeof(double));
				char ATran = 'n';
				char BTran = 't';
				const double alpha = 1.0;
				const double beta = 0.0;
				F77_CALL(dgemm)(&ATran, &BTran, &rowA, &rowB, &colA, &alpha, A, &rowA, B, &rowB, &beta, C, &rowA);

				return;
		}

		//C=Transpose(A).Transpose(B) .
		static void MatrixMultiTransAB(const double *A, int rowA, int colA, const double *B,  int rowB,  int colB,  double *C)
		{
				if(SEM_DEBUG) Rprintf("A'[%4d-by-%4d]*B'[%4d-by-%4d]=C[%4d-by-%4d]\n", rowA, colA, rowB, colB, colA, rowB);

				if(colA != colB) error(("The matrices are not conformable."));
				///			  C := alpha*op( A )*op( B ) + beta*C, 
				// R_ext/Blas.h
				// F77_NAME(dgemm)(const char *transa,  const char *transb,  const int *m, 
				//     const int *n,  const int *k,  const double *alpha, 
				//     const double *a,  const int *lda, 
				//     const double *b,  const int *ldb, 
				//     const double *beta,  double *c,  const int *ldc);
				memset(C, 0, colA*rowB*sizeof(double));
				char ATran = 't';
				char BTran = 't';
				const double alpha = 1.0;
				const double beta = 0.0;
				F77_CALL(dgemm)(&ATran, &BTran, &colA, &rowB, &rowA, &alpha, A, &rowA, B, &rowB, &beta, C, &rowA);

				return;
		}

		//The fowllowing function may be rewritten as similar as ddot.f
		static void sempdot(const int *n, const double *x, const int *incx, const double *y, const int *incy, double *z)
		{
				for(int i=0; i < *n; ++i) *z++ = *x++ * *y++;
				return;
		}


//column_wise
static double *Kronecker(const double *A,  const int &rowA, const int &colA, const double *B, const int &rowB, const int &colB)
{
		double *kron = new double[rowA*colA*rowB*colB];
		int iA, jA, iB, jB;
		int rowC = rowA*rowB;

		for(iA = 0; iA < colA; ++iA)
		{
				for(jA = 0; jA < rowA; ++jA)
				{
						double a = *A++;   // the matrix is column-wise.
						const double *b = &B[0];
						int M = iA * colB;
						int N = jA * rowB;
						for(iB = 0; iB < colB; ++iB)
						{
								for(jB = 0; jB < rowB; ++jB)
								{
										kron[(M+iB)*rowC+N+jB] = a * *b++;   // the matrix is column-wise
								}
						}
				}
		}

		return kron;
}

		/* 
#this function will transfer into C code for generating matrix A and P with new parameters.
# the index of par should be the same as start,  so we set par's names as start's names.
generate.AP <- function(par, start=start, model=model){
names(par) <- names(start)
A <- P <- matrix(0, model$m, model$m)
val <- ifelse(model$fixed, model$ram[, 5], par[model$sel.free])
A[model$arrows.1]<-val[model$one.head]
P[model$arrows.2t] <- P[model$arrows.2] <- val[!model$one.head]
A <- diag(model$m)-A  #please note here.We don't need A,  we always want I-A. 
AP <- list(A=A, P=P)
AP
}
		 * */
static void generate_AP(int n,  const double x[], double *A, double *P,  double *ImA, model_info *model)
{
		int n_val = length(model->fixed);
		int *fixed = new int[n_val];
		int *sel_free = new int[length(model->sel_free)];
		double *ram5 = new double[nrows(model->ram)];
		double *val = new double[n_val];

		Memcpy(fixed, INTEGER(AS_INTEGER(model->fixed)), n_val);
		Memcpy(sel_free, INTEGER(AS_INTEGER(model->sel_free)), length(model->sel_free));
		Memcpy(ram5, REAL(AS_NUMERIC(model->ram))+4*nrows(model->ram), nrows(model->ram));

		for(int i=0;i<length(model->fixed); ++i){
				val[i] = (fixed[i]==1) ? ram5[i] : x[sel_free[i]-1]; //fortran to C
		}

		int *one_head = new int[length(model->one_head)];
		double *val_one_head = new double[n_val];
		double *val_two_head = new double[n_val];

		Memcpy(one_head, INTEGER(AS_INTEGER(model->one_head)), length(model->one_head));

		int ind_one = 0;
		int ind_two = 0;
		for(int i=0; i< n_val; ++i){
				if(one_head[i]==1) {
						val_one_head[ind_one] = val[i];
						++ind_one;
				}
				else {
						val_two_head[ind_two] = val[i];
						++ind_two;
				}
		}

		int m = model->m;
		memset(A, 0, m*m*sizeof(double));
		memset(P, 0, m*m*sizeof(double));
		memset(ImA, 0, m*m*sizeof(double));

		//A and diag(m)-A
		int nA = length(model->arrows_1)/2;
		int nP = length(model->arrows_2)/2;
		int *tA = new int[max(nA*2, nP*2)];

		Memcpy(tA, INTEGER(AS_INTEGER(model->arrows_1)), nA*2);

		for(int i=0;i<nA;++i){
				int ir = tA[i]-1;  //fortran to C
				int il = tA[nA+i]-1;
				A[ir+il*m] = val_one_head[i];  //column-wise.
				ImA[ir+il*m] = -val_one_head[i];  //column-wise.
		}
		for(int i=0;i<m;++i) ImA[i+i*m] = 1.0 + ImA[i+i*m];

		//P
		Memcpy(tA, INTEGER(AS_INTEGER(model->arrows_2)), nP*2);
		int *tAt = new int[nP*2];
		Memcpy(tAt, INTEGER(AS_INTEGER(model->arrows_2t)), nP*2);
		for(int i=0;i<nP;++i){
				int ir = tA[i]-1;  //fortran to C
				int il = tA[nP+i]-1;
				P[ir+il*m] = val_two_head[i];  //column-wise.
				ir = tAt[i]-1;  //fortran to C
				il = tAt[nP+i]-1;
				P[ir+il*m] = val_two_head[i];  //column-wise.
		}

		delete[] fixed;
		delete[] sel_free;
		delete[] ram5;
		delete[] val;
		delete[] one_head;
		delete[] val_one_head;
		delete[] val_two_head;
		delete[] tA;
		delete[] tAt;

		return;
}


// Please note that R uses the column-wise to store matrices.
// this is for valid.data.pattern[i, ]
// input: valid.data.pattern
static int *SubMatrixRow(SEXP A,  const int &row, const int &col, const int &ithrow)
{
		int *transA=new int[row*col];  // for row-wise
		int *rowA = new int[col];  // the ith row of the matrix A
		Memcpy(transA, INTEGER(AS_INTEGER(A)), row*col);

		MatrixTranspose(transA, col, row);   //column-wise

		Memcpy(rowA, &transA[ithrow*col], col);

		delete[] transA;

		return(rowA);

}
// Please note that R uses the column-wise to store matrices.
// this is for objectiveFIML
// x <- data[pattern.number==i, sel, drop=FALSE]
// C[sel, sel]
static double *SubMatrix(const double *A, const int *rA,  const int *cA, const int &row, const int &col,  int &row_subA, int &col_subA)
{
		double *tmpA = new double[row * col];

		row_subA = 0;
		col_subA = 0;

		if(SEM_DEBUG)
				printMatrix(A, row, col, "data/C", 1);

		for(int i = 0; i < col; ++i)
		{
				if(cA[i] == 1) // keep these columnss
				{
						Memcpy(&tmpA[col_subA * row],  &A[i*row], row);
						++col_subA;
				}
		}
		if(SEM_DEBUG)
				printMatrix(tmpA, row, col_subA, "select columns", 1);

		MatrixTranspose(tmpA, col_subA, row);

		if(SEM_DEBUG)
				printMatrix(tmpA,  row, col_subA,  "select columns", 1);

		double *subA = new double[col_subA * row];
		for(int i = 0; i< row; ++i)
		{
				if(rA[i] == 1) // keep these rows
				{
						Memcpy(&subA[col_subA*row_subA], &tmpA[i*col_subA], col_subA);
						++row_subA;
				}
		}
		if(SEM_DEBUG)
				printMatrix(subA, row_subA, col_subA, "select rows", 0 );

		MatrixTranspose(subA, row_subA, col_subA);

		delete[] tmpA;

		return(subA);
}

// commutation matrix
static double *CommutationMatrix(const int &m, const int &n)
{
		int p = m*n;
		double *C = new double[p*p];
		memset(C, 0, p*p*sizeof(double));  //zero matrix.

		int r, c, i, j;
		r = 0;
		for (i=0; i < m; ++i)
		{
				c = i;
				for (j = 0; j < n; ++j)
				{
						C[r*p+c] = 1.0;
						++r;
						c += m;
				}
		}

		return(C);
}

// vec
static double *vec(const double *A, const int &row, const int &col) 
{// the matrix is in the C order,  if the matrix is in the Fortran/R order,  we don't need this.
		double *vecA = new double[row*col];

		Memcpy(vecA, A, row*col);

		MatrixTranspose(vecA, row, col);

		return(vecA);
}

// extend matrix
static double *ExtendMatrix(const double *selA, const int &row, const int &col, const int *sel, const int &nsel)
{// selA is in the Fortran/R order (the column wise), 
		double *A = new double[nsel*nsel];
		double *B = new double[row*nsel];

		memset(A, 0, nsel*nsel*sizeof(double));
		memset(B, 0, nsel*row*sizeof(double));

		int iselA, i;
		iselA = 0;

		for(i = 0; i < nsel; ++i)
		{
				if(sel[i])
				{
						Memcpy(&B[i*row], &selA[iselA*row], row);
						++ iselA;
				}
		} //nsel*row
		MatrixTranspose(B, nsel, row);

		iselA = 0; 
		for(i = 0; i < nsel; ++i)
		{
				if(sel[i])
				{
						Memcpy(&A[i*nsel], &B[iselA*nsel], nsel);
						++ iselA;
				}
		}

		delete[] B;

		return(A);
}

// this function will compute the Maximum Likelihood and gradients.
void objectiveML(int n, const double x[], double *f, double *g, double *h,  double *A,  double *P, double *C, function_info *state)
{
		R_CheckUserInterrupt();

		model_info *model = state->model;

		int m = model->m;
		int modeln = model->n;
		int maxmn = (m > modeln ? m : modeln);


		double *ImA = new double[m*m];

		generate_AP(n, x, A, P, ImA, model);
		if(SEM_DEBUG) {
				printMatrix(A, m, m, "Matrix A", 1);
				printMatrix(ImA, m, m, "Matrix (I-A)", 1);
				printMatrix(P, m, m, "Matrix P", 1);
		}


		double *invA,  *invC;
		double *C0;  //for matrix multiplication.


		invC = new double[maxmn*maxmn]; //After the compuation,  C is n-by-n.
		C0 = new double[maxmn*maxmn];  //After the compuation,  C is n-by-n.
		memset(C, 0, maxmn*maxmn*sizeof(double));
		memset(C0, 0, maxmn*maxmn*sizeof(double));
		memset(invC, 0, maxmn*maxmn*sizeof(double));

		///Please be careful that R uses column-wise to store matrix.
		invA = new double[m*m];


		Memcpy(invA, ImA, m*m);
		MatrixInverse(invA, m);
		if(SEM_DEBUG) printMatrix(invA, m, m, "invA", 1);

		if(SEM_DEBUG) {
				MatrixMulti(ImA, m, m, invA, m, m, C0);
				printMatrix(C0, m, m, "A %*% invA", 1);
		}

		MatrixMulti(REAL(model->J), nrows(model->J), ncols(model->J), invA, m, m, C0);
		if(SEM_DEBUG) printMatrix(C0, modeln, m, "J %*% I.Ainv", 1);

		MatrixMulti(C0, modeln, m, P, m, m, C);
		if(SEM_DEBUG) printMatrix(C, modeln, m, "J %*% I.Ainv %*% P", 1);

		MatrixMultiTransB(C, modeln, m, invA, m, m, C0);
		if(SEM_DEBUG) printMatrix(C0, modeln, m, "J %*% I.Ainv %*% P %*% t(I.Ainv)", 1);

		MatrixMultiTransB(C0, modeln, m, REAL(model->J),  modeln, m,  C);
		if(SEM_DEBUG) printMatrix(C, modeln, modeln, "J %*% I.Ainv %*% P %*% t(I, Ainv) %*% t(J)", 1);

		Memcpy(invC, C, modeln*modeln);
		MatrixInverse(invC, modeln);
		if(SEM_DEBUG) {
				printMatrix(invC, modeln, modeln, "Cinv", 1);
		}

		MatrixMulti(REAL(model->S), modeln, modeln, invC, modeln, modeln, C0);
		if(SEM_DEBUG) printMatrix(C0, modeln, modeln, "S %*% Cinv", 1);
		*f = MatrixTrace(C0, modeln, modeln)+log(MatrixDeterminant(C, modeln, modeln))-modeln-model->logdetS;

		//now we start to calculate the gradient.

		if(state->have_gradient)
		{
				double *grad_P, *grad_A;

				grad_P = new double[maxmn*maxmn];  //After the compuation,  grad_P is m-by-m.
				grad_A = new double[maxmn*maxmn];  //After the compuation,  grad_A is m-by-m.
				memset(grad_P, 0, maxmn*maxmn*sizeof(double));
				memset(grad_A, 0, maxmn*maxmn*sizeof(double));

				MatrixMultiTransAB(invA, m, m, REAL(model->J), modeln, m, C0);
				if(SEM_DEBUG) printMatrix(C0, m, modeln, "t(I.Ainv) %*% t(J)", 1);

				MatrixMulti(C0, m, modeln, invC, modeln, modeln, grad_P);
				if(SEM_DEBUG) printMatrix(grad_P, m, modeln, "t(I.Ainv) %*% t(J) %*% Cinv", 1);

				Memcpy(grad_A, C, modeln*modeln); //y, we use daxpy,  y=ax+y
				Memcpy(C0, REAL(model->S), modeln*modeln); //y, we use daxpy,  y=ax+y
				//F77_NAME(daxpy)(const int *n,  const double *alpha, 
				//								const double *dx,  const int *incx, 
				//								double *dy,  const int *incy);
				double alpha = -1.0;
				int incx = 1;
				int mm = modeln*modeln;

				F77_CALL(daxpy)(&mm,&alpha, C0,&incx, grad_A, &incx); //grad_A = -C0 + grad_A ==>  grad_A = -S + C
				if(SEM_DEBUG) printMatrix(grad_A, modeln, modeln, "(C-S)", 1);

				MatrixMulti(grad_P, m, modeln, grad_A, modeln, modeln, C0);
				if(SEM_DEBUG) printMatrix(C0, m, modeln, "t(I.Ainv) %*% t(J) %*% Cinv %*% (C-S)", 1);

				MatrixMulti(C0, m, modeln, invC, modeln, modeln, grad_P);
				if(SEM_DEBUG) printMatrix(grad_P, m, modeln, "t(I.Ainv) %*% t(J) %*% Cinv %*% (C-S) %*% Cinv", 1);

				MatrixMulti(grad_P, m, modeln, REAL(model->J), modeln, m, C0);
				if(SEM_DEBUG) printMatrix(C0, m, m, "t(I.Ainv) %*% t(J) %*% Cinv %*% (C-S) %*% Cinv %*% J", 1);

				MatrixMulti(C0, m, m, invA, m, m, grad_A);
				if(SEM_DEBUG) printMatrix(grad_A, m, m, "t(I.Ainv) %*% t(J) %*% Cinv %*% (C-S) %*% Cinv %*% J %*% I.Ainv ", 1);

				mm = m*m;
				sempdot(&mm, REAL(model->correct), &incx, grad_A, &incx, grad_P);
				if(SEM_DEBUG) printMatrix(grad_P, m, m, "correct * t(I.Ainv) %*% t(J) %*% Cinv %*% (C-S) %*% Cinv %*% J %*% I.Ainv ", 1);


				MatrixMulti(grad_P, m, m, P, m, m, C0);
				if(SEM_DEBUG) printMatrix(C0, m, m, "grad.P %*% P", 1);

				MatrixMultiTransB(C0, m, m,invA, m, m, grad_A );
				if(SEM_DEBUG) printMatrix(grad_A, m, m, "grad.A = grad.P %*% P %*% t(I.Ainv)", 1);

				//The following code will produce gradients based on grad_A and grad_P.
				//This is the implentation of "tapply" in R.
				double *A_grad, *P_grad;
				int nA;
				int nP;
				A_grad = new double[n];
				P_grad = new double[n];
				memset(A_grad, 0, n*sizeof(double));
				memset(P_grad, 0, n*sizeof(double));

				double *grad_Au,  *grad_Pu;
				nA = length(model->arrows_1_free)/2;
				nP = length(model->arrows_2_free)/2;
				grad_Au = new double[nA];
				grad_Pu = new double[nP];
				int *tA = new int[max(max(nA*2, nP*2), max(length(model->unique_free_1), length(model->unique_free_2)))];

				Memcpy(tA, INTEGER(AS_INTEGER(model->arrows_1_free)), length(model->arrows_1_free));
				for(int i=0;i<nA;++i){
						int ir = tA[i]-1;  //fortran to C
						int il = tA[nA+i]-1;
						grad_Au[i] = grad_A[ir+il*m];  //column-wise.
				}
				for(int i=0;i<nA;i++) {
						A_grad[model->arrows_1_seq[i]-1] += grad_Au[i];
				}

				Memcpy(tA, INTEGER(AS_INTEGER(model->arrows_2_free)), nP*2);
				for(int i=0;i<nP;++i){
						int ir = tA[i]-1;
						int il = tA[nP+i]-1;
						grad_Pu[i] = grad_P[ir+il*m];  //column-wise.
				}
				for(int i=0;i<nP;i++) {
						P_grad[model->arrows_2_seq[i]-1] += grad_Pu[i];
				}

				nA=length(model->unique_free_1);
				Memcpy(tA, INTEGER(AS_INTEGER(model->unique_free_1)), nA);

				for(int i=0;i<nA;++i) g[tA[i]-1] = A_grad[tA[i]-1];

				nP=length(model->unique_free_2);
				Memcpy(tA, INTEGER(AS_INTEGER(model->unique_free_2)), nP);
				for(int i=0;i<nP;++i) g[tA[i]-1] = P_grad[tA[i]-1];

				delete[] tA;
				delete[] grad_Pu;
				delete[] grad_Au;
				delete[] P_grad;
				delete[] A_grad;
				delete[] grad_A;
				delete[] grad_P;
		}


		delete[] invA ;
		delete[] C0;
		delete[] invC;
		delete[] ImA;

		if(SEM_DEBUG) Rprintf("Exit from objectiveML.\n");

		return;
}

void msem_objectiveML(int n, const double x[], double *f, double *g, double *h,  double *A,  double *P, double *C, double *ff, msem_function_info *m_state)
{
		R_CheckUserInterrupt();

		msem_model_info *m_model = m_state->model;
		function_info *state = new function_info; //(function_info *)R_alloc(1, sizeof(function_info));
		state->have_gradient = m_state->have_gradient;
		state->have_hessian = m_state->have_hessian;

		int G = m_model->G;
		int i;

		int indAP=0,  indC=0;

		*f = 0.0;
		if(state->have_gradient) memset(g, 0, n*sizeof(double));

		int sumN=0;
		double *grad = new double[n];
		int maxmn = 0;
		int maxmni;
		for(i = 0;i < G; ++i)
		{
				sumN += INTEGER(AS_INTEGER(m_model->N))[i];
				maxmni = (m_model->gmodel[i].n > m_model->gmodel[i].m ? m_model->gmodel[i].n:m_model->gmodel[i].m);
				maxmn = (maxmni > maxmn ? maxmni : maxmn);
		}

		double *C0 = new double[maxmn*maxmn];

		for(i = 0; i < G; ++i)
		{

				state->model = &m_model->gmodel[i];

				memset(grad, 0, n*sizeof(double));
				memset(C0, 0, maxmn*maxmn);
				objectiveML(n,  x, &ff[i], grad, h,  &A[indAP],  &P[indAP], C0, state);
				Memcpy(&C[indC], C0, state->model->n*state->model->n);
				indAP += state->model->m*state->model->m;   //update the index for A, P, C
				indC += state->model->n*state->model->n;

				*f += (state->model->N-(1-state->model->raw))*ff[i];

				if(state->have_gradient)
				{
						double alpha = (state->model->N-(1-state->model->raw))/(sumN-(1.0-state->model->raw)*G);
						int incx = 1;
						//F77_NAME(daxpy)(const int *n,  const double *alpha, 
						//								const double *dx,  const int *incx, 
						//								double *dy,  const int *incy); y = ax+y
						F77_CALL(daxpy)(&n,&alpha, grad,&incx, g, &incx); //grad.all = grad.all+((N[g]-!raw)/(sum(N)-(!raw)*G))*grad
				}

		}

		*f = *f/(sumN-(1-m_model->raw)*G);

		delete[] C0;
		delete[] grad;
		delete state;
		return;
}
// this function will compute  GLS.
void objectiveGLS(int n, const double x[], double *f, double *g, double *h,  double *A,  double *P, double *C, function_info *state)
{
		R_CheckUserInterrupt();

		model_info *model = state->model;

		int m = model->m;
		int modeln = model->n;
		int maxmn = (m > modeln ? m : modeln);


		double *ImA = new double[m*m];

		generate_AP(n, x, A, P, ImA, model);
		if(SEM_DEBUG) {
				printMatrix(A, m, m, "Matrix A", 1);
				printMatrix(ImA, m, m, "Matrix (I-A)", 1);
				printMatrix(P, m, m, "Matrix P", 1);
		}


		double *invA,  *invC;
		double *C0;  //for matrix multiplication.


		invC = new double[maxmn*maxmn]; //After the compuation,  C is n-by-n.
		C0 = new double[maxmn*maxmn];  //After the compuation,  C is n-by-n.
		memset(C, 0, maxmn*maxmn*sizeof(double));
		memset(C0, 0, maxmn*maxmn*sizeof(double));
		memset(invC, 0, maxmn*maxmn*sizeof(double));

		///Please be careful that R uses column-wise to store matrix.
		invA = new double[m*m];

		Memcpy(invA, ImA, m*m);
		MatrixInverse(invA, m);
		if(SEM_DEBUG) printMatrix(invA, m, m, "invA", 1);

		if(SEM_DEBUG) {
				MatrixMulti(ImA, m, m, invA, m, m, C0);
				printMatrix(C0, m, m, "A %*% invA", 1);
		}

		MatrixMulti(REAL(model->J), nrows(model->J), ncols(model->J), invA, m, m, C0);
		if(SEM_DEBUG) printMatrix(C0, modeln, m, "J %*% I.Ainv", 1);

		MatrixMulti(C0, modeln, m, P, m, m, C);
		if(SEM_DEBUG) printMatrix(C, modeln, m, "J %*% I.Ainv %*% P", 1);

		MatrixMultiTransB(C, modeln, m, invA, m, m, C0);
		if(SEM_DEBUG) printMatrix(C0, modeln, m, "J %*% I.Ainv %*% P %*% t(I.Ainv)", 1);

		MatrixMultiTransB(C0, modeln, m, REAL(model->J),  modeln, m,  C);
		if(SEM_DEBUG) printMatrix(C, modeln, modeln, "J %*% I.Ainv %*% P %*% t(I, Ainv) %*% t(J)", 1);


		double *grad_P;

		grad_P = new double[maxmn*maxmn];  //After the compuation,  grad_P is m-by-m.

		double alpha = -1.0;
		int incx = 1;
		int mm = modeln*modeln;

		Memcpy(invC, C, modeln*modeln); //y, we use daxpy,  y=ax+y
		Memcpy(C0, REAL(model->S), modeln*modeln); //y, we use daxpy,  y=ax+y
		//F77_NAME(daxpy)(const int *n,  const double *alpha, 
		//								const double *dx,  const int *incx, 
		//								double *dy,  const int *incy);
		F77_CALL(daxpy)(&mm,&alpha, invC,&incx, C0, &incx); //C0 = C0 - grad_A ==>  grad_A = S - C
		if(SEM_DEBUG) printMatrix(C0, modeln, modeln, "(S-C)", 1);

		Memcpy(grad_P, REAL(AS_NUMERIC(model->invS)), modeln*modeln);

		MatrixMulti(grad_P, modeln, modeln, C0, modeln, modeln, invC);  /* SS (grad_A) = invS * (S-C) */

		MatrixMulti(invC, modeln, modeln, invC, modeln, modeln, C0);  /* C0 = SS %*% SS */

		*f = 0.5 * MatrixTrace(C0, modeln, modeln);

		delete[] grad_P;
		delete[] invA;
		delete[] C0;
		delete[] invC;
		delete[] ImA;

		return;
}
//
// this function will compute  GLS.
void msem_objectiveGLS(int n, const double x[], double *f, double *g, double *h,  double *A,  double *P, double *C, double *ff, msem_function_info *m_state)
{
		R_CheckUserInterrupt();

		msem_model_info *m_model = m_state->model;
		function_info *state = new function_info; //(function_info *)R_alloc(1, sizeof(function_info));
		state->have_gradient = m_state->have_gradient;
		state->have_hessian = m_state->have_hessian;

		int G = m_model->G;
		int i;

		int indAP=0,  indC=0;

		*f = 0.0;
		if(state->have_gradient) memset(g, 0, n*sizeof(double));

		int sumN=0;
		double *grad = new double[n];
		int maxmn = 0;
		for(i = 0;i < G; ++i)
		{
				sumN += INTEGER(AS_INTEGER(m_model->N))[i];
				maxmn = (m_model->gmodel[i].n > m_model->gmodel[i].m ? m_model->gmodel[i].n:m_model->gmodel[i].m);
		}
		double *C0 = new double[maxmn*maxmn];

		for(i = 0; i < G; ++i)
		{

				state->model = &m_model->gmodel[i];

				memset(grad, 0, n*sizeof(double));
				memset(C0, 0, maxmn*maxmn*sizeof(double));
				objectiveGLS(n,  x, &ff[i], grad, h,  &A[indAP],  &P[indAP], C0, state);
				Memcpy(&C[indC], C0, state->model->n*state->model->n);
				indAP += state->model->m*state->model->m;   //update the index for A, P, C
				indC += state->model->n*state->model->n;

				*f += (state->model->N-(1-state->model->raw))*ff[i];

				if(state->have_gradient)
				{
						double alpha = (state->model->N-(1-state->model->raw))/(sumN-(1.0-state->model->raw)*G);
						int incx = 1;
						//F77_NAME(daxpy)(const int *n,  const double *alpha, 
						//								const double *dx,  const int *incx, 
						//								double *dy,  const int *incy); y = ax+y
						F77_CALL(daxpy)(&n,&alpha, grad,&incx, g, &incx); //grad.all = grad.all+((N[g]-!raw)/(sum(N)-(!raw)*G))*grad
				}

		}

		*f = *f/(sumN-(1-m_model->raw)*G);

		//				Rprintf("Number of Evaluations: %d [%f]\n", m_state->n_eval, *f);

		delete[] C0;
		delete[] grad;
		delete state;

		return;
}

// this function will deal with missing data.
void objectiveFIML(int n, const double x[], double *f, double *g, double *h,  double *A,  double *P, double *C, function_info *state)
{
		R_CheckUserInterrupt();

		model_info *model = state->model;

		int m = model->m;
		int modeln = model->n;
		int maxmn = (m > modeln ? m : modeln);

		double *ImA = new double[m*m];

		generate_AP(n, x, A, P, ImA, model);
		if(SEM_DEBUG) {
				printMatrix(A, m, m, "Matrix A", 1);
				printMatrix(ImA, m, m, "Matrix (I-A)", 1);
				printMatrix(P, m, m, "Matrix P", 1);
		}


		double *invA;
		double *C0;  //for matrix multiplication.
		double *JAinv;


		C0 = new double[maxmn*maxmn];  //After the compuation,  C is n-by-n.
		JAinv = new double[modeln*m];  //for gradient.n-by-m
		memset(C, 0, maxmn*maxmn*sizeof(double));
		memset(C0, 0, maxmn*maxmn*sizeof(double));

		///Please be careful that R uses column-wise to store matrix.
		invA = new double[m*m];

		Memcpy(invA, ImA, m*m);
		MatrixInverse(invA, m);
		if(SEM_DEBUG) 
				printMatrix(invA, m, m, "invA", 1);

		if(SEM_DEBUG) {
				MatrixMulti(ImA, m, m, invA, m, m, C0);
				printMatrix(C0, m, m, "A %*% invA", 1);
		}

		MatrixMulti(REAL(model->J), nrows(model->J), ncols(model->J), invA, m, m, C0);
		Memcpy(JAinv, C0, modeln*m);  // for gradient
		if(SEM_DEBUG) 
				printMatrix(C0, modeln, m, "J %*% I.Ainv", 1);

		MatrixMulti(C0, modeln, m, P, m, m, C);
		if(SEM_DEBUG) 
				printMatrix(C, modeln, m, "J %*% I.Ainv %*% P", 1);

		MatrixMultiTransB(C, modeln, m, invA, m, m, C0);
		if(SEM_DEBUG) 
				printMatrix(C0, modeln, m, "J %*% I.Ainv %*% P %*% t(I.Ainv)", 1);

		MatrixMultiTransB(C0, modeln, m, REAL(model->J),  modeln, m,  C);
		if(SEM_DEBUG) 
				printMatrix(C, modeln, modeln, "J %*% I.Ainv %*% P %*% t(I, Ainv) %*% t(J)", 1);

		*f = 0.0;
		int Npatterns = nrows(model->valid_data_patterns);  // number of patterns,  each row represents one pattern of the missing data.

		int Npattern_number = length(model->pattern_number);   // equal to nrows(model->data)
		int *pattern_number = new int[Npattern_number];

		double *dfdC = new double[modeln*modeln];  //for gradient df/dC
		memset(dfdC, 0, modeln*modeln*sizeof(double));

		for(int i = 0; i < Npatterns; ++i)
		{
				int *sel = SubMatrixRow(model->valid_data_patterns, Npatterns, ncols(model->valid_data_patterns), i);
				if(SEM_DEBUG) 
						printMatrix(sel, 1, ncols(model->valid_data_patterns), "sel", 0);

				Memcpy(pattern_number, INTEGER(AS_INTEGER(model->pattern_number)), Npattern_number);
				if(SEM_DEBUG) 
						printMatrix(pattern_number, 1, Npattern_number, "pattern.number", 0);

				for(int j = 0; j < Npattern_number; ++j)
				{
						pattern_number[j] = (pattern_number[j]==(i+1) ? 1 : 0);
				}
				if(SEM_DEBUG) 
						printMatrix(pattern_number, 1, Npattern_number, "pattern.number", 0);

				int row_subX, col_subX;
				double *X = SubMatrix(REAL(AS_NUMERIC(model->data)),  pattern_number,  sel, Npattern_number, ncols(model->data),  row_subX, col_subX);
				if(SEM_DEBUG) 
						printMatrix(X, row_subX, col_subX, "X", 1);

				int row_subC, col_subC;
				double *subC = SubMatrix(C, sel, sel, modeln, modeln, row_subC, col_subC);  //row_subC should be equal to col_subC;
				if(SEM_DEBUG) 
						printMatrix(subC, row_subC, col_subC, "C[sel, sel]", 1);

				double detC = MatrixDeterminant(subC, row_subC, col_subC);

				*f += row_subX * (log2Pi + log(detC));

				MatrixInverse(subC, row_subC);   // now subC is the inverse of C[sel, sel]
				if(SEM_DEBUG) 
						printMatrix(subC, row_subC, col_subC, "Inverse(C[sel, sel])", 1);

				double *CC = new double[row_subX*col_subC];
				double *CCX = new double[max(row_subX*max(col_subX, row_subX), row_subC*max(row_subC, col_subC))];

				MatrixMulti(X, row_subX, col_subX, subC, row_subC, col_subC, CC);
				if(SEM_DEBUG) 
						printMatrix(CC, row_subX, col_subC, "X %*% inv(C[sel, sel])", 1);

				MatrixMultiTransB(CC, row_subX, col_subC, X, row_subX, col_subX, CCX);
				if(SEM_DEBUG) 
						printMatrix(CCX, row_subX, row_subX, "X %*% inv(C[sel, sel]) %*% X^T", 1);

				*f += MatrixTrace(CCX, row_subX, row_subX);

				// gradient

				if(state->have_gradient)
				{
						int ndfdCi = row_subC*col_subC;
						double *dfdCi = new double[ndfdCi];
						memset(dfdCi, 0, ndfdCi*sizeof(double));

						double *subC_T = new double[ndfdCi];   //transpose (C[sel, sel])^(-T)
						Memcpy(subC_T, subC, ndfdCi);
						MatrixTranspose(subC_T, col_subC, row_subC); 
						if(SEM_DEBUG) 
								printMatrix(subC_T, col_subC, row_subC, "(inv(C[sel, sel]))^T", 1);

						int incx = 1;
						//F77_NAME(daxpy)(const int *n,  const double *alpha, 
						//								const double *dx,  const int *incx, 
						//								double *dy,  const int *incy); y = ax+y
						double alpha = static_cast<double>(row_subX);
						F77_CALL(daxpy)(&ndfdCi,&alpha, subC_T,&incx, dfdCi, &incx); 

						if(SEM_DEBUG) 
								printMatrix(dfdCi, 1, ndfdCi, "nrow(X)*t(vec(t(inv(C[sel, sel]))))", 1);

						int nX = row_subX*col_subX;
						double *X_T = new double[nX];
						Memcpy(X_T, X, nX);
						MatrixTranspose(X_T, col_subX, row_subX);   // X: the column-wise
						if(SEM_DEBUG) 
								printMatrix(X_T, col_subX, row_subX, "X^T", 1);

						double *IsubC = new double[ndfdCi];  //identity matrix 
						memset(IsubC, 0, ndfdCi*sizeof(double));
						for(int j = 0; j < row_subC; ++j) 
								IsubC[j*col_subC + j] = 1.0;  // in fact,  row_subC = col_subC.

						double *XkronI = Kronecker(X, row_subX, col_subX, IsubC, row_subC, col_subC);
						if(SEM_DEBUG) 
								printMatrix(XkronI, row_subX*row_subC, col_subX*col_subC, "X %x% diag(nrow(Cinv))", 1);

						MatrixMulti(X_T, 1, nX, XkronI, row_subX*row_subC, col_subX*col_subC, CCX);
						if(SEM_DEBUG) 
								printMatrix(CCX, 1, ndfdCi, "t(vec(t(X))) %*% (X %x% diag(nrow(Cinv)))", 1);

						double *CKronC = Kronecker(subC_T, col_subC, row_subC, subC, row_subC, col_subC );
						if(SEM_DEBUG) 
								printMatrix(CKronC, col_subC*row_subC, row_subC*col_subC, "t(Cinv) %x% Cinv", 1);

						MatrixMulti(CCX, 1, ndfdCi, CKronC, ndfdCi, ndfdCi, IsubC); 
						if(SEM_DEBUG) 
								printMatrix(IsubC, 1, ndfdCi, "t(vec(t(X))) %*% (X %x% diag(nrow(Cinv))) %*% t(Cinv) %x% Cinv", 1);

						alpha = -1.0;

						F77_CALL(daxpy)(&ndfdCi,&alpha, IsubC, &incx, dfdCi, &incx); //
						if(SEM_DEBUG) 
								printMatrix(dfdCi, 1, ndfdCi, "dfdCi", 1);

						double *dfdCiExtend = ExtendMatrix(dfdCi, row_subC, col_subC, sel, ncols(model->valid_data_patterns));  //modeln
						if(SEM_DEBUG) 
								printMatrix(dfdCiExtend, modeln, modeln, "dfdCiExtend", 1);

						alpha = 1.0/static_cast<double>(model->N);
						int nn = modeln*modeln;
						F77_CALL(daxpy)(&nn,&alpha, dfdCiExtend, &incx, dfdC, &incx); //
						if(SEM_DEBUG) 
								printMatrix(dfdC, 1, nn, "dfdC = dfdC + dfdCi", 1);

						delete[] dfdCi;
						delete[] subC_T;
						delete[] X_T;
						delete[] IsubC;
						delete[] XkronI;
						delete[] CKronC;
						delete[] dfdCiExtend;
				}


				delete[] sel;
				delete[] X;
				delete[] subC;
				delete[] CC;
				delete[] CCX;

		}

		*f /= static_cast<double>(model->N);

		if(state->have_gradient)
		{
				double *dCdP = Kronecker(JAinv, modeln, m, JAinv, modeln, m);
				if(SEM_DEBUG) 
						printMatrix(dCdP, modeln*modeln, m*m, "dCdP = (J %*% I.Ainv) %x% (J %*% I.Ainv)", 1);

				memset(ImA, 0, m*m*sizeof(double));
				for(int i=0; i < m; ++i) ImA[i*m+i] = 1.0;
				double *B = Kronecker(ImA, m, m, REAL(model->J), modeln, m);
				if(SEM_DEBUG) 
						printMatrix(B, modeln*m, m*m, "(diag(nrow(A)) %x% J)", 1);

				double *tinvA = new double[m*m];
				Memcpy(tinvA, invA, m*m);
				MatrixTranspose(tinvA, m, m);
				double *kinvA = Kronecker(tinvA, m, m, invA, m, m);
				if(SEM_DEBUG) 
						printMatrix(kinvA, m*m, m*m, "(t(I.Ainv) %x% I.Ainv)", 1);

				double *dBdA = new double[modeln*m*m*m]; 
				MatrixMulti(B, m*modeln, m*m, kinvA, m*m, m*m, dBdA);
				if(SEM_DEBUG) 
						printMatrix(dBdA, m*modeln, m*m, "(diag(nrow(A)) %x% J) %*% (t(I.Ainv) %x% I.Ainv)", 1);

				double *Tmn = CommutationMatrix(m, modeln);  //column-wise
				if(SEM_DEBUG) 
						printMatrix(Tmn, modeln*m, modeln*m, "Tmn", 1);   //zhenghua,  check

				double *dCdA1 = new double[modeln*m];
				MatrixMultiTransB(JAinv, modeln, m, P, m, m, dCdA1);
				if(SEM_DEBUG) 
						printMatrix(dCdA1, modeln, m, "(B %*% t(P))", 1);

				double *ImB = new double[modeln*modeln];
				memset(ImB, 0, modeln*modeln*sizeof(double));
				for(int i=0; i<modeln; ++i)
						ImB[i*modeln+i] = 1.0;
				double *dCdA2 = Kronecker(dCdA1, modeln, m, ImB, modeln, modeln);
				if(SEM_DEBUG) 
						printMatrix(dCdA2, modeln*modeln, modeln*m, "(B %*% t(P)) %x% diag(nrow(B))", 1);

				double *dCdA30 = Kronecker(ImB, modeln, modeln, JAinv, modeln, m);
				if(SEM_DEBUG) 
						printMatrix(dCdA30, modeln*modeln, modeln*m, "(diag(nrow(B)) %x% B)", 1);

				double *dCdA31 = new double[modeln*modeln*modeln*m];
				MatrixMulti(dCdA30, modeln*modeln, modeln*m, Tmn, modeln*m, modeln*m, dCdA31);
				if(SEM_DEBUG) 
						printMatrix(dCdA31, modeln*modeln, modeln*m, "(diag(nrow(B)) %x% B) %*% Tmn", 1);

				double *dCdA32 = Kronecker(P, m, m, ImB, modeln, modeln);
				if(SEM_DEBUG) 
						printMatrix(dCdA32, modeln*m, modeln*m, "P %x% diag(nrow(B))", 1);

				double *dCdA33 = new double[modeln*modeln*modeln*m];
				MatrixMulti(dCdA31, modeln*modeln, modeln*m, dCdA32, modeln*m, modeln*m, dCdA33);
				if(SEM_DEBUG) 
						printMatrix(dCdA33, modeln*modeln, modeln*m, "(diag(nrow(B)) %x% B) %*% Tmn %*% (P %x% diag(nrow(B))) ", 1);

				double alpha = 1.0;
				int mn = modeln*modeln*modeln*m;
				int incx = 1;
				F77_CALL(daxpy)(&mn,&alpha,dCdA2 , &incx, dCdA33, &incx); //
				if(SEM_DEBUG) 
						printMatrix(dCdA33, modeln*modeln, modeln*m, "(B %*% t(P)) %x% diag(nrow(B)) + (diag(nrow(B)) %x% B) %*% Tmn %*% (P %x% diag(nrow(B))) ", 1);


				double *dCdA = new double[modeln*modeln*m*m];
				MatrixMulti(dCdA33, modeln*modeln, modeln*m, dBdA, m*modeln, m*m, dCdA);
				if(SEM_DEBUG) 
						printMatrix(dCdA, modeln*modeln, m*m, "dCdA", 1);

				double *dfdP0 = new double[m*m];
				MatrixMulti(dfdC, 1, modeln*modeln, dCdP, modeln*modeln, m*m, dfdP0);
				if(SEM_DEBUG) 
						printMatrix(dfdP0, 1, m*m, "dfdP0", 1);

				double *vcorr = vec(REAL(model->correct), m, m);
				if(SEM_DEBUG) 
						printMatrix(vcorr, 1, m*m, "correct", 1);

				double *dfdP = new double[m*m];
				int mm = m*m;
				incx = 1;
				sempdot(&mm, vcorr, &incx, dfdP0, &incx, dfdP);
				if(SEM_DEBUG) 
						printMatrix(dfdP, 1, m*m, "dfdP", 1);

				double *dfdA = new double[m*m];
				MatrixMulti(dfdC, 1, modeln*modeln, dCdA, modeln*modeln, m*m, dfdA);
				if(SEM_DEBUG) 
						printMatrix(dfdA, 1, m*m, "dfdA", 1);

				//The following code will produce gradients based on grad_A and grad_P.
				//This is the implentation of "tapply" in R.
				double *A_grad, *P_grad;
				int nA;
				int nP;
				A_grad = new double[n];
				P_grad = new double[n];
				memset(A_grad, 0, n*sizeof(double));
				memset(P_grad, 0, n*sizeof(double));

				double *grad_Au,  *grad_Pu;
				nA = length(model->arrows_1_free)/2;
				nP = length(model->arrows_2_free)/2;
				grad_Au = new double[nA];
				grad_Pu = new double[nP];
				int *tA = new int[max(max(nA*2, nP*2), max(length(model->unique_free_1), length(model->unique_free_2)))];

				Memcpy(tA, INTEGER(AS_INTEGER(model->arrows_1_free)), length(model->arrows_1_free));
				for(int i=0;i<nA;++i){
						int ir = tA[i]-1;  //fortran to C
						int il = tA[nA+i]-1;
						grad_Au[i] = dfdA[ir+il*m];  //column-wise.
				}
				for(int i=0;i<nA;i++) {
						A_grad[model->arrows_1_seq[i]-1] += grad_Au[i];
				}

				Memcpy(tA, INTEGER(AS_INTEGER(model->arrows_2_free)), nP*2);
				for(int i=0;i<nP;++i){
						int ir = tA[i]-1;
						int il = tA[nP+i]-1;
						grad_Pu[i] = dfdP[ir+il*m];  //column-wise.
				}
				for(int i=0;i<nP;i++) {
						P_grad[model->arrows_2_seq[i]-1] += grad_Pu[i];
				}

				nA=length(model->unique_free_1);
				Memcpy(tA, INTEGER(AS_INTEGER(model->unique_free_1)), nA);

				for(int i=0;i<nA;++i) g[tA[i]-1] = A_grad[tA[i]-1];

				nP=length(model->unique_free_2);
				Memcpy(tA, INTEGER(AS_INTEGER(model->unique_free_2)), nP);
				for(int i=0;i<nP;++i) g[tA[i]-1] = P_grad[tA[i]-1];

				delete[] dCdP;
				delete[] B;
				delete[] tinvA;
				delete[] kinvA;
				delete[] dBdA;
				delete[] dCdA1;
				delete[] ImB;
				delete[] dCdA2;
				delete[] dCdA30;
				delete[] dCdA31;
				delete[] dCdA32;
				delete[] dCdA33;
				delete[] dCdA;
				delete[] dfdP0;
				delete[] vcorr;
				delete[] dfdP;
				delete[] dfdA;

				delete[] A_grad;
				delete[] P_grad;
				delete[] grad_Au;
				delete[] grad_Pu;
				delete[] tA;
		}

		delete[] ImA;
		delete[] C0;
		delete[] JAinv;
		delete[] invA;

		delete[] pattern_number;
		delete[] dfdC;


		return;
}

void objectivelogLik(int n, const double x[], double *f, double *g, double *h,  double *A,  double *P, double *C, function_info *state)
{
		R_CheckUserInterrupt();

		*A = csem_NaN; 
		*P = csem_NaN;
		*C = csem_NaN;

		sem_object *semObject = state->model->semObject;
		int ncolData = ncols(semObject->data);

		double *C0 = new double[ncolData*ncolData];

		memset(C0, 0, ncolData*ncolData*sizeof(double));

		int posn_intercept = semObject->posn_intercept - 1;  // fortran to C

		C0[posn_intercept*n+posn_intercept] = 1.0;


		int ncolTri = ncols(semObject->tri);
		int nrowTri = nrows(semObject->tri);
		int *tri = new int[length(semObject->tri)];
		Memcpy(tri, INTEGER(AS_INTEGER(semObject->tri)), length(semObject->tri));

		//printMatrix(tri, nrowTri, ncolTri, "tri", 1);

		for(int i = 0; i < ncolTri; ++i)
		{
				for(int j = 0; j < nrowTri; ++j)
						if(tri[i*nrowTri + j])
								C0[i*nrowTri + j] = *x++;
						else
								C0[i*nrowTri + j] = C0[j*nrowTri + i];

		}  // C <- C + t(C) - diag(diag(C))

		//printMatrix(C0, ncolData, ncolData, "C", 1);


		*f = 0.0;
		int Npatterns = nrows(semObject->valid_data_patterns);  // number of patterns,  each row represents one pattern of the missing data.

		int Npattern_number = length(semObject->pattern_number);   // equal to nrows(model->data)
		int *pattern_number = new int[Npattern_number];

		//	double *dfdC = new double[ncolTri*nrowTri];  //for gradient df/dC
		//	memset(dfdC, 0, ncolTri*nrowTri*sizeof(double));

		for(int i = 0; i < Npatterns; ++i)
		{
				int *sel = SubMatrixRow(semObject->valid_data_patterns, Npatterns, ncols(semObject->valid_data_patterns), i);
				if(SEM_DEBUG) 
						printMatrix(sel, 1, ncols(semObject->valid_data_patterns), "sel", 0);

				Memcpy(pattern_number, INTEGER(AS_INTEGER(semObject->pattern_number)), Npattern_number);
				if(SEM_DEBUG) 
						printMatrix(pattern_number, 1, Npattern_number, "pattern.number", 0);

				for(int j = 0; j < Npattern_number; ++j)
				{
						pattern_number[j] = (pattern_number[j]==(i+1) ? 1 : 0);
				}
				if(SEM_DEBUG) 
						printMatrix(pattern_number, 1, Npattern_number, "pattern.number", 0);

				int row_subX, col_subX;
				double *X = SubMatrix(REAL(AS_NUMERIC(semObject->data)),  pattern_number,  sel, Npattern_number, ncols(semObject->data),  row_subX, col_subX);
				if(SEM_DEBUG) 
						printMatrix(X, row_subX, col_subX, "X", 1);

				int row_subC, col_subC;
				double *subC = SubMatrix(C0, sel, sel, nrowTri, nrowTri, row_subC, col_subC);  //row_subC should be equal to col_subC;
				if(SEM_DEBUG) 
						printMatrix(subC, row_subC, col_subC, "C[sel, sel]", 1);

				double detC = MatrixDeterminant(subC, row_subC, col_subC);

				*f += row_subX * (log2Pi + log(detC));

				MatrixInverse(subC, row_subC);   // now subC is the inverse of C[sel, sel]
				if(SEM_DEBUG) 
						printMatrix(subC, row_subC, col_subC, "Inverse(C[sel, sel])", 1);

				double *CC = new double[row_subX*col_subC];
				double *CCX = new double[max(row_subX*max(col_subX, row_subX), row_subC*max(row_subC, col_subC))];

				MatrixMulti(X, row_subX, col_subX, subC, row_subC, col_subC, CC);
				if(SEM_DEBUG) 
						printMatrix(CC, row_subX, col_subC, "X %*% inv(C[sel, sel])", 1);

				MatrixMultiTransB(CC, row_subX, col_subC, X, row_subX, col_subX, CCX);
				if(SEM_DEBUG) 
						printMatrix(CCX, row_subX, row_subX, "X %*% inv(C[sel, sel]) %*% X^T", 1);

				*f += MatrixTrace(CCX, row_subX, row_subX);


				delete[] sel;
				delete[] X;
				delete[] subC;
				delete[] CC;
				delete[] CCX;

		}


		delete[] pattern_number;
		delete[] C0;
		delete[] tri;
		return;
}
//
// this function will compute  GLS.
void msem_objectiveFIML(int n, const double x[], double *f, double *g, double *h,  double *A,  double *P, double *C, double *ff, msem_function_info *m_state)
{
		R_CheckUserInterrupt();

		msem_model_info *m_model = m_state->model;
		function_info *state = new function_info; //(function_info *)R_alloc(1, sizeof(function_info));
		state->have_gradient = m_state->have_gradient;
		state->have_hessian = m_state->have_hessian;

		int G = m_model->G;
		int i;

		int indAP=0,  indC=0;

		*f = 0.0;
		if(state->have_gradient) memset(g, 0, n*sizeof(double));

		int sumN=0;
		double *grad = new double[n];
		int maxmn = 0;
		for(i = 0;i < G; ++i)
		{
				sumN += INTEGER(AS_INTEGER(m_model->N))[i];
				maxmn = (m_model->gmodel[i].n > m_model->gmodel[i].m ? m_model->gmodel[i].n:m_model->gmodel[i].m);
		}
		double *C0 = new double[maxmn*maxmn];

		for(i = 0; i < G; ++i)
		{

				state->model = &m_model->gmodel[i];

				memset(grad, 0, n*sizeof(double));
				memset(C0, 0, maxmn*maxmn*sizeof(double));
				objectiveFIML(n,  x, &ff[i], grad, h,  &A[indAP],  &P[indAP], C0, state);
				Memcpy(&C[indC], C0, state->model->n*state->model->n);
				indAP += state->model->m*state->model->m;   //update the index for A, P, C
				indC += state->model->n*state->model->n;

				*f += (state->model->N-(1-state->model->raw))*ff[i];

				if(state->have_gradient)
				{
						double alpha = (state->model->N-(1-state->model->raw))/(sumN-(1.0-state->model->raw)*G);
						int incx = 1;
						//F77_NAME(daxpy)(const int *n,  const double *alpha, 
						//								const double *dx,  const int *incx, 
						//								double *dy,  const int *incy); y = ax+y
						F77_CALL(daxpy)(&n,&alpha, grad,&incx, g, &incx); //grad.all = grad.all+((N[g]-!raw)/(sum(N)-(!raw)*G))*grad
				}

		}

		*f = *f/(sumN-(1-m_model->raw)*G);

		//				Rprintf("Number of Evaluations: %d [%f]\n", m_state->n_eval, *f);

		delete[] C0;
		delete[] grad;
		delete state;

		return;
}


/*
	 "start" = start, 
	 "opts" = opts, 
	 "S" = model$S, 
	 "logdetS" = model$logdetS, 
	 "invS" = model$invS, 
	 "N" = model$N, 
	 "m" = model$m, 
	 "n" = model$n, 
	 "t" = model$t, 
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
//

//column_wise: 1

SEXP csemSolve( SEXP args )
{
		R_CheckUserInterrupt();

		theenv = getListElement(args, "csem.environment");
		thefun = getListElement(args, "print.f");  // for print a Real vector

		csem_NaN = std::numeric_limits<double>::quiet_NaN();

		SEXP solution;
		int num_prot = 0;

		PROTECT(solution=args);

		//showArgs1(args);
		// Define objective functions and their properties.
		const int num_objs = 5;
		const string name_objs[num_objs] =
		{
				"objectiveML", 
				"objectiveGLS", 
				"objectiveFIML", 
				"objectivelogLik", 
				"test_objective"
		};
		const myfcn_p objectiveFun[num_objs]=
		{
				(myfcn_p) objectiveML,   //objective, gradient (iagflg[0]=1),  no hessian (iahflg[0]=0)
				(myfcn_p) objectiveGLS,  //objective, no gradient (iagflg[1]=0), no hessian (iahflg[1]=0) 
				(myfcn_p) objectiveFIML,  //objective, gradient (iagflg[2]=1), no hessian (iahflg[1]=0) 
				(myfcn_p) objectivelogLik,  //objective,  no gradient (iaglfg[1]=0i), no hessian (iahflg[1]=0)
				(myfcn_p) test_objective  //objective, gradient, hessian
		};
		const int iagflg[num_objs]={1, 0, 1, 0, 1};  //gradients
		const int iahflg[num_objs]={0, 0, 0, 0, 1};  //hessian
		int gradient = 1;                 //0: we don't use gradient to optimize.

		int obj_ind = 0;  //default objective function.

		SEXP st;
		st= getListElement(args, "objective");
		for (int i = 0; i < num_objs; ++i) {
				if(name_objs[i].compare(CHAR(STRING_ELT(st, 0))) == 0) {
						obj_ind = i;
						break;
				}
		}
		if(SEM_DEBUG) printSEXP(st, "Objective Function");

		st = getListElement(args, "gradient");
		gradient = INTEGER(st)[0];

		int optimize;   //0: only compute the objective function,  gradients and  Hessian if it is provided.
		st = getListElement(args, "opt.flg");
		optimize = INTEGER(st)[0];

		model_info *model = NULL;

		model = new model_info;


		if(obj_ind < 3 )
		{
				// model 

				st = getListElement(args, "logdetS");
				model->logdetS = REAL(st)[0];
				st = getListElement(args, "N");
				model->N = INTEGER(st)[0];
				st = getListElement(args, "t");
				model->t = INTEGER(st)[0];
				st = getListElement(args, "n");
				model->n = INTEGER(st)[0];
				st = getListElement(args, "m");
				model->m = INTEGER(st)[0];
				PROTECT(model->ram = getListElement(args, "ram"));
				PROTECT(model->sel_free = getListElement(args, "sel.free"));
				PROTECT(model->arrows_1 = getListElement(args, "arrows.1"));
				PROTECT(model->arrows_1_free = getListElement(args, "arrows.1.free"));
				PROTECT(model->one_head = getListElement(args, "one.head"));
				PROTECT(model->arrows_2t = getListElement(args, "arrows.2t"));
				PROTECT(model->arrows_2 = getListElement(args, "arrows.2"));
				PROTECT(model->arrows_2_free = getListElement(args, "arrows.2.free"));
				PROTECT(model->unique_free_1 = getListElement(args, "unique.free.1"));
				PROTECT(model->unique_free_2 = getListElement(args, "unique.free.2"));
				PROTECT(model->param_names = getListElement(args, "param.names"));
				PROTECT(model->var_names = getListElement(args, "var.names"));
				PROTECT(model->one_free = getListElement(args, "one.free"));
				PROTECT(model->two_free = getListElement(args, "two.free"));

				PROTECT(model->S = getListElement(args, "S"));
				PROTECT(model->invS = getListElement(args, "invS"));
				PROTECT(model->fixed = getListElement(args, "fixed"));
				PROTECT(model->J = getListElement(args, "J"));
				PROTECT(model->correct = getListElement(args, "correct"));
				num_prot += 19;

				PROTECT(model->data = getListElement(args, "data"));
				PROTECT(model->pattern_number = getListElement(args, "pattern.number"));
				PROTECT(model->valid_data_patterns = getListElement(args, "valid.data.patterns"));
				num_prot += 3;



				st = getListElement(args, "raw");
				model->raw = INTEGER(st)[0];

				st = getListElement(args, "arrows.1.seq");
				model->arrows_1_seq = (int *)R_alloc(length(st), sizeof(int)); 
				Memcpy(model->arrows_1_seq, INTEGER(AS_INTEGER(st)), length(st));

				st = getListElement(args, "arrows.2.seq");
				model->arrows_2_seq = (int *)R_alloc(length(st), sizeof(int)); 
				Memcpy(model->arrows_2_seq, INTEGER(AS_INTEGER(st)), length(st));

				//Print if debug
				if(SEM_DEBUG){
						printSEXP(model->S, "\nMatrix S");
						printSEXP(model->invS, "\nMatrix invS");
						printSEXP(model->J, "\nMatrix J");
						printSEXP(model->fixed, "\nVector fixed");
						printSEXP(model->correct, "\nMatrix correct");
				}

				//initial values x0 and typsize
				double *x0 = new double[model->t];
				double *typsiz = new double[model->t];
				//typsiz = (double *)R_alloc(model->t, sizeof(double));
				int ind;
				double sum = 0.0;

				SEXP sx0 = getListElement(args, "start");
				if(LENGTH(sx0) != model->t) error(("The number of variables are not consistent!\n"));

				SEXP stypsiz = getListElement(args, "typsize");
				Memcpy(typsiz, REAL(AS_NUMERIC(stypsiz)), model->t);

				//x0 = (double *)R_alloc(model->t, sizeof(double));

				for(ind=0;ind <model->t;++ind){
						R_CheckUserInterrupt();
						x0[ind]=REAL(sx0)[ind];
						sum += (x0[ind]/typsiz[ind])*(x0[ind]/typsiz[ind]);
				}

				//options for optimization

				double stepmax=1000*sqrt(sum);
				stepmax = stepmax > 1000.0? stepmax: 1000.0;
				int hessian; 
				double fscale;
				double steptol;
				double gradtol;
				int iterlim;
				int ndigit;
				int print_level;
				int check_analyticals;
				int msg=0;
				int msg_print[]={8, 0, 16};

				SEXP opts;

				PROTECT(opts = getListElement(args, "options"));
				setApplicationOptions(hessian, fscale, steptol, stepmax, iterlim, ndigit, print_level, check_analyticals,gradtol, opts );
				UNPROTECT(1);

				if(SEM_DEBUG) {
						Rprintf("hessian: [%d]\n", hessian);
						Rprintf("iterlim: [%d]\n", iterlim);
						Rprintf("ndigit: [%d]\n", ndigit);
						Rprintf("print.level: [%d]\n", print_level);
						Rprintf("check.analyticals: [%d]\n", check_analyticals);
						Rprintf("fscale: [%f]\n", fscale);
						Rprintf("steptol: [%f]\n", steptol);
						Rprintf("stepmax: [%f]\n", stepmax);
						Rprintf("gradtol: [%f]\n", gradtol);
				}

				msg = 1+msg_print[print_level];
				if (check_analyticals==0) msg += 2 + 4;


				solution = csemnlm(x0, model->t, iagflg[obj_ind]&gradient, iahflg[obj_ind], hessian, typsiz, fscale, msg, ndigit, gradtol, 
								stepmax, steptol,  iterlim, model, (myfcn_p) objectiveFun[obj_ind], optimize);

				delete[] x0;
				delete[] typsiz;
		}
		else
		{
				//for matrix A,  P,  C
				model->m = 1;
				model->n = 1;

				sem_object *semObject = new sem_object;

				st = getListElement(args, "t");
				semObject->t = INTEGER(st)[0];

				st = getListElement(args, "posn.intercept");
				semObject->posn_intercept = INTEGER(st)[0];

				PROTECT(semObject->data = getListElement(args, "data"));
				PROTECT(semObject->pattern_number = getListElement(args, "pattern.number"));
				PROTECT(semObject->valid_data_patterns = getListElement(args, "valid.data.patterns"));
				PROTECT(semObject->tri = getListElement(args, "tri"));
				num_prot += 4;

				model->semObject = semObject;

				//printSEXP(semObject->valid_data_patterns, "valid.data.patterns");
				//printSEXP(semObject->pattern_number, "pattern.number");

				//initial values x0 and typsize
				double *x0 = new double[semObject->t];
				double *typsiz = new double[semObject->t];

				int ind;
				double sum = 0.0;

				SEXP sx0 = getListElement(args, "start");
				if(LENGTH(sx0) != semObject->t) error(("The number of variables are not consistent!\n"));

				SEXP stypsiz = getListElement(args, "typsize");
				Memcpy(typsiz, REAL(AS_NUMERIC(stypsiz)), semObject->t);


				for(ind=0;ind <semObject->t;++ind){
						R_CheckUserInterrupt();
						x0[ind]=REAL(sx0)[ind];
						sum += (x0[ind]/typsiz[ind])*(x0[ind]/typsiz[ind]);
				}

				//options for optimization

				double stepmax=1000*sqrt(sum);
				stepmax = stepmax > 1000.0? stepmax: 1000.0;
				int hessian; 
				double fscale;
				double steptol;
				double gradtol;
				int iterlim;
				int ndigit;
				int print_level;
				int check_analyticals;
				int msg=0;
				int msg_print[]={8, 0, 16};

				SEXP opts;

				PROTECT(opts = getListElement(args, "options"));
				setApplicationOptions(hessian, fscale, steptol, stepmax, iterlim, ndigit, print_level, check_analyticals,gradtol, opts );
				UNPROTECT(1);

				msg = 1+msg_print[print_level];
				if (check_analyticals==0) msg += 2 + 4;

				solution = csemnlm(x0, semObject->t, iagflg[obj_ind]&gradient, iahflg[obj_ind], hessian, typsiz, fscale, msg, ndigit, gradtol, 
								stepmax, steptol,  iterlim, model, (myfcn_p) objectiveFun[obj_ind], optimize);

				delete semObject;
				delete[] x0;
				delete[] typsiz;
		}

		UNPROTECT(num_prot);

		delete model;

		return(solution);

}

SEXP cmsemSolve( SEXP args )
{
		R_CheckUserInterrupt();

		theenv = getListElement(args, "csem.environment");
		thefun = getListElement(args, "print.f");  // for print a Real vector

		csem_NaN = std::numeric_limits<double>::quiet_NaN();

		SEXP solution;
		int num_prot = 0;

		PROTECT(solution=args);

		//	showArgs1(args);

		// Define objective functions and their properties.
		const int num_objs = 4;
		const string name_objs[num_objs] =
		{
				"objectiveML", 
				"objectiveGLS", 
				"objectiveFIML", 
				"test_objective"
		};
		const msem_fcn_p objectiveFun[num_objs]=
		{
				(msem_fcn_p) msem_objectiveML,   //objective, gradient (iagflg[0]=1),  no hessian (iahflg[0]=0)
				(msem_fcn_p) msem_objectiveGLS,  //objective, no gradient (iagflg[1]=0), no hessian (iahflg[1]=0) 
				(msem_fcn_p) msem_objectiveFIML,  //objective, no gradient (iagflg[1]=0), no hessian (iahflg[1]=0) 
				(msem_fcn_p) msem_test_objective  //objective, gradient, hessian
		};
		const int iagflg[num_objs]={1, 0, 0, 1};  //gradients
		const int iahflg[num_objs]={0, 0, 0, 1};  //hessian
		int gradient=1;

		int obj_ind = 0;  //default objective function.

		SEXP st;
		st= getListElement(args, "objective");
		for (int i = 0; i < num_objs; ++i) {
				if(name_objs[i].compare(CHAR(STRING_ELT(st, 0))) == 0) {
						obj_ind = i;
						break;
				}
		}
		if(SEM_DEBUG) printSEXP(st, "Objective Function");

		st = getListElement(args, "gradient");
		gradient = INTEGER(st)[0];

		int optimize;   //0: only compute the objective function,  gradients and  Hessian if it is provided.
		st = getListElement(args, "opt.flg");
		optimize = INTEGER(st)[0];

		// model 
		msem_model_info *model;

		model = new msem_model_info;


		st = getListElement(args, "G");
		model->G = INTEGER(AS_INTEGER(st))[0];

		PROTECT(model->logdetS = getListElement(args, "logdetS"));
		PROTECT(model->N = getListElement(args, "N"));
		st = getListElement(args, "t");
		model->t = INTEGER(st)[0];
		PROTECT(model->n = getListElement(args, "n"));
		PROTECT(model->m = getListElement(args, "m"));
		PROTECT(model->ram = getListElement(args, "ram"));
		PROTECT(model->sel_free = getListElement(args, "sel.free"));
		PROTECT(model->arrows_1 = getListElement(args, "arrows.1"));
		PROTECT(model->arrows_1_free = getListElement(args, "arrows.1.free"));
		PROTECT(model->one_head = getListElement(args, "one.head"));
		PROTECT(model->arrows_2t = getListElement(args, "arrows.2t"));
		PROTECT(model->arrows_2 = getListElement(args, "arrows.2"));
		PROTECT(model->arrows_2_free = getListElement(args, "arrows.2.free"));
		PROTECT(model->unique_free_1 = getListElement(args, "unique.free.1"));
		PROTECT(model->unique_free_2 = getListElement(args, "unique.free.2"));
		PROTECT(model->param_names = getListElement(args, "param.names"));
		PROTECT(model->var_names = getListElement(args, "var.names"));
		PROTECT(model->one_free = getListElement(args, "one.free"));
		PROTECT(model->two_free = getListElement(args, "two.free"));

		PROTECT(model->S = getListElement(args, "S"));
		PROTECT(model->invS = getListElement(args, "invS"));
		PROTECT(model->fixed = getListElement(args, "fixed"));
		PROTECT(model->J = getListElement(args, "J"));
		PROTECT(model->correct = getListElement(args, "correct"));
		num_prot += 23;

		if(obj_ind == 2) 
		{ // objectiveFIML
				PROTECT(model->data = getListElement(args, "data"));
				PROTECT(model->pattern_number = getListElement(args, "pattern.number"));
				PROTECT(model->valid_data_patterns = getListElement(args, "valid.data.patterns"));
				num_prot += 3;
		}

		st = getListElement(args, "raw");
		model->raw = INTEGER(st)[0];

		PROTECT(model->arrows_1_seq = getListElement(args, "arrows.1.seq"));
		PROTECT(model->arrows_2_seq = getListElement(args, "arrows.2.seq"));
		num_prot += 2;

		//produce pointer for each group's model.
		model->gmodel = new model_info[model->G];
		for(int i = 0; i < model->G; ++i)
		{
				model_info *gmodel = &model->gmodel[i];

				PROTECT(gmodel->S=getListElement(model->S, i));
				gmodel->logdetS = REAL(AS_NUMERIC(model->logdetS))[i];
				gmodel->N = INTEGER(AS_INTEGER(model->N))[i];
				gmodel->m = INTEGER(AS_INTEGER(model->m))[i];
				gmodel->n = INTEGER(AS_INTEGER(model->n))[i];
				PROTECT(gmodel->fixed = getListElement(model->fixed, i));
				PROTECT(gmodel->ram = getListElement(model->ram, i));
				PROTECT(gmodel->sel_free = getListElement(model->sel_free, i));
				gmodel->t = length(model->sel_free);
				PROTECT(gmodel->arrows_1 = getListElement(model->arrows_1, i));
				PROTECT(gmodel->arrows_1_free = getListElement(model->arrows_1_free, i));
				PROTECT(gmodel->one_head = getListElement(model->one_head, i));
				PROTECT(gmodel->arrows_2t = getListElement(model->arrows_2t, i));
				PROTECT(gmodel->arrows_2 = getListElement(model->arrows_2, i));
				PROTECT(gmodel->arrows_2_free = getListElement(model->arrows_2_free, i));
				PROTECT(gmodel->unique_free_1 = getListElement(model->unique_free_1, i));
				PROTECT(gmodel->unique_free_2 = getListElement(model->unique_free_2, i));
				PROTECT(gmodel->J = getListElement(model->J, i));
				PROTECT(gmodel->correct = getListElement(model->correct, i));
				st = getListElement(model->arrows_1_seq, i);
				gmodel->arrows_1_seq = (int *)R_alloc(length(st), sizeof(int)); 
				Memcpy(gmodel->arrows_1_seq, INTEGER(AS_INTEGER(st)), length(st));
				st = getListElement(model->arrows_2_seq, i);
				gmodel->arrows_2_seq = (int *)R_alloc(length(st), sizeof(int)); 
				Memcpy(gmodel->arrows_2_seq, INTEGER(AS_INTEGER(st)), length(st));
				gmodel->raw = model->raw;
				num_prot += 14;

				//inverse of S for GLS
				if(obj_ind == 1)   //objectiveGLS
				{
						int nrow = nrows(gmodel->S);
						int ncol = ncols(gmodel->S);
						double *invS = new double[nrow*ncol];
						Memcpy(invS, REAL(AS_NUMERIC(gmodel->S)), nrow*ncol);
						MatrixInverse(invS, nrow);
						PROTECT(gmodel->invS=generateMatrix(invS, nrow, ncol));
						num_prot++;
				}

				// for objectiveFIML
				if(obj_ind == 2) 
				{
						PROTECT(gmodel->data = getListElement(model->data, i));
						PROTECT(gmodel->valid_data_patterns = getListElement(model->valid_data_patterns, i));
						PROTECT(gmodel->pattern_number = getListElement(model->pattern_number, i));
						num_prot += 3;
				}
		}

		//Print if debug
		if(SEM_DEBUG){
				printSEXP(model->S, "\nMatrix S");
				printSEXP(model->invS, "\nMatrix invS");
				printSEXP(model->J, "\nMatrix J");
				printSEXP(model->fixed, "\nVector fixed");
				printSEXP(model->correct, "\nMatrix correct");
		}

		//initial values x0 and typsize
		double *x0 = new double[model->t];
		double *typsiz = new double[model->t];
		//typsiz = (double *)R_alloc(model->t, sizeof(double));
		int ind;
		double sum = 0.0;

		SEXP sx0 = getListElement(args, "start");
		if(LENGTH(sx0) != model->t) error(("The number of variables are not consistent!\n"));

		SEXP stypsiz = getListElement(args, "typsize");
		Memcpy(typsiz, REAL(AS_NUMERIC(stypsiz)), model->t);

		//x0 = (double *)R_alloc(model->t, sizeof(double));

		for(ind=0;ind <model->t;++ind){
				R_CheckUserInterrupt();
				x0[ind]=REAL(sx0)[ind];
				sum += (x0[ind]/typsiz[ind])*(x0[ind]/typsiz[ind]);
		}


		//options for optimization

		double stepmax=1000*sqrt(sum);
		stepmax = stepmax > 1000.0? stepmax: 1000.0;
		int hessian; 
		double fscale;
		double steptol;
		double gradtol;
		int iterlim;
		int ndigit;
		int print_level;
		int check_analyticals;
		int msg=0;
		int msg_print[]={8, 0, 16};

		SEXP opts;

		PROTECT(opts = getListElement(args, "options"));
		setApplicationOptions(hessian, fscale, steptol, stepmax, iterlim, ndigit, print_level, check_analyticals,gradtol, opts );
		UNPROTECT(1);

		if(SEM_DEBUG) {
				Rprintf("hessian: [%d]\n", hessian);
				Rprintf("iterlim: [%d]\n", iterlim);
				Rprintf("ndigit: [%d]\n", ndigit);
				Rprintf("print.level: [%d]\n", print_level);
				Rprintf("check.analyticals: [%d]\n", check_analyticals);
				Rprintf("fscale: [%f]\n", fscale);
				Rprintf("steptol: [%f]\n", steptol);
				Rprintf("stepmax: [%f]\n", stepmax);
				Rprintf("gradtol: [%f]\n", gradtol);
		}

		msg = 1+msg_print[print_level];
		if (check_analyticals==0) msg += 2 + 4;

		solution = cmsemnlm(x0, model->t, iagflg[obj_ind]&gradient, iahflg[obj_ind], hessian, typsiz, fscale, msg, ndigit, gradtol, 
						stepmax, steptol,  iterlim, model, (msem_fcn_p) objectiveFun[obj_ind], optimize);

		UNPROTECT(num_prot);
		delete[] model->gmodel;
		delete model;
		delete[] x0;
		delete[] typsiz;

		return(solution);

}

}  /* extern "C"  */
