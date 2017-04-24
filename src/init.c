#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP cmsemSolve(SEXP);
extern SEXP csemSolve(SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"cmsemSolve", (DL_FUNC) &cmsemSolve, 1},
    {"csemSolve",  (DL_FUNC) &csemSolve,  1},
    {NULL, NULL, 0}
};

void R_init_sem(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
 
