#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP CanonicalToRealForR(SEXP, SEXP, SEXP);
extern SEXP ComputeMLEForR(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP RealToCanonicalForR(SEXP, SEXP);
extern SEXP ReductionStepForR(SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"CanonicalToRealForR", (DL_FUNC) &CanonicalToRealForR, 3},
    {"ComputeMLEForR",      (DL_FUNC) &ComputeMLEForR,      5},
    {"RealToCanonicalForR", (DL_FUNC) &RealToCanonicalForR, 2},
    {"ReductionStepForR",   (DL_FUNC) &ReductionStepForR,   4},
    {NULL, NULL, 0}
};

void R_init_MLEcens(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
