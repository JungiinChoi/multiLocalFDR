/* created by calling tools::package_native_routine_registration_skeleton('fmlogcondens',,,FALSE) from R */
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME:
   Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void calcExactIntegralC(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void calcKernelDens(void *, void *, void *, void *, void *, void *);
extern void newtonBFGSLC(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void newtonBFGSLInitC(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void printAVXInfo();
extern void recalcParamsC(void *, void *, void *, void *, void *, void *, void *);

static const R_CMethodDef CEntries[] = {
    {"calcExactIntegralC", (DL_FUNC) &calcExactIntegralC, 10},
    {"calcKernelDens",     (DL_FUNC) &calcKernelDens,      6},
    {"newtonBFGSLC",       (DL_FUNC) &newtonBFGSLC,       18},
    {"newtonBFGSLInitC",   (DL_FUNC) &newtonBFGSLInitC,   13},
    {"printAVXInfo",       (DL_FUNC) &printAVXInfo,        0},
    {"recalcParamsC",      (DL_FUNC) &recalcParamsC,       7},
    {NULL, NULL, 0}
};

void R_init_fmlogcondens(DllInfo *dll)
{   
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
	R_forceSymbols(dll, TRUE);
}

