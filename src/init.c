#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* .Call calls */
extern SEXP log_dmvncholC(SEXP, SEXP, SEXP, SEXP);
extern SEXP log_dmvnmodcholC(SEXP, SEXP, SEXP, SEXP);
extern SEXP mu_score_mvncholC(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP mu_score_mvnmodcholC(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP lamdiag_score_mvncholC(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP innov_score_mvnmodcholC(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP lambda_score_mvncholC(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP phi_score_mvnmodcholC(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef callMethods[]  = {
    {"log_dmvncholC", (DL_FUNC) &log_dmvncholC, 4},
    {"log_dmvnmodcholC", (DL_FUNC) &log_dmvnmodcholC, 4},
    {"mu_score_mvncholC", (DL_FUNC) &mu_score_mvncholC, 5},
    {"mu_score_mvnmodcholC", (DL_FUNC) &mu_score_mvnmodcholC, 5},
    {"lamdiag_score_mvncholC", (DL_FUNC) &lamdiag_score_mvncholC, 5},
    {"innov_score_mvnmodcholC", (DL_FUNC) &innov_score_mvnmodcholC, 5},
    {"lambda_score_mvncholC", (DL_FUNC) &lambda_score_mvncholC, 6},
    {"phi_score_mvnmodcholC", (DL_FUNC) &phi_score_mvnmodcholC, 6},
    {NULL, NULL, 0}
};

void R_init_bamlssMVN(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, callMethods, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}


