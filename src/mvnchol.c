#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>
#include <Rinternals.h>

/*
Function: log_dmvnchol

Input:
Y observations (matrix)
PAR parameters (matrix)
N sample size (integer)
K dimension of multivariate distribution
MJ number of columns with location parameters
DJ number of columns diagonal elements of L^{-1}
TJ number of columns element in lower triangle of L^{-1}
*/
SEXP log_dmvncholC(SEXP Y, SEXP PAR, SEXP N, SEXP K, SEXP MJ, SEXP DJ, SEXP TJ)
{
/* Pointers for obs and par */
  double *Yptr = REAL(Y);
  double *PARptr = REAL(PAR);

/* integers from input */
  int n = asInteger(N);
  int k = asInteger(K);
  int mj = asInteger(MJ);
  int dj = asInteger(DJ);
  int tj = asInteger(TJ);

/* loop counters */
  int i, j;

/* inverse of L matrix for single obs */
  SEXP Linv;
  PROTECT(Linv = allocMatrix(REALSXP, k, k));
  double *Linvptr = REAL(Linv);

/* residuals for single obs */
  SEXP ymu;
  PROTECT(ymu = allocVector(REALSXP, k));
  double *ymuptr = REAL(ymu);

/* return value: vector of log-densities */
  SEXP d;
  PROTECT(d = allocVector(REALSXP, n));
  double *dptr = REAL(d);

/* initialize three terms of log-density */
  double lpi = -0.5 * (float)k * log(2.0 * 3.14159265358979323846);
  double det = 0.0;
  double nrm = 0.0;

/* HERE comes the code!!!
   I start with returning a constant vector.
   For checking if the declaration section gives any errors!
*/
  for(i = 0; i < n; i++) {
/* compute term2 the determinante */
    det = 0.0;
    for(j = 0; j < k; j++) {
      det += PARptr[i + n * (j + mj)];
    }

    dptr[i] = lpi + det;
  }

/* return */
  UNPROTECT(3);
  return d;
}



