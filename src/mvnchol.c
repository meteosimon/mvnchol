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

Output: Vector of log-densities
*/
SEXP log_dmvncholC(SEXP Y, SEXP PAR, SEXP N, SEXP K)
{
/* Pointers for obs and par */
  double *Yptr = REAL(Y);
  double *PARptr = REAL(PAR);

/* integers from input */
  int n = asInteger(N);
  int k = asInteger(K);

/* loop counters */
  int i, j, l;

/* residuals for single obs */
  SEXP ymu;
  PROTECT(ymu = allocVector(REALSXP, k));
  double *ymuptr = REAL(ymu);

/* inverse of L matrix times residuals for single obs */
  SEXP Linvymu;
  PROTECT(Linvymu = allocVector(REALSXP, k));
  double *Linvymuptr = REAL(Linvymu);

/* return value: vector of log-densities */
  SEXP d;
  PROTECT(d = allocVector(REALSXP, n));
  double *dptr = REAL(d);

/* initialize three terms of log-density */
  double lpi = -0.5 * (float)k * log(2.0 * 3.14159265358979323846);
  double det = 0.0;
  double nrm = 0.0;

/* HERE comes the computation of the log-density */
  for(i = 0; i < n; i++) {
/* compute term2 the determinante */
    det = 0.0;
    for(j = 0; j < k; j++) {
      det += PARptr[i + n * (j + k)];
    }

/* compute term3 the norm */
    /* y - mu */
    for(j = 0; j < k; j++) {
      ymuptr[j] = Yptr[i + j * n] - PARptr[i + j * n];
    }

    /* inverse of L times vector of residuals */
    for(j = 0; j < k; j++) {
      Linvymuptr[j] = PARptr[i + n * (j + k)] * ymuptr[j];
      for(l = 0; l < j; l++) {
	/* The l in the index subsetting of PARptr is wrong, but
	   must be replaced to grep the right element of the lower
	   triangular. */
        Linvymuptr[j] += PARptr[i + n * (l + k + k)] * ymuptr[l];
      }
    }

    /* L2-norm */
    nrm = 0.0;
    for(j = 0; j < k; j++) {
      nrm += pow(Linvymuptr[j], 2.0);
    }

/* sum up terms of log likelihood */
    dptr[i] = lpi + det - 0.5 * nrm;
  }

/* return */
  UNPROTECT(3);
  return d;
}



