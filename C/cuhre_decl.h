/*
	decl.h
		Type declarations
		this file is part of Cuhre
		last modified 8 Apr 09 th
*/
//Compilation note for R interface: add ifndef
#ifndef __cuhre_decl_h__
#define __cuhre_decl_h__

#include "common_stddecl.h"

typedef struct {
  real avg, err;
  count bisectdim;
} Result;

typedef const Result cResult;


typedef struct {
  real avg, err, lastavg, lasterr;
  real weightsum, avgsum;
  real guess, chisum, chisqsum, chisq;
} Totals;

typedef const Totals cTotals;


typedef struct {
  real lower, upper;
} Bounds;

typedef const Bounds cBounds;


typedef struct {
  real *x, *f;
  void *first, *last;
  real errcoeff[3];
  count n;
} Rule;

typedef const Rule cRule;
//typedef void (*Integrand)(ccount *, ctreal *, ccount *, ctreal *lower, ctreal *upper, ctreal prdbounds, real *, SEXP rho, SEXP globf, SEXP globdim);


#define CUHRETYPEDEFREGION \
  typedef struct region { \
    count div; \
    Result result[MAXNCOMP]; \
    Bounds bounds[MAXNDIM]; \
  } Region;



#endif
