/*
	decl.h
		Type declarations
		this file is part of Suave
		last modified 30 Aug 07 th
*/

//Compilation note for R interface: add ifndef
#ifndef __suave_decl_h__
#define __suave_decl_h__

#include "common_stddecl.h"

#define MINSAMPLES 10

#define NBINS 64

typedef unsigned char bin_t;
/* Note: bin_t must be wide enough to hold the numbers 0..NBINS */

typedef const bin_t cbin_t;

typedef real Grid[NBINS];

typedef const Grid cGrid;

typedef struct {
  real avg, err, sigsq, chisq;
} Result;

typedef const Result cResult;


typedef struct {
  real lower, upper, mid;
  Grid grid;
} Bounds;

typedef const Bounds cBounds;


#define SUAVETYPEDEFREGION \
  typedef struct region { \
    struct region *next; \
    count div, df; \
    number n; \
    Result result[MAXNCOMP]; \
    Bounds bounds[MAXNDIM]; \
    real fluct[MAXNCOMP][MAXNDIM][2]; \
    real w[]; \
  } Region


// typedef void (*Integrand)(ccount, ctreal *, ccount, ctreal *lower,
//   ctreal *upper, ctreal rdbounds, real *, ctreal *);

#endif
