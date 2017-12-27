#ifndef __vegas_util_h__
#define __vegas_util_h__
/*
	util.c
		Utility functions
		this file is part of Vegas
		last modified 2 Mar 06 th
*/


#include "vegas_decl.h"
/* globals */
// Integrand integrand_;
 count ndim_, ncomp_;
 number neval_;


 Grid *gridptr_[MAXGRIDS];
 count griddim_[MAXGRIDS];


#define SamplesAlloc(p, n) \
  MemAlloc(p, (n)*((ndim_ + ncomp_ + 1)*sizeof(real) + ndim_*sizeof(bin_t)))

#ifdef DEBUG
#include "common_debug.h"
#endif
#endif
