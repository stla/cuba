#ifndef __suave_util_h__
#define __suave_util_h__
/*
	util.c
		Utility functions
		this file is part of Suave
		last modified 9 Feb 05 th
*/


#include "suave_decl.h"
/* globals */
 real *lower_, *upper_; //  rdbounds_;
 count ndim_, ncomp_;
 int nregions_;
 number neval_;
//  Integrand integrand_;

#define RegionAlloc(p, n, nnew) \
  MemAlloc(p, sizeof(Region) + \
              (n)*(ndim_ + ncomp_ + 1)*sizeof(real) + \
              (nnew)*ndim_*sizeof(bin_t))


#ifdef DEBUG
#include "common_debug.h"
#endif
#endif
