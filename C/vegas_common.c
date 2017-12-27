/*
	common.c
		include most of the modules
		this file is part of Vegas
		last modified 14 Feb 05 th
*/


#include "common_stddecl.h"


 bool vegasBadDimension(cint ndim, cint flags)
{
#if NDIM > 0
  if( ndim > NDIM ) return true;
#endif
  return ndim < SOBOL_MINDIM || (!PSEUDORNG && ndim > SOBOL_MAXDIM);
}


 bool vegasBadComponent(cint ncomp)
{
#if NCOMP > 0
  if( ncomp > NCOMP ) return true;
#endif
  return ncomp < 1;
}

