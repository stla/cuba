
/*
	common.c
		includes most of the modules
		this file is part of Suave
		last modified 14 Feb 05 th
*/


#include "common_stddecl.h"


bool suaveBadDimension(cint ndim, cint flags)
{
#if NDIM > 0
  if( ndim > NDIM ) return true;
#endif
  return ndim < SOBOL_MINDIM || (!PSEUDORNG && ndim > SOBOL_MAXDIM);
}


bool suaveBadComponent(cint ncomp)
{
#if NCOMP > 0
  if( ncomp > NCOMP ) return true;
#endif
  return ncomp < 1;
}
