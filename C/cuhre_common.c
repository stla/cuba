/*
	common.c
		includes most of the modules
		this file is part of Cuhre
*/

#include "common_stddecl.h"
#include "cuhre_util.h"

bool cuhreBadDimension(ccount ndim)
{
#if NDIM > 0
  if( ndim > NDIM ) return true;
#endif
return ndim < 1;
}

bool cuhreBadComponent(cint ncomp)
{
#if NCOMP > 0
  if( ncomp > NCOMP ) return true;
#endif
  return ncomp < 1;
}
