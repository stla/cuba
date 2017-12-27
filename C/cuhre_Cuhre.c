/*
	Cuhre.c
		Adaptive integration using cubature rules
		by Thomas Hahn
*/

#include "cuhre_util.h"

extern bool cuhreBadDimension(ccount ndim);
extern bool cuhreBadComponent(cint ncomp);
extern int cuhreIntegrate2(ctreal epsrel, ctreal epsabs,
              cint flags, number mineval, cnumber maxeval, ccount key,
			        real *integral, real *erreur, real *prob,
			        double* fun(double*), Glob *globdim);

/*********************************************************************/
void EXPORT(Cuhre)(ccount ndim, ccount ncomp,
  ctreal epsrel, ctreal epsabs,
  cint flags, cnumber mineval, cnumber maxeval,
  ccount key,
  count *pnregions, number *pneval, int *pfail,
  real *integral, real *erreur, real *prob,
	double* fun(double*), Glob *globdim)
{

  globdim->ndim_ = ndim;
  globdim->ncomp_ = ncomp;

  if( cuhreBadComponent(ncomp) || cuhreBadDimension(ndim) ) *pfail = -1;
  else {
    globdim->neval_ = 0;
    *pfail = cuhreIntegrate2( epsrel, Max(epsabs, NOTZERO),
              flags, mineval, maxeval, key,
			        integral, erreur, prob,
			        fun, globdim);
    *pnregions = globdim->nregions_;
    *pneval =globdim->neval_;
  }
}

/*********************************************************************/
void (cuhre)(ccount ndim, ccount ncomp,
      ctreal epsrel, ctreal epsabs,
      cint flags, cnumber mineval, cnumber maxeval,
      ccount key,
      count *pnregions, number *pneval, int *pfail,
	    real *integral, real *erreur, real *prob,
	    double* fun(double*), Glob *globdim)
{
  EXPORT(Cuhre)(ndim, ncomp,
    epsrel, epsabs,
    flags, mineval, maxeval,
    key,
    pnregions, pneval, pfail,
		integral, erreur, prob, fun, globdim);
}

/*********************************************************************/
/* Hcuhre: interface between Haskell and cuhre                       */
/*********************************************************************/
void Hcuhre(int ndim, int ncomp, double* fun(double*),
      real *lower, real *upper, // real rdbounds,
      double epsrel, double epsabs,
      int flags, int mineval, int maxeval,
	    int key,
      int *pnregions, int *pneval, int *pfail,
      double *integral, double *erreur, double *prob)
{
  Glob globdim;
  globdim.lower_ = lower;
  globdim.upper_ = upper;

 /* call cuhre */
  cuhre(ndim, ncomp,
    epsrel, epsabs,
	  flags, mineval, maxeval,
	  key,
	  pnregions,
	  pneval, pfail,
	  integral, erreur, prob,
	  fun, &globdim);

} // End Hcuhre
