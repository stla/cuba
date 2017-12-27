/*
	Suave.c
		Subregion-adaptive Vegas Monte-Carlo integration
		by Thomas Hahn
*/

#include "common_stddecl.h"
#include "struct_Random.h"
#include "suave_util.h"
extern bool suaveBadDimension(cint ndim, cint flags);
extern bool suaveBadComponent(cint ncomp);
extern  int suaveIntegrate(  ctreal epsrel, ctreal epsabs,
  cint flags, cnumber mineval, cnumber maxeval,
  cnumber nnew, ctreal flatness,
		      real *integral, real *erreur, real *prob, double* fun(double*));

/* Compilation note for R interface: modif #define Print(s) puts(s); fflush(stdout) */
//#define Print(s) Rprintf(s)

/*********************************************************************/
void EXPORT(Suave)(ccount ndim, ccount ncomp,
  double* fun(double*),
  ctreal epsrel, ctreal epsabs,
  cint flags, cnumber mineval, cnumber maxeval,
  cnumber nnew, ctreal flatness,
  int *pnregions, int *pneval, int *pfail,
  real *integral, real *erreur, real *prob)
{
  ndim_ = ndim;
  ncomp_ = ncomp;

  if( suaveBadComponent(ncomp) || suaveBadDimension(ndim, flags) ) *pfail = -1;
  else {
    neval_ = 0;
    nregions_ = 0;
    *pfail = suaveIntegrate(epsrel, Max(epsabs, NOTZERO),
              flags, mineval, maxeval, nnew, flatness,
              integral, erreur, prob, fun);
    *pnregions = nregions_;
    *pneval = neval_;
  }
}

/*********************************************************************/
void (suave)(ccount ndim, ccount ncomp,
  double* fun(double*),
  ctreal epsrel, ctreal epsabs,
  cint flags, cnumber mineval, cnumber maxeval,
  cnumber nnew, ctreal flatness,
  int *pnregions, int *pneval, int *pfail,
  real *integral, real *erreur, real *prob)
{
  EXPORT(Suave)(ndim, ncomp,
    fun,
    epsrel, epsabs,
    flags, mineval, maxeval,
    nnew, flatness,
    pnregions, pneval, pfail,
    integral, erreur, prob);
}

/*********************************************************************/
/* Hsuave: interface between Haskell and suave                       */
/*********************************************************************/
void Hsuave(int ndim, int ncomp, double* fun(double*),
            real *lower, real *upper, //real rdbounds,
            double epsrel, double epsabs,
            int seed, int flags, int mineval, int maxeval,
            int nnew, double flatness,
            int *pnregions, int *pneval, int *pfail,
            double *integral, double *erreur, double *prob)
{
  lower_ = lower;
  upper_ = upper;

  //  if (NA_INTEGER != *pmersenneseed)
      SUFFIX(mersenneseed)= seed;

 /* call suave */
  suave(ndim, ncomp, fun,
  epsrel, epsabs,
	flags, mineval, maxeval,
  nnew, flatness,
	pnregions,
	pneval, pfail,
	integral, erreur, prob);

} // End Rsuave
