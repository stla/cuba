/*
	Vegas.c
		Vegas Monte-Carlo integration
		by Thomas Hahn
*/

#include "struct_Random.h"
#include "vegas_util.h"

extern bool vegasBadDimension(cint ndim, cint flags);
extern bool vegasBadComponent(cint ncomp);

extern int vegasIntegrate(ctreal *lower, ctreal *upper,
	ctreal epsrel, ctreal epsabs,
  cint flags, cnumber mineval, cnumber maxeval,
  cnumber nstart, cnumber nincrease,
	real *integral, real *erreur, real *prob, double* fun(double*));

/* Compilation note for R interface: remove
 int EXPORT(vegasnbatch) = 1000;
 int EXPORT(vegasgridno) = 0;*/
 int EXPORT(vegasnbatch);
 int EXPORT(vegasgridno);
 char EXPORT(vegasstate)[MAXSTATESIZE] = "";

/* Compilation note for R interface: modif #define Print(s) puts(s); fflush(stdout) */
//#define Print(s) Rprintf(s)


/*********************************************************************/
void EXPORT(Vegas)(ccount ndim, ccount ncomp,
  double* fun(double*),
  ctreal *lower, ctreal *upper,
  ctreal epsrel, ctreal epsabs,
  cint flags, cnumber mineval, cnumber maxeval,
  cnumber nstart, cnumber nincrease,
  number *pneval, int *pfail,
  real *integral, real *erreur, real *prob)
{
  ndim_ = ndim;
  ncomp_ = ncomp;

  if( vegasBadComponent(ncomp) || vegasBadDimension(ndim, flags) ) *pfail = -1;
  else {
    neval_ = 0;
    *pfail = vegasIntegrate(lower, upper,
		  				epsrel, epsabs,
      				flags, mineval, maxeval, nstart, nincrease,
      				integral, erreur, prob, fun);
    *pneval = neval_;
  }
}

/*********************************************************************/
 void (vegas)(ccount ndim, ccount ncomp,
  double* fun(double*),
  ctreal *lower, ctreal *upper,
  ctreal epsrel, ctreal epsabs,
  cint flags, cnumber mineval, cnumber maxeval,
  cnumber nstart, cnumber nincrease,
  int *pneval, int *pfail,
  real *integral, real *erreur, real *prob)
{
  /* make sure the filename is null-terminated */
  if( *EXPORT(vegasstate) ) {
    char *p;
    EXPORT(vegasstate)[sizeof(EXPORT(vegasstate)) - 1] = 0;
    if( (p = strchr(EXPORT(vegasstate), ' ')) ) *p = 0;
  }

  EXPORT(Vegas)(ndim, ncomp,
		fun,
		lower, upper,
		epsrel, epsabs,
    flags, mineval, maxeval,
    nstart, nincrease,
    pneval, pfail,
    integral, erreur, prob);
}

/*********************************************************************/
/* Hvegas :  interface between Haskell and vegas                     */
/*********************************************************************/
void Hvegas(int ndim, int ncomp,
	    double* fun(double*),
	    real *lower, real *upper,
	    double epsrel, double epsabs,
	    int seed, int nbatch,
	    int gridno,
	    int flags, int mineval, int maxeval,
	    int nstart, int nincrease,
	    char *state,
	    int *pneval, int *pfail,
	    double *integral, double *erreur, double *prob)
{
  // 5/11/2012  if (strlen(*state) >0) {
    strncpy(vegasstate, state, 128);
// 5/11/2012   }

//  if (NA_INTEGER != *pmersenneseed)
    SUFFIX(mersenneseed) = seed;
  // if (NA_INTEGER != *pvegasnbatch)
     EXPORT(vegasnbatch)= nbatch;
  // else
  //  EXPORT(vegasnbatch)=1000;

  // if (NA_INTEGER != *pvegasgridno)
     EXPORT(vegasgridno)= gridno;
  // else
  //  EXPORT(vegasgridno)=0;

  /* call vegas */
  vegas(ndim, ncomp, fun,
		lower, upper,
  	epsrel, epsabs,
		flags, mineval, maxeval,
  	nstart, nincrease,
		pneval, pfail,
		integral, erreur, prob);

} // End Hvegas
