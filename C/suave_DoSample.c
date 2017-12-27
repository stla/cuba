#include "common_stddecl.h"
#include "suave_util.h"
/*********************************************************************/
void Sintegrand(ccount ndim, ctreal xx[], ccount ncomp,
       ctreal *lower, ctreal *upper, real ff[],
		       ctreal *weight, double* fun(double*))
{
  // SEXP args, argw, s, t, resultsxp;
  int i;
  double* args = (double*) malloc(ndim * sizeof(double));
  double rdbounds = 1;
  for (i =0; i<ndim; i++){
    args[i] = xx[i] * (upper[i] - lower[i]) + lower[i];
    rdbounds *= upper[i] - lower[i];
  }
  double* result = fun(args);
  free(args);
//  REAL(argw)[ 0]=*weight;

 for (i =0; i<ncomp;  i++) {
   ff[i] = result[i] * rdbounds;
 }
 free(result);
} // End Sintegrand

 void suaveDoSample(number n, ctreal *w, ctreal *x, real *f, double* fun(double*))
{
  neval_ += n;
  while( n-- ) {
    Sintegrand(ndim_, x, ncomp_, lower_, upper_, f, w++, fun);
//    integrand_(ndim_, x, ncomp_, lower_, upper_, rdbounds_, f, w++);
    x += ndim_;
    f += ncomp_;
  }
}
