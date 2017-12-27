#include "common_stddecl.h"
#include "vegas_util.h"

/*********************************************************************/
void VIntegrand(ccount ndim, ctreal xx[], ccount ncomp,
       ctreal *lower, ctreal *upper, real ff[],
		       ctreal *weight, double* fun(double*))
{
  int i;
  double* args = (double*) malloc(ndim * sizeof(double));
  double rdbounds = 1;
  for (i =0; i<ndim; i++){
    args[i] = xx[i] * (upper[i] - lower[i]) + lower[i];
    rdbounds *= upper[i] - lower[i];
  }
  double* result = fun(args);
  free(args);
  for (i =0; i<ncomp;  i++) {
    ff[i] = result[i] * rdbounds;
  }
  free(result);
} // End VIntegrand


 // void vegasDoSample(number n, ctreal *w, ctreal *x,
 //   ctreal *lower, ctreal *upper, ctreal prdbounds, real *f, double* fun(double*))
void vegasDoSample(number n, ctreal *w, ctreal *x, ctreal *lower, ctreal *upper,
  real *f, double* fun(double*))
{
  neval_ += n;
  while( n-- ) {
    VIntegrand(ndim_, x, ncomp_, lower, upper, f, w++, fun);
//    integrand_(&ndim_, x, &ncomp_,  lower, upper, prdbounds, f, w++);
    x += ndim_;
    f += ncomp_;
  }


}
