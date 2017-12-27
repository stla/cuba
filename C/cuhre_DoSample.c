#include "common_stddecl.h"
#include "cuhre_util.h"

/*********************************************************************/
/* The function RIntegrand calls the R user function */
/*********************************************************************/
 void Cintegrand(ccount *ndim, ctreal xx[],
		       ccount  *ncomp,
		       //ctreal *lower, ctreal *upper,
          // ctreal prdbounds, // bizarre c'est dans globdim
		       real ff[], double* fun(double*), Glob *globdim)
{
  int i;
  double* xxx = (double*) malloc(*ndim * sizeof(double));
  double rdbounds = 1.0;
  for (i =0; i<*ndim; i++) {
    xxx[i] = xx[i] * (globdim->upper_[i] - globdim->lower_[i]) + globdim->lower_[i];
    rdbounds *= globdim->upper_[i] - globdim->lower_[i];
  }
  double* result = fun(xxx);
  free(xxx);

  for (i =0; i<*ncomp;  i++) {
    ff[i] = result[i] * rdbounds;
  }
  free(result);
} // End theIntegrand


/*********************************************************************/
 void cuhreDoSample2(count n, ctreal *x,  real *f,
		    double* fun(double*), Glob *globdim)
{

  globdim->neval_ += n;

  while( n-- ) {
    Cintegrand(&(globdim->ndim_), x, &(globdim->ncomp_),
                  //globdim->lower_, globdim->upper_,
                  //globdim->prdbounds_,
	                f, fun, globdim);
    x += globdim->ndim_;
    f += globdim->ncomp_;
  }
}
