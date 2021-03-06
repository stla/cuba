/*
	Sample.c
		the sampling step of Suave
		this file is part of Suave
		last modified 9 Feb 05 th
*/

// To rescale into the user unit

#include "suave_decl.h"
#include "suave_util.h"

extern void GetRandom(real *x, count ndim);
extern void suaveDoSample(number n, ctreal *w, ctreal *x, real *f, double* fun(double*));

#define RESCALE(a, d) (a * (upper_[d] - lower_[d]) + lower_[d])

typedef struct {
  real sum, sqsum;
  real weight, weightsum, avg, avgsum;
  real guess, chisum, chisqsum;
} Cumulants;

/*********************************************************************/

 void suaveSample(cnumber nnew, void *voidregion,
  real *lastw, real *lastx, real *lastf, cint flags, double* fun(double*))
{
  SUAVETYPEDEFREGION;

  Region *const region = (Region *)voidregion;
  count comp, dim, df;
  number n;
  Cumulants cumul[NCOMP];
  char **ss=NULL, *s="";
  ccount chars = 128*(region->div + 1);

  ctreal jacobian = 1/ldexp((real)nnew, region->div);
  real *w = lastw, *f = lastx;
  bin_t *bin = (bin_t *)(lastf + nnew*ncomp_);


  for( n = nnew; n; --n ) {
    real weight = jacobian;

    GetRandom(f, ndim_);

    for( dim = 0; dim < ndim_; ++dim ) {
      cBounds *b = &region->bounds[dim];
      ctreal pos = *f*NBINS;
      ccount ipos = (count)pos;
      ctreal prev = (ipos == 0) ? 0 : b->grid[ipos - 1];
      ctreal diff = b->grid[ipos] - prev;
      *f++ = b->lower + (prev + (pos - ipos)*diff)*(b->upper - b->lower);
      *bin++ = ipos;
      weight *= diff*NBINS;
    }

    *w++ = weight;
  }

  suaveDoSample(nnew, lastw, lastx, lastf, fun);

  *(w - 1) = -*(w - 1);
  lastw = w;
  w = region->w;
  region->n = lastw - w;

  if( VERBOSE > 2 ) {
    char *p0;
     MemAlloc(ss,  ndim_*64 + ncomp_*(sizeof(char *) + chars));
    s = (char *)(ss + ncomp_);
    p0 = s + ndim_*64;
    for( comp = 0; comp < ncomp_; ++comp ) {
      ss[comp] = p0;
      p0 += chars;
    }
  }

  Zap(cumul);
  df = n = 0;

  while( w < lastw ) {
    cbool final = (*w < 0);
    ctreal weight = fabs(*w++);
    ++n;

    for( comp = 0; comp < ncomp_; ++comp ) {
      Cumulants *c = &cumul[comp];

      ctreal wfun = weight*(*f++);
      c->sum += wfun;
      c->sqsum += Sq(wfun);

      if( final ) {
        if( n > 1 ) {
          real wi = Weight(c->sum, c->sqsum, n);
          c->weightsum += c->weight = wi;
          c->avgsum += c->avg = wi*c->sum;

          if( VERBOSE > 2 ) {
            ctreal sig = sqrt(1/wi);
            ss[comp] += (df == 0) ?
              sprintf(ss[comp], "\n[" COUNT "] "
                REEL " +- " REEL " (" NUMBER ")", comp + 1,
                c->sum, sig, n) :
              sprintf(ss[comp], "\n    "
                REEL " +- " REEL " (" NUMBER ")",
                c->sum, sig, n);
          }

          if( df == 0 ) c->guess = c->sum;
          else {
            c->chisum += wi *= c->sum - c->guess;
            c->chisqsum += wi*c->sum;
          }
        }

        c->sum = c->sqsum = 0;
      }
    }

    if( final ) ++df, n = 0;
  }

  region->df = --df;

  for( comp = 0; comp < ncomp_; ++comp ) {
    Result *r = &region->result[comp];
    Cumulants *c = &cumul[comp];
    ctreal sigsq = 1/c->weightsum;
    ctreal avg = sigsq*c->avgsum;

    if( LAST ) {
      r->sigsq = 1/c->weight;
      r->avg = r->sigsq*c->avg;
    }
    else {
      r->sigsq = sigsq;
      r->avg = avg;
    }
    r->err = sqrt(r->sigsq);

    r->chisq = (sigsq < .9*NOTZERO) ? 0 : c->chisqsum - avg*c->chisum;
      /* This catches the special case where the integrand is constant
         over the entire region.  Unless that constant is zero, only the
         first set of samples will have zero variance, and hence weight
         (n - 1) 1e30 (see above).  All other sets have been sampled
         from a non-constant weight function and therefore inevitably
         show some variance.  This is an artificial effect, brought about
         by the fact that the constancy of the integrand in the region is
         seen only in this subdivision, and can degrade the chi-square
         score quite a bit.  If the constancy was determined from more
         than two samples (hence .9*NOTZERO), the chi-squares from the
         other sets are removed here. */
  }

  if( VERBOSE > 2 ) {
    char *p = s;
    char *p0 = p + ndim_*64;

    for( dim = 0; dim < ndim_; ++dim ) {
      cBounds *b = &region->bounds[dim];
      p += sprintf(p,
        (dim == 0) ? "\nRegion (" REALF ") - (" REALF ")" :
                     "\n       (" REALF ") - (" REALF ")",
        RESCALE(b->lower,dim),  RESCALE(b->upper,dim));
    }

    for( comp = 0; comp < ncomp_; ++comp ) {
      cResult *r = &region->result[comp];
      p += sprintf(p, "%s  \tchisq " REEL " (" COUNT " df)\n",
        p0, r->chisq, df);
      p0 += chars;
    }

    printf(s); // Print
    free(ss);
  }
}
