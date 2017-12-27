#include "cuhre_decl.h"
#include "cuhre_util.h"

extern void cuhreDoSample2(count n, ctreal *x,  real *f,
			  double* fun(double*), Glob *globdim);

// Compilation note for R interface: move Rule.c into Rule.h
/*
	Rule.c
		integration with cubature rules
		code lifted with minor modifications from DCUHRE
		by J. Berntsen, T. Espelid, and A. Genz
		this file is part of Divonne
*/

#define RESCALE(a, d,  lower_, upper_) (a * (upper_[d] - lower_[d]) + lower_[d])

/*********************************************************************/
void cuhreRuleFree(Rule *rule)
{
  free(rule->first);
}

/*********************************************************************/
real *cuhreExpandFS(cBounds *b, real *g, real *x, count ndim_ )
{
  count dim, ndim = ndim_;

next:
  /* Compute centrally symmetric sum for permutation of G */

  for( dim = 0; dim < ndim; ++dim )
    *x++ = (.5 + g[dim])*b[dim].lower + (.5 - g[dim])*b[dim].upper;

  for( dim = 0; dim < ndim; ) {
    g[dim] = -g[dim];
    if( g[dim++] < 0 ) goto next;
  }

  /* Find next distinct permutation of G and loop back for next sum.
     Permutations are generated in reverse lexicographic order. */

  for( dim = 1; dim < ndim; ++dim ) {
    ctreal gd = g[dim];
    if( g[dim - 1] > gd ) {
      count i, ix = 0, j = dim, dx = dim - 1;
      for( i = 0; i < --j; ++i ) {
        ctreal tmp = g[i];
        g[i] = g[j];
        g[j] = tmp;
        if( tmp <= gd ) --dx;
        if( g[i] > gd ) ix = i;
      }
      if( g[dx] <= gd ) dx = ix;
      g[dim] = g[dx];
      g[dx] = gd;
      goto next;
    }
  }

  /* Restore original order to generators */
  for( dim = 0; dim < --ndim; ++dim ) {
    ctreal tmp = g[dim];
    g[dim] = g[ndim];
    g[ndim] = tmp;
  }

  return x;
}

/*********************************************************************/


 void cuhreSample2(cRule *rule, void *voidregion, cint flags,
		  double* fun(double*), Glob *globdim)
{
  CUHRETYPEDEFREGION;
  CUHRETYPEDEFSET;
  count ndim_ = globdim->ndim_;
  count ncomp_  = globdim->ncomp_;

  Region *const region = (Region *)voidregion;
  ctreal vol = ldexp(1., -region->div);

  real *x = rule->x, *f = rule->f;
  Set *first = (Set *)rule->first, *last = (Set *)rule->last, *s;
  ctreal *errcoeff = rule->errcoeff;
  ctreal ratio = Sq(first[2].gen[0]/first[1].gen[0]);

  ccount offset = 2*ndim_*ncomp_;
  count dim, comp, rul, n, maxdim = 0;
  real maxrange = 0;

  for( dim = 0; dim < ndim_; ++dim ) {
    cBounds *b = &region->bounds[dim];
    ctreal range = b->upper - b->lower;
    if( range > maxrange ) {
      maxrange = range;
      maxdim = dim;
    }
  }

  for( s = first; s <= last; ++s )
    if( s->n ) x = cuhreExpandFS(region->bounds, s->gen, x, ndim_);

  cuhreDoSample2(rule->n, rule->x, f, fun, globdim);

  for( comp = 0; comp < ncomp_; ++comp ) {
    Result *r = &region->result[comp];
    real sum[nrules];
    ctreal *f1 = f;
    ctreal base = *f1*2*(1 - ratio);
    real maxdiff = 0;
    count bisectdim = maxdim;

    for( dim = 0; dim < ndim_; ++dim ) {
      ctreal *fp = f1 + ncomp_;
      ctreal *fm = fp + ncomp_;
      ctreal fourthdiff = fabs(base +
        ratio*(fp[0] + fm[0]) - (fp[offset] + fm[offset]));
      f1 = fm;
      if( fourthdiff > maxdiff ) {
        maxdiff = fourthdiff;
        bisectdim = dim;
      }
    }
    r->bisectdim = bisectdim;

    f1 = f++;
    Zap(sum);
    for( s = first; s <= last; ++s )
      for( n = s->n; n; --n ) {
        ctreal fun = *f1;
        f1 += ncomp_;
        for( rul = 0; rul < nrules; ++rul )
          sum[rul] += fun*s->weight[rul];
      }

    /* Search for the null rule, in the linear space spanned by two
       successive null rules in our sequence, which gives the greatest
       error estimate among all normalized (1-norm) null rules in this
       space. */

    for( rul = 1; rul < nrules - 1; ++rul ) {
      real maxerr = 0;
      for( s = first; s <= last; ++s )
        maxerr = Max(maxerr,
          fabs(sum[rul + 1] + s->scale[rul]*sum[rul])*s->norm[rul]);
      sum[rul] = maxerr;
    }

    r->avg = vol*sum[0];
    r->err = vol*(
      (errcoeff[0]*sum[1] <= sum[2] && errcoeff[0]*sum[2] <= sum[3]) ?
        errcoeff[1]*sum[1] :
        errcoeff[2]*Max(Max(sum[1], sum[2]), sum[3]) );
  }

  if( VERBOSE > 2 ) {
    char si[64*NDIM + 128*NCOMP], *p = si;

    for( dim = 0; dim < ndim_; ++dim ) {
      cBounds *b = &region->bounds[dim];
      // Rescale the region bounds
      p += sprintf(p,
        (dim == 0) ? "Region (" REALF ") - (" REALF ")" :
                     "\n       (" REALF ") - (" REALF ")",
		   RESCALE(b->lower,dim, globdim->lower_, globdim->upper_),
			 RESCALE(b->upper,dim, globdim->lower_, globdim->upper_));
    }

    for( comp = 0; comp < ncomp_; ++comp ) {
      cResult *r = &region->result[comp];
      p += sprintf(p, "\n[" COUNT "] "
		   REEL " +- " REEL "\n",
		   comp + 1, r->avg, r->err);
    }

    printf(si); // Print
  }
}
