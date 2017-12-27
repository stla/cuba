
/*
	Integrate.c
		integrate over the unit hypercube
		this file is part of Cuhre
		last modified 8 Apr 09 th
*/
#define POOLSIZE 1024


#include "common_ChiSquare.h"
#include "cuhre_decl.h"
#include "cuhre_util.h"
#include "cuhre_Rule.h"
extern void cuhreSample2(cRule *rule, void *voidregion, cint flags,
			double* fun(double*), Glob *globdim);
extern void cuhreRuleFree(Rule *rule);

extern void decodflags(cint flags,
		int *smooth,
		int *pseudorandom,
		int *final,
		int *verbose);


/*********************************************************************/
 int cuhreIntegrate2(  ctreal epsrel, ctreal epsabs,
  cint flags, number mineval, cnumber maxeval, ccount key,
  real *integral, real *erreur, real *prob,
		      double* fun(double*), Glob *globdim)
{

  CUHRETYPEDEFREGION;
  typedef struct pool {
    struct pool *next;
    Region region[POOLSIZE];
  } Pool;

  count dim, comp, ncur, nregions, ipool, npool;
  int fail = 1;
  Rule rule;
  Totals totals[MAXNCOMP];
  Pool *cur = NULL, *pool;
  Region *region;

  if( VERBOSE > 1 ) {
    char s[512];
    //  smooth ignored here
		int smooth, pseudorandom, final, verbose;
    decodflags( flags,
			&smooth,
			&pseudorandom,
			&final,
			&verbose);
    sprintf(s, "Cuhre input parameters:\n"
      "  ndim " COUNT "\n  ncomp " COUNT "\n"
      "  rel.tol " REEL "\n  abs.tol " REEL "\n"
      "  pseudo.random  %d\n  final %d\n  verbose %d\n  min.eval " NUMBER "\n  max.eval " NUMBER "\n"
      "  key " COUNT "\n",
      globdim->ndim_, globdim->ncomp_,
      epsrel, epsabs,
      pseudorandom, final, verbose, mineval, maxeval,
      key);
    printf(s); // Print
  }

#ifdef MLVERSION
  if( setjmp(abort_) ) goto abort;
#endif

  if( key == 13 && globdim->ndim_ == 2 ) cuhreRule13Alloc(&rule);
  else if( key == 11 && globdim->ndim_ == 3 ) cuhreRule11Alloc(&rule);
  else if( key == 9 ) cuhreRule9Alloc(&rule, globdim->ndim_);
  else if( key == 7 ) cuhreRule7Alloc(&rule, globdim->ndim_);
  else {
    if( globdim->ndim_ == 2 ) cuhreRule13Alloc(&rule);
    else if( globdim->ndim_ == 3 ) cuhreRule11Alloc(&rule);
    else cuhreRule9Alloc(&rule, globdim->ndim_);
  }
  Alloc(rule.x, rule.n*(globdim->ndim_ + globdim->ncomp_));
  rule.f = rule.x + rule.n*globdim->ndim_;

  mineval = IMax(mineval, rule.n + 1);

 Alloc(cur, 1);
  cur->next = NULL;
  ncur = 1;

  region = cur->region;
  region->div = 0;
  for( dim = 0; dim < globdim->ndim_; ++dim ) {
    Bounds *b = &region->bounds[dim];
    b->lower = 0;
    b->upper = 1;
  }

  cuhreSample2(&rule, region, flags, fun, globdim);

  for( comp = 0; comp < globdim->ncomp_; ++comp ) {
    Totals *tot = &totals[comp];
    Result *r = &region->result[comp];
    tot->avg = tot->lastavg = tot->guess = r->avg;
    tot->err = tot->lasterr = r->err;
    tot->weightsum = 1/Max(Sq(r->err), NOTZERO);
    tot->avgsum = tot->weightsum*r->avg;
    tot->chisq = tot->chisqsum = tot->chisum = 0;
  }

  for( nregions = 1; ; ++nregions ) {
    count maxcomp, bisectdim;
    real maxratio, maxerr;
    Result result[MAXNCOMP];
    Region *regionL, *regionR;
    Bounds *bL, *bR;

    if( VERBOSE ) {
      char s[128 + 128*MAXNCOMP], *p = s;

      p += sprintf(p,
        "Iteration " COUNT ":  " NUMBER " integrand evaluations so far",
        nregions,  globdim->neval_);

      for( comp = 0; comp <  globdim->ncomp_; ++comp ) {
        cTotals *tot = &totals[comp];
        p += sprintf(p, "\n[" COUNT "] "
          REEL " +- " REEL "  \tchisq " REEL " (" COUNT " df)",
          comp + 1, tot->avg, tot->err, tot->chisq, nregions - 1);
      }
      p += sprintf(p, "\n");

      printf(s);
    }

    maxratio = -INFTY;
    maxcomp = 0;
    for( comp = 0; comp <  globdim->ncomp_; ++comp ) {
      ctreal ratio = totals[comp].err/MaxErr(totals[comp].avg);
      if( ratio > maxratio ) {
        maxratio = ratio;
        maxcomp = comp;
      }
    }

    if( maxratio <= 1 &&  globdim->neval_ >= mineval ) {
      fail = 0;
      break;
    }

    if(  globdim->neval_ >= maxeval ) break;

    maxerr = -INFTY;
    regionL = cur->region;
    npool = ncur;
    for( pool = cur; pool; npool = POOLSIZE, pool = pool->next )
      for( ipool = 0; ipool < npool; ++ipool ) {
        Region *regioni = &pool->region[ipool];
        ctreal err = regioni->result[maxcomp].err;
        if( err > maxerr ) {
          maxerr = err;
          regionL = regioni;
        }
      }

    if( ncur == POOLSIZE ) {
      Pool *prev = cur;
      Alloc(cur, 1);
      cur->next = prev;
      ncur = 0;
    }
    regionR = &cur->region[ncur++];

    regionR->div = ++regionL->div;
    Copy(result, regionL->result, globdim->ncomp_);
    Copy(regionR->bounds, regionL->bounds, globdim->ndim_);

    bisectdim = result[maxcomp].bisectdim;
    bL = &regionL->bounds[bisectdim];
    bR = &regionR->bounds[bisectdim];
    bL->upper = bR->lower = .5*(bL->upper + bL->lower);

    cuhreSample2(&rule, regionL, flags, fun, globdim);
    cuhreSample2(&rule, regionR, flags, fun, globdim);

    for( comp = 0; comp < globdim->ncomp_; ++comp ) {
      cResult *r = &result[comp];
      Result *rL = &regionL->result[comp];
      Result *rR = &regionR->result[comp];
      Totals *tot = &totals[comp];
      real diff, err, w, avg, sigsq;

      tot->lastavg += diff = rL->avg + rR->avg - r->avg;

      diff = fabs(.25*diff);
      err = rL->err + rR->err;
      if( err > 0 ) {
        ctreal c = 1 + 2*diff/err;
        rL->err *= c;
        rR->err *= c;
      }
      rL->err += diff;
      rR->err += diff;
      tot->lasterr += rL->err + rR->err - r->err;

      tot->weightsum += w = 1/Max(Sq(tot->lasterr), NOTZERO);
      sigsq = 1/tot->weightsum;
      tot->avgsum += w*tot->lastavg;
      avg = sigsq*tot->avgsum;
      tot->chisum += w *= tot->lastavg - tot->guess;
      tot->chisqsum += w*tot->lastavg;
      tot->chisq = tot->chisqsum - avg*tot->chisum;

      if( LAST ) {
        tot->avg = tot->lastavg;
        tot->err = tot->lasterr;
      }
      else {
        tot->avg = avg;
        tot->err = sqrt(sigsq);
      }
    }
  }

  for( comp = 0; comp < globdim->ncomp_; ++comp ) {
    cTotals *tot = &totals[comp];
    integral[comp] = tot->avg;
    erreur[comp] = tot->err;
    prob[comp] = ChiSquare(tot->chisq, nregions - 1);
  }

#ifdef MLVERSION
  if( REGIONS ) {
    MLPutFunction(stdlink, "List", 2);
    MLPutFunction(stdlink, "List", nregions);

    npool = ncur;
    for( pool = cur; pool; npool = POOLSIZE, pool = pool->next )
      for( ipool = 0; ipool < npool; ++ipool ) {
        Region const *region = &pool->region[ipool];
        real lower[NDIM], upper[NDIM];

        for( dim = 0; dim <globdim->ndim_; ++dim ) {
          cBounds *b = &region->bounds[dim];
          lower[dim] = b->lower;
          upper[dim] = b->upper;
        }

        MLPutFunction(stdlink, "Cuba`Cuhre`region", 3);
        MLPutRealList(stdlink, lower, globdim->ndim_);
        MLPutRealList(stdlink, upper, globdim->ndim_);

        MLPutFunction(stdlink, "List", globdim->ncomp_);
        for( comp = 0; comp < globdim->ncomp_; ++comp ) {
          cResult *r = &region->result[comp];
          real res[] = {r->avg, r->err};
          MLPutRealList(stdlink, res, Elements(res));
        }
      }
  }
#endif

#ifdef MLVERSION
abort:
#endif

  while( (pool = cur) ) {
    cur = cur->next;
     free(pool);
  }

   free(rule.x);
  cuhreRuleFree(&rule);

  globdim->nregions_ = nregions;


  return fail;
}
