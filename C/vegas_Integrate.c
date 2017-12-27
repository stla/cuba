
#include "common_ChiSquare.h"
#include "vegas_decl.h"
#include "vegas_util.h"

extern void IniRandom(cnumber n, cint flags, count ndim);
extern void GetRandom(real *x, count ndim);
extern  void SkipRandom(cnumber n, count ndim);
extern void PutGrid(Grid *grid);
  extern void GetGrid(Grid *grid);
  extern void vegasRefineGrid(Grid grid, Grid margsum, cint flags);
// extern void vegasDoSample(number n, ctreal *w, ctreal *x,
// 			  ctreal *lower, ctreal *upper, ctreal prdbounds, real *f, double* fun(double*));
extern void vegasDoSample(number n, ctreal *w, ctreal *x, ctreal *lower, ctreal *upper,
			  real *f, double* fun(double*));

extern void decodflags(cint flags,
		int *smooth,
		int *pseudorandom,
		int *final,
		int *verbose);

extern char EXPORT(vegasstate)[MAXSTATESIZE];
extern  int EXPORT(vegasnbatch);
extern   int EXPORT(vegasgridno);


// #include <R_ext/Utils.h> // to allow interruptions

/*
	Integrate.c
		integrate over the unit hypercube
		this file is part of Vegas
		last modified 17 Dec 07 th
*/


/************************************************************* */
 int vegasIntegrate(  ctreal *lower, ctreal *upper,
		       ctreal epsrel, ctreal epsabs,
  cint flags, cnumber mineval, cnumber maxeval,
  cnumber nstart, cnumber nincrease,
  real *integral, real *erreur, real *prob, double* fun(double*))
{
  real *sample;
  count dim, comp;
  int fail = 1;
  struct {
    count niter;
    number nsamples, neval;
    Cumulants cumul[MAXNCOMP];
    Grid grid[MAXNDIM];
  } state;
  int statemsg = VERBOSE;
  struct stat st;

  if( VERBOSE > 1 ) {
    char s[512];
    int smooth, pseudorandom, final, verbose;
    decodflags( flags,
		&smooth,
		&pseudorandom,
		&final,
		&verbose);
    sprintf(s, "Vegas input parameters:\n"
      "  ndim " COUNT "\n  ncomp " COUNT "\n"
      "  rel.tol " REEL "\n  abs.tol " REEL "\n"
      "  smooth %d\n  pseudo.random  %d\n  final %d\n  verbose %d\n  min.eval " NUMBER "\n  max.eval " NUMBER "\n"
      "  nstart " NUMBER "\n  nincrease " NUMBER "\n"
      "  vegas.gridno %d\n  vegas.state \"%s\"\n",
      ndim_, ncomp_,
      epsrel, epsabs,
      smooth, pseudorandom, final, verbose, mineval, maxeval,
      nstart, nincrease,
      EXPORT(vegasgridno), EXPORT(vegasstate));
    printf(s); // Print
  }

#ifdef MLVERSION
  if( setjmp(abort_) ) goto abort;
#endif

  IniRandom(2*maxeval, flags, ndim_);

  if( *EXPORT(vegasstate) && stat(EXPORT(vegasstate), &st) == 0 &&
      st.st_size == sizeof(state) && (st.st_mode & 0400) ) {
    cint h = open(EXPORT(vegasstate), O_RDONLY);
    read(h, &state, sizeof(state));
    close(h);
    SkipRandom(neval_ = state.neval, ndim_);

    if( VERBOSE ) {
      char s[256];
      sprintf(s, "\nRestoring state from %s.", EXPORT(vegasstate));
      printf(s); // Print
    }
  }
  else {
    state.niter = 0;
    state.nsamples = nstart;
    Zap(state.cumul);
    GetGrid(state.grid);
  }

  SamplesAlloc(sample, EXPORT(vegasnbatch));

  /* main iteration loop */

  for( ; ; ) {
//R_CheckUserInterrupt(); // permettre a l'utilisateur d'interrompre
    number nsamples = state.nsamples;
    ctreal jacobian = 1./nsamples;
    Grid margsum[NCOMP][NDIM];

    Zap(margsum);

    for( ; nsamples > 0; nsamples -= EXPORT(vegasnbatch) ) {
      cnumber nbatch = IMin(EXPORT(vegasnbatch), nsamples);
      real *w = sample;
      real *x = w + nbatch;
      real *f = x + nbatch*ndim_;
      real *lastf = f + nbatch*ncomp_;
      bin_t *bin = (bin_t *)lastf;

      while( x < f ) {
        real weight = jacobian;

        GetRandom(x, ndim_);

        for( dim = 0; dim < ndim_; ++dim ) {
          ctreal pos = *x*NBINS;
          ccount ipos = (count)pos;
          ctreal prev = (ipos == 0) ? 0 : state.grid[dim][ipos - 1];
          ctreal diff = state.grid[dim][ipos] - prev;
          *x++ = prev + (pos - ipos)*diff;
          *bin++ = ipos;
          weight *= diff*NBINS;
        }

        *w++ = weight;
      }

//      vegasDoSample(nbatch, sample, w, lower, upper, prdbounds, f);
      vegasDoSample(nbatch, sample, w, lower, upper, f, fun); // bizarre, j'échangerais sample et w

      w = sample;
      bin = (bin_t *)lastf;

      while( f < lastf ) {
        ctreal weight = *w++;

        for( comp = 0; comp < ncomp_; ++comp ) {
          real wfun = weight*(*f++);
          if( wfun ) {
            Cumulants *c = &state.cumul[comp];
            Grid *m = margsum[comp];

            c->sum += wfun;
            c->sqsum += wfun *= wfun;
            for( dim = 0; dim < ndim_; ++dim )
              m[dim][bin[dim]] += wfun;
          }
        }

        bin += ndim_;
      }
    }

    fail = 0;

    /* compute the integral and error values */

    for( comp = 0; comp < ncomp_; ++comp ) {
      Cumulants *c = &state.cumul[comp];
      real avg, sigsq;
      real w = Weight(c->sum, c->sqsum, state.nsamples);

      sigsq = 1/(c->weightsum += w);
      avg = sigsq*(c->avgsum += w*c->sum);

      c->avg = LAST ? (sigsq = 1/w, c->sum) : avg;
      c->err = sqrt(sigsq);
      fail |= (c->err > MaxErr(c->avg));

      if( state.niter == 0 ) c->guess = c->sum;
      else {
        c->chisum += w *= c->sum - c->guess;
        c->chisqsum += w*c->sum;
      }
      c->chisq = c->chisqsum - avg*c->chisum;

      c->sum = c->sqsum = 0;
    }

    if( VERBOSE ) {
      char s[128 + 128*NCOMP], *p = s;

      p += sprintf(p,
        "Iteration " COUNT ":  " NUMBER " integrand evaluations so far",
        state.niter + 1, neval_);

      for( comp = 0; comp < ncomp_; ++comp ) {
        cCumulants *c = &state.cumul[comp];
        p += sprintf(p, "\n[" COUNT "] "
          REEL " +- " REEL "  \tchisq " REEL " (" COUNT " df)\n",
          comp + 1, c->avg, c->err, c->chisq, state.niter);
      }

      printf(s); // Print
    }

    if( fail == 0 && neval_ >= mineval ) {
      if( *EXPORT(vegasstate) ) unlink(EXPORT(vegasstate));
      break;
    }

    if( neval_ >= maxeval && *EXPORT(vegasstate) == 0 ) break;

    if( ncomp_ == 1 )
      for( dim = 0; dim < ndim_; ++dim )
        vegasRefineGrid(state.grid[dim], margsum[0][dim], flags);
    else {
      for( dim = 0; dim < ndim_; ++dim ) {
        Grid wmargsum;
        Zap(wmargsum);
        for( comp = 0; comp < ncomp_; ++comp ) {
          real w = state.cumul[comp].avg;
          if( w != 0 ) {
            ctreal *m = margsum[comp][dim];
            count bin;
            w = 1/Sq(w);
            for( bin = 0; bin < NBINS; ++bin )
              wmargsum[bin] += w*m[bin];
          }
        }
        vegasRefineGrid(state.grid[dim], wmargsum, flags);
      }
    }

    ++state.niter;
    state.nsamples += nincrease;

    if( *EXPORT(vegasstate) ) {
      cint h = creat(EXPORT(vegasstate), 0666);
      if( h != -1 ) {
        state.neval = neval_;
        write(h, &state, sizeof(state));
        close(h);

        if( statemsg ) {
          char s[256];
          sprintf(s, "\nSaving state to %s.\n", EXPORT(vegasstate));
          printf(s); // Print
          statemsg = false;
        }
      }
      if( neval_ >= maxeval ) break;
    }
  }

  for( comp = 0; comp < ncomp_; ++comp ) {
    cCumulants *c = &state.cumul[comp];
    integral[comp] = c->avg;
    erreur[comp] = c->err;
    prob[comp] = ChiSquare(c->chisq, state.niter);
  }

#ifdef MLVERSION
abort:
#endif

  free(sample);
  PutGrid(state.grid);

  return fail;
}
