/*
	Grid.c
		utility functions for the Vegas grid
		this file is part of Suave
		last modified 15 Feb 08 th
*/

#include "suave_decl.h"
//#include "suave_util.h"

 void suaveRefineGrid(Grid grid, Grid margsum, cint flags)
{
  real avgperbin, thisbin, newcur, delta;
  Grid imp, newgrid;
  int bin, newbin;

  /* smooth the f^2 value stored for each bin */
  real prev = margsum[0];
  real cur = margsum[1];
  real norm = margsum[0] = .5*(prev + cur);
  for( bin = 1; bin < NBINS - 1; ++bin ) {
    ctreal s = prev + cur;
    prev = cur;
    cur = margsum[bin + 1];
    norm += margsum[bin] = (s + cur)/3.;
  }
  norm += margsum[NBINS - 1] = .5*(prev + cur);

  if( norm == 0 ) return;
  norm = 1/norm;

  /* compute the importance function for each bin */
  avgperbin = 0;
  for( bin = 0; bin < NBINS; ++bin ) {
    real impfun = 0;
    if( margsum[bin] > 0 ) {
      ctreal r = margsum[bin]*norm;
      avgperbin += impfun = pow((r - 1)/log(r), 1.5);
    }
    imp[bin] = impfun;
  }
  avgperbin /= NBINS;

  /* redefine the size of each bin */
  cur = newcur = 0;
  thisbin = 0;
  bin = -1;
  for( newbin = 0; newbin < NBINS - 1; ++newbin ) {
    while( thisbin < avgperbin ) {
      thisbin += imp[++bin];
      prev = cur;
      cur = grid[bin];
    }
    thisbin -= avgperbin;
    delta = (cur - prev)*thisbin;
    newgrid[newbin] = SHARPEDGES ?
      cur - delta/imp[bin] :
      (newcur = Max(newcur + 0x1p-48,
        cur - 2*delta/(imp[bin] + imp[IDim(bin - 1)])));
  }
  Copy(grid, newgrid, NBINS - 1);
  grid[NBINS - 1] = 1;
}

/*********************************************************************/

 void Reweight(Bounds *b,
  ctreal *w, ctreal *f, ctreal *lastf, cResult *total, cint flags, cint ndim, cint ncomp)
{
  Grid margsum[ndim];
  real scale[ncomp];
  cbin_t *bin = (cbin_t *)lastf;
  count dim, comp;

  if( ncomp == 1 ) scale[0] = 1;
  else {
    for( comp = 0; comp < ncomp; ++comp )
      scale[comp] = (total[comp].avg == 0) ? 0 : 1/total[comp].avg;
  }

  Zap(margsum);

  while( f < lastf ) {
    real fsq = 0;
    for( comp = 0; comp < ncomp; ++comp )
      fsq += Sq(*f++*scale[comp]);
    fsq *= Sq(*w++);
    if( fsq != 0 )
      for( dim = 0; dim < ndim; ++dim )
        margsum[dim][bin[dim]] += fsq;
    bin += ndim;
  }

  for( dim = 0; dim < ndim; ++dim )
    suaveRefineGrid(b[dim].grid, margsum[dim], flags);
}

/*********************************************************************/

 void StretchGrid(cGrid grid, Grid gridL, Grid gridR)
{
  real prev = 0, cur, step, x;
  count bin = 0;

  while( bin < NBINS ) {
    cur = grid[bin++];
    if( cur >= .5 ) break;
    prev = cur;
  }

  step = (bin - (cur - .5)/(cur - prev))/NBINS;

  prev = x = 0;
  cur = *grid;

  for( bin = 0; bin < NBINS; ++bin ) {
    x += step;
    if( x > 1 ) {
      --x;
      prev = cur;
      cur = *++grid;
    }
    gridL[bin] = 2*(prev + (cur - prev)*x);
  }

  step = 1 - step;
  for( bin = 0; bin < NBINS - 1; ++bin ) {
    x += step;
    if( x > 1 ) {
      --x;
      prev = cur;
      cur = *++grid;
    }
    gridR[bin] = 2*(prev + (cur - prev)*x) - 1;
  }
  gridR[NBINS - 1] = 1;
}
