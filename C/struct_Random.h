#ifndef __struct_random_h__
#define  __struct_random_h__
#include "common_stddecl.h"


struct {
  real norm;
  number v[SOBOL_MAXDIM][30], prev[SOBOL_MAXDIM];
  number seq;
} sobol_;

typedef unsigned int state_t;
 unsigned int SUFFIX(mersenneseed);

struct {
  state_t state[MERSENNE_N];
  count next;
} mersenne_;
#endif
