/***********************************************************/
/* Common subroutines for R interface */
/***********************************************************/
#include "common_stddecl.h"
/***********************************************************/
void decodflags(cint flags, 
		int *smooth,
		int *pseudorandom,
		int *final,
		int *verbose)
{
/* decode the argument "flags" */
  int mask =3; /* mask=11 */
  *verbose= (mask & flags);
  mask=4; /* mask=100 */
  *final = (mask & flags) >> 2;
  mask = mask <<1;
  *pseudorandom = (mask & flags) >>3;
  mask = mask <<1;
  *smooth = (mask & flags) >>4;
} // end decodflags

/*********************************************************************/
