#ifndef __erf_h__
#define __erf_h__
/*
	Erf.c
		Gaussian error function
		= 2/Sqrt[Pi] Integrate[Exp[-t^2], {t, 0, x}]
		Code from Takuya Ooura's gamerf2a.f
		http://www.kurims.kyoto-u.ac.jp/~ooura/gamerf.html
*/

// Compilation note for R interface: add include
#include "common_stddecl.h"

static real Erfc(ctreal x)
{
  static ctreal c[] = {
    2.96316885199227378e-01, 6.12158644495538758e-02,
    1.81581125134637070e-01, 5.50942780056002085e-01,
    6.81866451424939493e-02, 1.53039662058770397e+00,
    1.56907543161966709e-02, 2.99957952311300634e+00,
    2.21290116681517573e-03, 4.95867777128246701e+00,
    1.91395813098742864e-04, 7.41471251099335407e+00,
    9.71013284010551623e-06, 1.04765104356545238e+01,
    1.66642447174307753e-07, 1.48455557345597957e+01,
    6.10399733098688199e+00, 1.26974899965115684e+01 };
  real y = x*x;
  y = exp(-y)*x*(
    c[0]/(y + c[1]) + c[2]/(y + c[3]) +
    c[4]/(y + c[5]) + c[6]/(y + c[7]) +
    c[8]/(y + c[9]) + c[10]/(y + c[11]) +
    c[12]/(y + c[13]) + c[14]/(y + c[15]) );
  if( x < c[16] ) y += 2/(exp(c[17]*x) + 1);
  return y;
}


static real Erf(ctreal x)
{
  static ctreal c[] = {
    1.12837916709551257e+00,
   -3.76126389031833602e-01,
    1.12837916706621301e-01,
   -2.68661698447642378e-02,
    5.22387877685618101e-03,
   -8.49202435186918470e-04 };
  real y = fabs(x);
  if( y > .125 ) {
    y = 1 - Erfc(y);
    return (x > 0) ? y : -y;
  }
  y *= y;
  return x*(c[0] + y*(c[1] + y*(c[2] +
    y*(c[3] + y*(c[4] + y*c[5])))));
}
#endif