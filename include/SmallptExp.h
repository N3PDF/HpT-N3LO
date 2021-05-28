/*
 * =====================================================================================
 *
 *       Filename:
 *
 *    Description:
 *
 *        Version:  1.0
 *        Created:
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:
 *   Organization:  N3PDF
 *
 * =====================================================================================
 */

#pragma once

#include <cmath>
#include <gsl/gsl_math.h>
#include <iostream>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <vector>

#include "./AnomalousDim.h"
#include "./ComplexDefs.h"
#include "./MellinFunc.h"

#include "higgs-fo/params.h"

class SmallptExp {
public:
  SmallptExp(int order, int channel, void *params);
  virtual ~SmallptExp();

  // Attribute that compute the expanded results
  std::complex<long double> SmallptExpExpr(std::complex<long double> N,
                                           long double pt);

private:
  AnomDimensions AD;

  int NC, NF, CA, ORD, CHANNEL;
  long double LF, LR, LQ;
  long double QS2, MH2, MUF2, MUR2, SROOT;
  long double CF, aass, SIGMA0;

  // Beta functions
  long double Beta0;
  long double Beta1;
  long double Beta2;

  // Cusp Anomalous Dimensions
  long double Apt1g;
  long double Apt2g;
  long double Bpt1g;
  long double Bpt2g;

  // Fourier Inverse of ln b
  long double LC1(long double xp);
  long double LC2(long double xp);
  long double LC3(long double xp);
  long double LC4(long double xp);

  // Ceofficient functions
  // LO
  std::complex<long double> C1GG(std::complex<long double> N);
  std::complex<long double> C1GQ(std::complex<long double> N);
  std::complex<long double> C2GG(std::complex<long double> N);

  // Zeta functions
  const long double zeta2 = gsl_sf_zeta_int(2);
  const long double zeta3 = gsl_sf_zeta_int(3);
  const long double zeta4 = gsl_sf_zeta_int(4);

  // EulerGamma
  // TODO: Move this to global
  const long double EulerGamma = 0.5772156649;
};
