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

#include <gsl/gsl_math.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include "./ComplexDefs.h"
#include "./MellinFunc.h"
#include "higgs-fo/params.h"

class HighEnergyExp {
 public:
  HighEnergyExp(int order, int channel, void *params);
  virtual ~HighEnergyExp();

  // Attribute that compute the expanded results
  std::complex<long double> HighEnergyExpExpr(std::complex<long double> N,
                                              long double pt);
  long double HighEnergyExpExprX(long double x,
                                              long double pt);

  // Coefficients in x space
  long double C1ggx(long double x, long double xp);
  long double C2ggx(long double x, long double xp);

  // Coefficients in N space
  std::complex<long double> C1gg(std::complex<long double> NN, long double xp);
  std::complex<long double> C2gg(std::complex<long double> NN, long double xp);

 private:
  int NC, NF, CA, ORD, CHANNEL;
  long double LF, LR, LQ;
  long double MH2, MUF2, MUR2, SROOT;
  long double CF, aass, SIGMA0;
};
