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
  /* std::vector<std::complex<long double>> HighEnergyExpExpr( */
  /*     std::complex<long double> N, long double xp */
  /* ); */
  std::complex<long double> HighEnergyExpExpr(std::complex<long double> N,
                                              long double pt);

  // Matching Coefficients
  long double Hth1gggH(long double xp);
  long double Hth1gqqH(long double xp);
  long double Hth1qqgH(long double xp);

  std::complex<long double> C1gg(std::complex<long double> NN, long double xp);
  std::complex<long double> C2gg(std::complex<long double> NN, long double xp);
  // LO pt-distribution
  std::complex<long double> LOgggH(std::complex<long double> NN,
                                   long double xp);
  std::complex<long double> LOgqqH(std::complex<long double> NN,
                                   long double xp);
  std::complex<long double> LOqqgH(std::complex<long double> NN,
                                   long double xp);

  // Matching function G)
  long double GOgggH(long double);
  long double GOgqqH(long double);
  long double GOqqgH(long double);

 private:
  int NC, NF, CA, ORD, CHANNEL;
  long double LF, LR, LQ;
  long double MH2, MUF2, MUR2, SROOT;
  long double CF, aass, SIGMA0;

  // Beta functions
  long double Beta0;
  long double Beta1;
  long double Beta2;

  // Cusp Anomalous Dimensions
  long double Ath1g;
  long double Ath2g;
  long double Ath1q;
  long double Ath2q;
  long double Bth1g;
  long double Bth1q;

  // Sigma functions
  // gg->g
  long double Sigma22ggg(long double xp);
  long double Sigma21ggg(long double xp);
  long double Sigma20ggg(long double xp);
  // gq->g
  long double Sigma22gqg(long double xp);
  long double Sigma21gqg(long double xp);
  long double Sigma20gqg(long double xp);
  // qq->g
  long double Sigma22qqg(long double xp);
  long double Sigma21qqg(long double xp);
  long double Sigma20qqg(long double xp);

  // Zeta functions
  const long double zeta2 = gsl_sf_zeta_int(2);
  const long double zeta3 = gsl_sf_zeta_int(3);
  const long double zeta4 = gsl_sf_zeta_int(4);

  // EulerGamma
  // TODO: Move this to global
  const long double EulerGamma = 0.5772156649;
};
