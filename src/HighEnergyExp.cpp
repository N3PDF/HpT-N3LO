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

#include "../include/HighEnergyExp.h"

#include <higgs-fo/params.h>

#include <cmath>
#include <ostream>

HighEnergyExp::HighEnergyExp(int order, int channel, void *params) {
  PhysParams param = *reinterpret_cast<PhysParams *>(params);

  ORD = order;
  CHANNEL = channel;
  NC = static_cast<long double>(param.nc);
  NF = static_cast<long double>(param.nf);

  CA = NC;
  CF = (NC * NC - 1.) / (2. * NC);

  MH2 = static_cast<long double>(std::pow(param.mh, 2));
  MUR2 = static_cast<long double>(std::pow(param.mur, 2));
  MUF2 = static_cast<long double>(std::pow(param.muf, 2));

  aass = static_cast<long double>(param.alphas);
  SROOT = static_cast<long double>(param.sroot);
  SIGMA0 = static_cast<long double>(param.sigma0);
}

HighEnergyExp::~HighEnergyExp() {}

//==========================================================================================//
//                               gg-channel NLO pt-distributions //
//------------------------------------------------------------------------------------------//

long double HighEnergyExp::C1ggx(long double x, long double xp) {
  long double result = 2 * CA / M_PIl / xp;
  return result;
}

long double HighEnergyExp::C2ggx(long double x, long double xp) {
  long double result =
      4 * std::pow(CA, 2) / std::pow(M_PIl, 2) * (log(xp) / xp) * std::log(x);
  return result;
}

std::complex<long double> HighEnergyExp::C1gg(std::complex<long double> NN,
                                              long double xp) {
  std::complex<long double> result;
  result = 2 * CA / M_PIl / xp / NN;
  return result;
}

std::complex<long double> HighEnergyExp::C2gg(std::complex<long double> NN,
                                              long double xp) {
  std::complex<long double> result;
  result = 4 * std::pow(CA, 2) / std::pow(M_PIl, 2) * (log(xp) / xp) /
           std::pow(NN, 2);
  return result;
}

//==========================================================================================//
//                             Expanded Resummed Expression //
//------------------------------------------------------------------------------------------//

// TODO: change this into a vector
std::complex<long double> HighEnergyExp::HighEnergyExpExpr(
    std::complex<long double> N, long double pt) {
  // TODO: re-check definition MH2 vs. Qs2
  long double xp = std::pow(pt, 2) / MH2;
  std::complex<long double> zero(0., 0.);
  std::complex<long double> result(0., 0.);

  switch (ORD) {
      // TODO: Match CASES with HpT-MON
    case (0):  // order as^1
    {
      if ((CHANNEL == 0) || (CHANNEL == 5)) {
        result += aass * C1gg(N, xp);
      }
      if ((CHANNEL == 1) || (CHANNEL == 5)) {
        result += zero;
      }
      if ((CHANNEL == 2) || (CHANNEL == 5)) {
        result += zero;
      }
    } break;
    case (1):  // order as^2
    {
      if ((CHANNEL == 0) || (CHANNEL == 5)) {
        result += std::pow(aass, 2) * C2gg(N, xp);
      }
      if ((CHANNEL == 1) || (CHANNEL == 5)) {
        return zero;
      }
      if ((CHANNEL == 2) || (CHANNEL == 5)) {
        return zero;
      }
    } break;
  }

  return 2. * pt / MH2 * SIGMA0 * result;
}

// TODO: change this into a vector
long double HighEnergyExp::HighEnergyExpExprX(long double x, long double pt) {
  // TODO: re-check definition MH2 vs. Qs2
  long double xp = std::pow(pt, 2) / MH2;
  long double zero = 0.;
  long double result = 0.;

  switch (ORD) {
      // TODO: Match CASES with HpT-MON
    case (0):  // order as^1
    {
      if ((CHANNEL == 0) || (CHANNEL == 5)) {
        result += aass * C1ggx(x, xp);
      }
      if ((CHANNEL == 1) || (CHANNEL == 5)) {
        result += zero;
      }
      if ((CHANNEL == 2) || (CHANNEL == 5)) {
        result += zero;
      }
    } break;
    case (1):  // order as^2
    {
      if ((CHANNEL == 0) || (CHANNEL == 5)) {
        result += std::pow(aass, 2) * C2ggx(x, xp);
      }
      if ((CHANNEL == 1) || (CHANNEL == 5)) {
        return zero;
      }
      if ((CHANNEL == 2) || (CHANNEL == 5)) {
        return zero;
      }
    } break;
  }

  return 2. * pt / MH2 * SIGMA0 * result;
}
