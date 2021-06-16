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

  LF = std::log(MH2 / MUF2);
  LR = std::log(MH2 / MUR2);

  aass = static_cast<long double>(
      param.alphas);  // TODO; long double-check PI normalization
  SROOT = static_cast<long double>(param.sroot);
  SIGMA0 = static_cast<long double>(param.sigma0);

  // Compute Beta Functions
  Beta0 = (11. * CA - 2. * NF) / (12. * M_PIl);
  Beta1 = ((17. * CA * CA - 5. * CA * NF - 3. * CF * NF) * 2. / 3.) /
          (16. * M_PIl * M_PIl);
  Beta2 = ((2857. / 54. * CA * CA * CA +
            (CF * CF - 205. / 18. * CF * CA - 1415. / 54. * CA * CA) * NF +
            (11. / 9. * CF + 79. / 54. * CA) * NF * NF)) /
          std::pow(4. * M_PIl, 3);

  // One & too loop Cusp Anomalous Dimensions
  // In the following, the 1/PI is contained in the definitions
  // of the cusp Anomalous dimensions not in the normalization of as
  Bth1g = -Beta0;
  Ath1g = CA / M_PIl;
  Ath1q = CF / M_PIl;
  Bth1q = -3. / 4. * CF / M_PIl;
  Ath2g = (CA / 2. * (CA * (67. / 18. - zeta2) - 5. / 9. * NF)) /
          std::pow(M_PIl, 2);
  Ath2q = (CF / 2. * (CA * (67. / 18. - zeta2) - 5. / 9. * NF)) /
          std::pow(M_PIl, 2);
}

HighEnergyExp::~HighEnergyExp() {}

//==========================================================================================//
//                               gg-channel NLO pt-distributions //
//------------------------------------------------------------------------------------------//

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
