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

#include "../include/ThresXspace.h"

#include <higgs-fo/params.h>

#include <cmath>
#include <ostream>

ThresXspace::ThresXspace(int order, int channel, void *params) {
  PhysParams param = *reinterpret_cast<PhysParams *>(params);

  ORD = order;
  CHANNEL = channel;
  NC = param.nc;
  NF = param.nf;

  CA = NC;
  CF = (NC * NC - 1.) / (2. * NC);

  MH2 = std::pow(param.mh, 2);
  MUR2 = std::pow(param.mur, 2);
  MUF2 = std::pow(param.muf, 2);

  LF = std::log(MH2 / MUF2);
  LR = std::log(MH2 / MUR2);

  aass = param.alphas;
  SROOT = param.sroot;
  SIGMA0 = param.sigma0;

  Beta0 = (11. * CA - 2. * NF) / 6.;
}

ThresXspace::~ThresXspace() {}

//==========================================================================================//
//                                  LO pt-distributions //
//------------------------------------------------------------------------------------------//

// gg->g
double ThresXspace::LOgggH(double x, double xp) {
  // Definitons
  double axp = std::pow(-std::sqrt(xp) + std::sqrt(1. + xp), 2);
  double den = std::sqrt(std::pow(1. - axp * x, 2.) - 4. * axp * x * xp);
  double tau = x * axp;

  double CLOgggH =
      (2. * NC * std::pow(1. - tau + std::pow(tau, 2), 2)) / (den * xp) +
      (-4. * NC * std::pow(1. - tau, 2) * tau +
       2. * NC * std::pow(tau, 2) * xp) /
          den;
  return aass / M_PI * CLOgggH;
}

// gq->q
double ThresXspace::LOgqqH(double x, double pt) { return 0.; }

// qq->g
double ThresXspace::LOqqgH(double x, double pt) { return 0.; }

//==========================================================================================//
//                             Expanded Resummed Expression //
//------------------------------------------------------------------------------------------//

// TODO: change this into a vector
double ThresXspace::ThresXspaceExpr(double x, double N, double pt) {
  // TODO: re-check definition MH2 vs. Qs2
  double xp = std::pow(pt, 2) / MH2;
  double zero = 0;
  double result = 0.;

  switch (ORD) {
      // TODO: Match CASES with HpT-MON
    case (0):  // order as^1
    {
      if ((CHANNEL == 0) || (CHANNEL == 5)) {
        result += LOgggH(x, xp);
      }
      if ((CHANNEL == 1) || (CHANNEL == 5)) {
        result += LOgqqH(x, xp);
      }
      if ((CHANNEL == 2) || (CHANNEL == 5)) {
        result += LOqqgH(x, xp);
      }
    } break;
    case (1):  // order as^2
    {
      if ((CHANNEL == 0) || (CHANNEL == 5)) {
        double axp = std::pow(-std::sqrt(xp) + std::sqrt(1. + xp), 2);

        // Theta functions
        double theta1 = psi0one - psi0half;
        double theta2 = 1. / 2. * std::pow(psi0half, 2) - 1. / 2. * psi1half +
                        1. / 2. * (psi1one + std::pow(psi0one, 2)) -
                        psi0one * psi0half;

        // g functions
        double g2 = 3. * NC;
        double g1 = -Beta0 + 6. * NC * theta1 +
                    2. * NC * std::log((1. - axp) / 2.) -
                    2. * NC * std::log(axp * xp) + 2. * LF +
                    4. * std::log((1. / axp - 1.) / 2.);
        double g0 = -(Beta0 * theta1) + 6. * NC * theta2 +
                    2. * NC * theta1 * std::log((1. - axp) / 2.) -
                    2. * NC * theta1 * std::log(axp * xp) +
                    2. * (2. * NC * theta1 + Beta0) * LF *
                        (LF + std::log((1 / axp) - 1.) / 2.);
        double SIGMAGG =
            g2 * std::pow(std::log(1. - x), 2) + g1 * std::log(1. - x) + g0;
        result += aass / (2. * M_PI) * LOgggH(x, xp) * SIGMAGG;
      }
      if ((CHANNEL == 1) || (CHANNEL == 5)) {
        return zero;
      }
      if ((CHANNEL == 2) || (CHANNEL == 5)) {
        return zero;
      }
    } break;
  }
  result *= std::pow(x, N - 1);  // Mellin factor
  return 2. * pt / MH2 * SIGMA0 * result;
}
