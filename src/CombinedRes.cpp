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

#include "../include/CombinedRes.h"

#include <higgs-fo/params.h>

#include <cmath>

CombinedRes::CombinedRes(int order, int channel, std::string pdfname,
                         void *params) {
  PhysParams param = *reinterpret_cast<PhysParams *>(params);

  MH2 = static_cast<long double>(std::pow(param.mh, 2));
  MUR2 = static_cast<long double>(std::pow(param.mur, 2));
  MUF2 = static_cast<long double>(std::pow(param.muf, 2));

  ORD = order;
  EXACT_ORD = order == 0 ? 0 : order - 1;

  SMALLPT = new SmallptExp(order, channel, params);
  THRESHOLD = new ThresExp(order, channel, params);
  MELLIN = new MellinTrans(order, channel, pdfname, params);
  MELLINPARTONIC = new CrossHiggs(EXACT_ORD, channel, pdfname, params);
}

CombinedRes::~CombinedRes() {
  delete MELLIN;
  delete MELLINPARTONIC;
  delete SMALLPT;
  delete THRESHOLD;
}

std::complex<long double> CombinedRes::Matching(std::complex<long double> N,
                                                long double pt, int scheme) {
  // TODO: re-check definition MH2 vs. Qs2
  long double xp = std::pow(pt, 2) / MH2;
  if (scheme == 0)  // small-pt only
  {
    return (0.);
  } else if (scheme == 1)  // threshold only
  {
    return (1.);
  } else  // combined
  {
    long double k = 2.;
    long double m = 9.75;
    return (std::pow(N, m) * std::pow(xp, k)) /
           (1. + std::pow(N, m) * std::pow(xp, 2));
  }
}

std::complex<long double> CombinedRes::CombinedResExpr(
    std::complex<long double> N, long double pt, int scheme) {
  double pp = static_cast<double>(pt);
  // take only real part. Does not work for complex
  double nn = static_cast<double>(N.real());
  std::complex<long double> mres;
  std::vector<double> ResultsMellin;
  std::vector<double> zero(2, 0.0);

  // Compute exact FO from HpT-MON
  if (ORD == 0) {
    ResultsMellin = zero;
  } else {
    ResultsMellin = MELLINPARTONIC->partonichiggsdpt(pp, nn);
  }
  std::vector<long double> ExactMellin(ResultsMellin.begin(),
                                       ResultsMellin.end());

  // Compute approximation from resummations
  std::complex<long double> SptMellin = SMALLPT->SmallptExpExpr(N, pt);
  std::complex<long double> ThresMellin = THRESHOLD->ThresExpExpr(N, pt);
  std::complex<long double> xThresMellin = MELLIN->xSpaceThres(N, pt);

  std::complex<long double> ExactMellinCmpx(ExactMellin[0], 0.);
  mres = (1. - Matching(N, pt, scheme)) * SptMellin +
         Matching(N, pt, scheme) * ThresMellin;
         /* Matching(N, pt, scheme) * xThresMellin; */

  return ExactMellinCmpx + mres;
}
