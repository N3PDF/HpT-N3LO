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

#include "../include/ThresExp.h"

#include <higgs-fo/params.h>

#include <cmath>
#include <ostream>

ThresExp::ThresExp(int order, int channel, void *params) {
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

ThresExp::~ThresExp() {}

//==========================================================================================//
//                                  LO pt-distributions //
//------------------------------------------------------------------------------------------//

// gg->g
std::complex<long double> ThresExp::LOgggH(std::complex<long double> NN,
                                           long double xp) {
  // CHECKED //
  // TODO: verify correspondence with small-pt
  std::complex<long double> CLOgggH;
  std::complex<long double> N = NN;
  std::complex<long double> half(0.5, 0.);
  std::complex<long double> xprad(
      std::pow(std::sqrt(1. + xp) - std::sqrt(xp), 4), 0.);

  CLOgggH =
      2. * aass * CA / std::sqrt(M_PIl) * 1. / xp *
      std::exp(LogGamma(N) - LogGamma(N + 0.5)) *
      (Hyp2F1(half, N, N + 0.5, xprad) -
       2. * (1. + xp) / (std::pow(std::sqrt(1. + xp) + std::sqrt(xp), 2.)) * N /
           (N + 0.5) * Hyp2F1(half, N + 1., N + 1.5, xprad) +
       ((1. + xp) * (3. + xp)) /
           (std::pow(std::sqrt(1. + xp) + std::sqrt(xp), 4.)) * N * (N + 1.) /
           ((N + 0.5) * (N + 1.5)) * Hyp2F1(half, N + 2., N + 2.5, xprad) -
       2. * (1. + xp) / (std::pow(std::sqrt(1. + xp) + std::sqrt(xp), 6.)) * N *
           (N + 1.) * (N + 2.) / ((N + 0.5) * (N + 1.5) * (N + 2.5)) *
           Hyp2F1(half, N + 3., N + 3.5, xprad) +
       1. / (std::pow(std::sqrt(1. + xp) + std::sqrt(xp), 8)) * N * (N + 1.) *
           (N + 2.) * (N + 3.) /
           ((N + 0.5) * (N + 1.5) * (N + 2.5) * (N + 3.5)) *
           Hyp2F1(half, N + 4., N + 4.5, xprad));

  return CLOgggH;
}

// gq->q
std::complex<long double> ThresExp::LOgqqH(std::complex<long double> NN,
                                           long double xp) {
  std::complex<long double> CLOgqqH;
  // TODO: verify correspondence with small-pt
  std::complex<long double> N = NN;
  std::complex<long double> half(0.5, 0.);
  std::complex<long double> xprad(
      std::pow(std::sqrt(1. + xp) - std::sqrt(xp), 4.), 0.);

  CLOgqqH =
      aass * CF / std::sqrt(M_PIl) * 1. / xp *
      std::exp(LogGamma(N) - LogGamma(N + 0.5)) *
      (Hyp2F1(half, N, N + 0.5, xprad) -
       (4. + 3. * xp) / (std::pow(std::sqrt(1. + xp) + std::sqrt(xp), 2.)) * N /
           (N + 0.5) * Hyp2F1(half, N + 1., N + 1.5, xprad) +
       3. * (1. + xp) / (std::pow(std::sqrt(1. + xp) + std::sqrt(xp), 4.)) * N *
           (N + 1.) / ((N + 0.5) * (N + 1.5)) *
           Hyp2F1(half, N + 2., N + 2.5, xprad) -
       1. / (std::pow(std::sqrt(1. + xp) + std::sqrt(xp), 6.)) * N * (N + 1.) *
           (N + 2.) / ((N + 0.5) * (N + 1.5) * (N + 2.5)) *
           Hyp2F1(half, N + 3., N + 3.5, xprad));

  return CLOgqqH;
}

// qq->g
std::complex<long double> ThresExp::LOqqgH(std::complex<long double> NN,
                                           long double xp) {
  std::complex<long double> CLOqqgH;
  // TODO: verify correspondence with small-pt
  std::complex<long double> N = NN;
  std::complex<long double> half(0.5, 0.);
  std::complex<long double> xprad(
      std::pow(std::sqrt(1. + xp) - std::sqrt(xp), 4.), 0.);

  CLOqqgH =
      2. * aass * CF * CF / std::sqrt(M_PIl) * 1. /
      (std::pow(std::sqrt(1. + xp) + std::sqrt(xp), 2.)) *
      std::exp(LogGamma(N) - LogGamma(N + 0.5)) *
      (Hyp2F1(half, N, N + 0.5, xprad) -
       2. * (1. + xp) / (std::pow(std::sqrt(1. + xp) + std::sqrt(xp), 2.)) * N /
           (N + 0.5) * Hyp2F1(half, N + 1., N + 1.5, xprad) +
       1. / (std::pow(std::sqrt(1. + xp) + std::sqrt(xp), 4.)) * N * (N + 1.) /
           ((N + 0.5) * (N + 1.5)) * Hyp2F1(half, N + 2., N + 2.5, xprad));

  return CLOqqgH;
}

//==========================================================================================//
//                  Sigma functions (factor of ln N-enhanced terms) //
//------------------------------------------------------------------------------------------//

// gg->g
long double ThresExp::Sigma22ggg(long double xp) { return 3. * Ath1g / 8.; }

long double ThresExp::Sigma21ggg(long double xp) {
  long double xQp = std::sqrt(1. + xp) + std::sqrt(xp);
  return -Bth1g / 2. - Ath1g / 2. * std::log(std::sqrt(xp) / xQp) -
         2. * Ath1g * LF;
}

long double ThresExp::Sigma20ggg(long double xp) {
  return -3. / 2. * Ath1g * zeta2;
}

//==========================================================================================//
//                  Sigma functions (factor of ln N-enhanced terms) //
//------------------------------------------------------------------------------------------//
// gg->g
long double ThresExp::GOgggH(long double xp) {
  const long double g1gg =
      67. / 36. * CA - 5. / 18. * NF + CA * gsl_sf_zeta(2.) -
      Beta0 * M_PIl * std::log(xp / (1. + xp)) -
      1. / 8. * CA * std::pow(std::log(xp / (1. + xp)), 2) +
      2. * CA * gsl_sf_dilog(1. - std::sqrt(xp / (1. + xp))) +
      CA * std::log(1. - std::sqrt(xp / (1. + xp))) * std::log(xp / (1. + xp)) -
      0.5 * CA * std::log(1. + std::sqrt(xp / (1. + xp))) *
          std::log(xp / (1. + xp)) +
      0.5 * CA * std::pow(std::log(1. + std::sqrt(xp / (1. + xp))), 2) +
      2. * Beta0 * M_PIl *
          std::pow(std::log(1. + std::sqrt(xp / (1. + xp))), 2) +
      CA * gsl_sf_dilog((2. * std::sqrt(xp)) /
                        (std::sqrt(1. + xp) + std::sqrt(xp))) -
      ((CA - NF) *
       (std::sqrt(xp) * std::sqrt(1. + xp) * (1. + xp) - 2. * xp - xp * xp)) /
          (6. * (1. + 8. * xp + 9. * xp * xp));
  return g1gg / M_PIl;
}

// gq->q
long double ThresExp::GOgqqH(long double xp) {
  const long double g1gq =
      -7. / 4. * CF + 134. / 36. * CA - 20. / 36. * NF -
      8. * CF * gsl_sf_zeta(2.) + 12. * CA * gsl_sf_zeta(2.) -
      4. * Beta0 * M_PIl * std::log(xp / (1. + xp)) +
      3. / 2. * CF * std::log(xp / (1. + xp)) -
      0.5 * CA * std::pow(std::log(xp / (1. + xp)), 2) +
      4. * (CF + CA) * gsl_sf_dilog(1. - std::sqrt(xp / (1. + xp))) +
      (2. * (CA - CF) * (1. + 3. * xp + 3. * std::sqrt(xp * (1. + xp)))) /
          (2. * std::sqrt(xp * (1. + xp)) + 1. + 3. * xp) +
      8. * Beta0 * M_PIl * std::log(1. + std::sqrt(xp / (1. + xp))) -
      3. * CF * std::log(1. + std::sqrt(xp / (1. + xp))) +
      2. * CF * std::log(1. - std::sqrt(xp / (1. + xp))) *
          std::log(xp / (1. + xp)) +
      2. * CA * std::log(1. - std::sqrt(xp / (1. + xp))) *
          std::log(xp / (1. + xp)) -
      2. * CF * std::log(1. + std::sqrt(xp / (1. + xp))) *
          std::log(xp / (1. + xp)) -
      2. * CF * std::pow(std::log(1. + std::sqrt(xp / (1. + xp))), 2) +
      4. * CF *
          gsl_sf_dilog(2. * std::sqrt(xp) /
                       (std::sqrt(1. + xp) + std::sqrt(xp)));
  return g1gq / M_PIl;
}

// qq->g
long double ThresExp::GOqqgH(long double xp) {
  const long double g1qq =
      -9. / 2. * CF + 79. / 12. * CA - 5. / 6. * NF +
      12. * CF * gsl_sf_zeta(2.) - 10. * CA * gsl_sf_zeta(2.) -
      (CF - CA) * std::sqrt(1. + xp) / std::sqrt(xp) +
      4. * CF * gsl_sf_dilog(1. - std::sqrt(xp / (1. + xp))) -
      3. / 4. * CF * std::log(xp / (1. + xp)) -
      Beta0 * M_PIl * std::log(xp / (1. + xp)) +
      0.25 * CA * std::pow(std::log(xp / (1. + xp)), 2) -
      0.5 * CF * std::pow(std::log(xp / (1. + xp)), 2) +
      2. * CF * std::log(1. - std::sqrt(xp / (1. + xp))) *
          std::log(xp / (1. + xp)) +
      1.5 * CF * std::log(1. + std::sqrt(xp / (1. + xp))) +
      2. * Beta0 * M_PIl * std::log(1. + std::sqrt(xp / (1. + xp))) +
      CA * std::pow(std::log(1. + std::sqrt(xp / (1. + xp))), 2) -
      CA * std::log(1. + std::sqrt(xp / (1. + xp))) * std::log(xp / (1. + xp)) +
      2. * CA *
          gsl_sf_dilog(2. * std::sqrt(xp) /
                       (std::sqrt(1. + xp) + std::sqrt(xp)));
  return g1qq / M_PIl;
}

//==========================================================================================//
//                             Expanded Resummed Expression //
//------------------------------------------------------------------------------------------//

// TODO: change this into a vector
std::complex<long double> ThresExp::ThresExpExpr(std::complex<long double> N,
                                                 long double pt) {
  // TODO: re-check definition MH2 vs. Qs2
  long double xp = std::pow(pt, 2) / MH2;
  std::complex<long double> zero(0., 0.);
  std::complex<long double> result(0., 0.);

  std::complex<long double> Nbar = N * std::exp(EulerGamma);
  std::complex<long double> LNbar = 2. * aass * Beta0 * std::log(Nbar);

  switch (ORD) {
    // TODO: Match CASES with HpT-MON
  case (0): // order as^1
  {
    if ((CHANNEL == 0) || (CHANNEL == 5)) {
      result += LOgggH(N, xp);
    }
    if ((CHANNEL == 1) || (CHANNEL == 5)) {
      result += LOgqqH(N, xp);
    }
    if ((CHANNEL == 2) || (CHANNEL == 5)) {
      result += LOqqgH(N, xp);
    }
  } break;
  case (1): // order as^2
  {
    if ((CHANNEL == 0) || (CHANNEL == 5)) {
      std::complex<long double> SIGMAGG =
          Sigma22ggg(xp) * std::pow(LNbar, 2) +
          Sigma21ggg(xp) * LNbar; /* +Sigma20ggg(xp); */
      result += aass * LOgggH(N, xp) * (GOgggH(xp) + SIGMAGG);
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
