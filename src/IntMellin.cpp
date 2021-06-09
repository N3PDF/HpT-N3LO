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

#include "../include/IntMellin.h"

#include <higgs-fo/params.h>

MellinTrans::MellinTrans(int order, int channel, std::string pdfname,
                         void *params)
    : xThres(order, channel, params) {
  PhysParams param = *reinterpret_cast<PhysParams *>(params);

  MH2 = static_cast<long double>(std::pow(param.mh, 2));
  MUR2 = static_cast<long double>(std::pow(param.mur, 2));
  MUF2 = static_cast<long double>(std::pow(param.muf, 2));
}

MellinTrans::~MellinTrans() {}

struct IntStruct {
  MellinTrans *TransMellin;
  double N;
  double pt;
};

double MellinTrans::ExtractThresMom(double x, double N, double pt) {
  return xThres.ThresXspaceExpr(x, N, pt);
}

double xThresIntegrand(double x, void *p) {
  IntStruct par = *reinterpret_cast<IntStruct *>(p);
  return par.TransMellin->ExtractThresMom(x, par.N, par.pt);
}

std::complex<long double> MellinTrans::xSpaceThres(std::complex<long double> N,
                                                   long double pt) {
  int method = 0;
  double reslt, error;
  double nn = static_cast<double>(N.real());

  // pass parameters
  IntStruct finalparams;
  finalparams.N = nn;
  finalparams.pt = static_cast<double>(pt);
  finalparams.TransMellin = this;

  // Integration
  reslt = ExtrIntegration::IntegrateOverx(method, xThresIntegrand, &finalparams,
                                          &error);

  std::complex<long double> finreslt(static_cast<long double>(reslt), 0.);
  return finreslt;
}
