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

#include "../include/Integration.h"

struct IntData {
  std::function<double(double, void *)> fun;
  void *param;
};


int xIntegrand(int *ndim, double *x, int *ncomp, double *y, void *p) {
  IntData inparam = *reinterpret_cast<IntData *>(p);
  y[0] = inparam.fun(x[0], inparam.param);
  return 0;
}


double ExtrIntegration::IntegrateOverx(int method, double(Func)(double, void *),
                                       void *pp, double *error) {
  int fail;
  double prec = 1e-8;
  double *res = NULL;
  double *err = NULL;
  double *prb = NULL;

  res = new double[1];
  err = new double[1];
  prb = new double[1];

  int ndim = 2;
  int ncomp = 1;
  int verbose = 0;

  IntData UserData;
  UserData.fun = Func;
  UserData.param = pp;

  integration::CubaIntegrator(method, xIntegrand, ndim, ncomp, prec, &res,
                              &fail, &err, &prb, &UserData, verbose);

  *error = err[0];
  return res[0];
}
