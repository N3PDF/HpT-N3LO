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
#include <gsl/gsl_sf_psi.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include "higgs-fo/params.h"

class ThresXspace {
 public:
  ThresXspace(int order, int channel, void *params);
  virtual ~ThresXspace();

  double ThresXspaceExpr(double x, double N, double pt);

  // LO pt-distribution
  double LOgggH(double x, double pt);
  double LOgqqH(double x, double pt);
  double LOqqgH(double x, double pt);

 private:
  int NC, NF, CA, ORD, CHANNEL;
  double LF, LR, LQ;
  double MH2, MUF2, MUR2, SROOT;
  double CF, aass, SIGMA0;

  double Beta0;

  // Polygamma functions
  double psi0one = gsl_sf_psi_n(0, 1);
  double psi1one = gsl_sf_psi_n(1, 1);
  double psi0half = gsl_sf_psi_n(0, 0.5);
  double psi1half = gsl_sf_psi_n(1, 0.5);
};
