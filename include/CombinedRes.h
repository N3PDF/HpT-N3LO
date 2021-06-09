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

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include "./ComplexDefs.h"
#include "./SmallptExp.h"
#include "./ThresExp.h"
<<<<<<< HEAD
#include "./HighEnergyExp.h"
=======
#include "./IntMellin.h"
>>>>>>> f00a03b (Mellin transform of the x-space threshold resummation)
#include "higgs-fo/partonic.h"

class CombinedRes {
 public:
  CombinedRes(int order, int channel, std::string pdfname, void *params);
  virtual ~CombinedRes();

  // Attribute that compute the expanded results
  std::complex<long double> CombinedResExpr(std::complex<long double> N,
                                            long double pt, int scheme);

 private:
  int ORD, EXACT_ORD;
  /* long double LF, LR, LQ; */
  long double QS2, MH2, MUF2, MUR2;
  /* long double aass, SIGMA0; */

  // Init. resummation classes
  SmallptExp *SMALLPT;
  ThresExp *THRESHOLD;
  HighEnergyExp *HIGHENERGY;
  MellinTrans *MELLIN;

  // Init. exact FO class
  CrossHiggs *MELLINPARTONIC;

  // Matching function
  std::complex<long double> Matching(std::complex<long double> N,
                                     long double pt, int scheme);
};
