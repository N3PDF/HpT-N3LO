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
#include "./Integration.h"
#include "./ThresXspace.h"

class MellinTrans {
 public:
  MellinTrans(int order, int channel, std::string pdfname, void *params);
  virtual ~MellinTrans();

  double ExtractThresMom(double x, double N, double pt);
  std::complex<long double> xSpaceThres(std::complex<long double> N,
										long double pt);

 private:
  int NF, ORD, CHANNEL;
  double MH2, MUR2, MUF2;

  ThresXspace xThres;
};
