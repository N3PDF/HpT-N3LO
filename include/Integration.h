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
#include "higgs-fo/integration.h"

namespace ExtrIntegration {
double IntegrateOverx(int method, double(Func)(double, void *), void *pp,
                      double *error);
}
