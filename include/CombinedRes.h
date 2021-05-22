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

#include <stdio.h>
#include <math.h>
#include <cmath>
#include <stdlib.h>
#include <string>
#include <vector>
#include <iostream>

#include "./ComplexDefs.h"
#include "./ThresExp.h"
#include "./SmallptExp.h"

#include "higgs-fo/partonic.h"

class CombinedRes{
    public:
        CombinedRes(int order, int channel, std::string pdfname, void *params);
        virtual ~CombinedRes();

        // Attribute that compute the expanded results
        std::complex<long double> CombinedResExpr(std::complex<long double> N, long double pt, int scheme);

    private:
        int EXACT_ORD;
        /* long double LF, LR, LQ; */
        long double QS2, MH2, MUF2, MUR2;
        /* long double aass, SIGMA0; */

        // Init. resummation classes
        SmallptExp *SMALLPT;
        ThresExp *THRESHOLD;

        // Init. exact FO class
        CrossHiggs *MELLINPARTONIC;


        // Matching function
        std::complex<long double> Matching(std::complex<long double> N, long double pt, int scheme);
};
