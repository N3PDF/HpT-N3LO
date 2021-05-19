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
#include <gsl/gsl_math.h>

#include "./AnomalousDim.h"
#include "./MellinFunc.h"
#include "./ComplexDefs.h"

#include "higgs-fo/params.h"

class SmallptExp{
    public:
        SmallptExp(int order, int channel, void *params);
        virtual ~SmallptExp();

        // Attribute that compute the expanded results
        std::complex<double> SmallptExpExpr(std::complex<double> N, double xp);

    private:
        AnomDimensions AD;
        int NC, NF, CA, ORD, CHANNEL;
        double LF, LR, LQ;
        double QS2, MH2, MUF2, MUR2, SROOT;
        double CF, aass, SIGMA0;


        // Beta functions
        double Beta0;
        double Beta1;
        double Beta2;

        // Cusp Anomalous Dimensions
        double Apt1g;
        double Apt2g;
        double Bpt1g;
		double Bpt2g;


        // Fourier Inverse of ln b
        double LC1(double xp);
        double LC2(double xp);
        double LC3(double xp);
        double LC4(double xp);

        // Ceofficient functions
		// LO
		std::complex<double> C1GG(std::complex<double> N);
		std::complex<double> C1GQ(std::complex<double> N);


        // Zeta functions
        const double zeta2=gsl_sf_zeta_int(2);
        const double zeta3=gsl_sf_zeta_int(3);
        const double zeta4=gsl_sf_zeta_int(4);


        // EulerGamma
        // TODO: Move this to global
        const double EulerGamma=0.5772156649;
};
