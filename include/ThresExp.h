#pragma once

#include <stdio.h>
#include <math.h>
#include <cmath>
#include <stdlib.h>
#include <string>
#include <vector>
#include <iostream>
#include <gsl/gsl_math.h>


#include "./MellinFunc.h"
#include "./ComplexDefs.h"

#include "higgs-fo/params.h"

class ThresExp{
    public:
        ThresExp(int order, int channel, void *params);
        virtual ~ThresExp();

        // Attribute that compute the expanded results
        std::vector<std::complex<double>> ThresExpExpr(
            std::complex<double> N, double xp
        );

        // Matching Coefficients
        double Hth1gggH(double xp);
        double Hth1gqqH(double xp);
        double Hth1qqgH(double xp);

        // LO pt-distribution
        std::complex<double> LOgggH(std::complex<double> N, double xp);
        std::complex<double> LOgqqH(std::complex<double> N, double xp);
        std::complex<double> LOqqgH(std::complex<double> N, double xp);

    private:
        int NC, NF, CA, ORD, CHANNEL;
        double MH2, MUF2, MUR2, SROOT;
        double CF, aass, SIGMA0;


        // Beta functions
        double Beta0;
        double Beta1;
        double Beta2;

        // Cusp Anomalous Dimensions
        double Ath1g;
        double Ath2g;
        double Ath1q;
        double Ath2q;
        double Bth1g;
        double Bth1q;

        // Zeta functions
        const double zeta2=gsl_sf_zeta_int(2);
        const double zeta3=gsl_sf_zeta_int(3);
        const double zeta4=gsl_sf_zeta_int(4);


        // EulerGamma
        // TODO: Move this to global
        const double EulerGamma=0.5772156649;
};
