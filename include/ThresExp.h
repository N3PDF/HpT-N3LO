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


#include "./MellinFunc.h"
#include "./ComplexDefs.h"

#include "higgs-fo/params.h"

class ThresExp{
    public:
        ThresExp(int order, int channel, void *params);
        virtual ~ThresExp();

        // Attribute that compute the expanded results
        /* std::vector<std::complex<long double>> ThresExpExpr( */
        /*     std::complex<long double> N, long double xp */
        /* ); */
        std::complex<long double> ThresExpExpr(std::complex<long double> N, long double pt);

        // Matching Coefficients
        long double Hth1gggH(long double xp);
        long double Hth1gqqH(long double xp);
        long double Hth1qqgH(long double xp);

        // LO pt-distribution
        std::complex<long double> LOgggH(std::complex<long double> N, long double xp);
        std::complex<long double> LOgqqH(std::complex<long double> N, long double xp);
        std::complex<long double> LOqqgH(std::complex<long double> N, long double xp);

    private:
        int NC, NF, CA, ORD, CHANNEL;
        long double LF, LR, LQ;
        long double MH2, MUF2, MUR2, SROOT;
        long double CF, aass, SIGMA0;


        // Beta functions
        long double Beta0;
        long double Beta1;
        long double Beta2;

        // Cusp Anomalous Dimensions
        long double Ath1g;
        long double Ath2g;
        long double Ath1q;
        long double Ath2q;
        long double Bth1g;
        long double Bth1q;

        // Sigma functions
        // gg->g
        long double Sigma22ggg(long double xp);
        long double Sigma21ggg(long double xp);
        long double Sigma20ggg(long double xp);
        // gq->g
        long double Sigma22gqg(long double xp);
        long double Sigma21gqg(long double xp);
        long double Sigma20gqg(long double xp);
        // qq->g
        long double Sigma22qqg(long double xp);
        long double Sigma21qqg(long double xp);
        long double Sigma20qqg(long double xp);

        // Zeta functions
        const long double zeta2=gsl_sf_zeta_int(2);
        const long double zeta3=gsl_sf_zeta_int(3);
        const long double zeta4=gsl_sf_zeta_int(4);


        // EulerGamma
        // TODO: Move this to global
        const long double EulerGamma=0.5772156649;
};
