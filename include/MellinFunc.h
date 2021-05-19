#pragma once

#include <cmath>
#include <iostream>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_zeta.h>
#include <gsl/gsl_sf_log.h>

#include "HarmonicSum.h"
#include "ComplexDefs.h"

using namespace std;

class MellinFunc
{
    public:
        MellinFunc();
        virtual ~MellinFunc();

        std::complex<double> D0(std::complex<double> N);                  // [1/(1-z)]_+
        std::complex<double> D1(std::complex<double> N);                  // [Log(1-z)/(1-z)]_+
        std::complex<double> D2(std::complex<double> N);                  // [Log(1-z)^2/(1-z)]_+
        std::complex<double> D3(std::complex<double> N);                  // [Log(1-z)^3/(1-z)]_+

        std::complex<double> Li2z(std::complex<double> N);                // Li2(z)
        std::complex<double> Li2mz(std::complex<double> N);               // Li2(-z)
        std::complex<double> S12z(std::complex<double> N);                // S12(z)
        std::complex<double> S12mz(std::complex<double> N);               // S12(mz)
        std::complex<double> S12z2(std::complex<double> N);               // S12(z^2)
        std::complex<double> Li3z(std::complex<double> N);                // Li3(z)
        std::complex<double> Li3mz(std::complex<double> N);               // Li3(mz)
        std::complex<double> Li2zLogz(std::complex<double> N);            // Li2(z)Log(z)
        std::complex<double> Logminus(std::complex<double> N);            // Log(1-z)
        std::complex<double> Logplus(std::complex<double> N);             // Log(1+z)
        std::complex<double> Logz(std::complex<double> N);                // Log(z)
        std::complex<double> Logz2(std::complex<double> N);               // Log(z)^2
        std::complex<double> Logz3(std::complex<double> N);               // Log(z)^3
        std::complex<double> LogzLogminus(std::complex<double> N);        // Log(z)Log(1-z)
        std::complex<double> LogzLogminus2(std::complex<double> N);       // Log(z)Log(1-z)^2
        std::complex<double> Logz2Logminus(std::complex<double> N);       // Log(z)^2Log(1-z)
        std::complex<double> Logz2minus(std::complex<double> N);          // log(z)^2/(1-z)
        std::complex<double> Logminus2(std::complex<double> N);           // Log(1-z)^2
        std::complex<double> Logminus3(std::complex<double> N);           // Log(1-z)^3
        std::complex<double> LogzLogminusminus(std::complex<double> N);   // Log(z)Log(1-z)/(1-z)
        std::complex<double> Logzminus(std::complex<double> N);           // log(z)/(1-z)
        std::complex<double> Logzminusplus(std::complex<double> N);       // log(z)/((1-z)(1+z))
        std::complex<double> plus(std::complex<double> N);                // 1/(1+z)
        std::complex<double> Li2minusminus(std::complex<double> N);       // Li2(1-z)/(1-z)
        std::complex<double> S12zregminus(std::complex<double> N);        // (S12(z)-Zeta(3))/(1-z)
        std::complex<double> Li2zregminus(std::complex<double> N);        // (Li2(z)-Zeta(2))/(1-z)
        std::complex<double> Li3zregminus(std::complex<double> N);        // (Li3(z)-Zeta(3))/(1-z)
        std::complex<double> Li2zregLogminusminus(std::complex<double> N); //(Li2(z)-Zeta(2))Log(1-z)/(1-z)
        std::complex<double> Li3zplus(std::complex<double> N);            // Li3(z)/(1+z)
        std::complex<double> S12z2plus(std::complex<double> N);           // S12(z^2)/(1+z)
        std::complex<double> S12mzplus(std::complex<double> N);           // S12(-z)/(1+z)
        std::complex<double> S12zplus(std::complex<double> N);            // S12(z)/(1+z)
        std::complex<double> Li3mzplus(std::complex<double> N);           // Li3(-z)/(1+z)
        std::complex<double> Li2zLogminus(std::complex<double> N);        // Li2(z)Log(1-z)
        std::complex<double> Li2zLogplus(std::complex<double> N);         // Li2(z)Log(1+z)
        std::complex<double> Li2mzLogplus(std::complex<double> N);        // Li2(-z)Log(1+z)
        std::complex<double> Li2mzLogz(std::complex<double> N);           // Li2(-z)Log(z)
        std::complex<double> Li2zLogzplusminus(std::complex<double> N);   // Li2(z)Log(z)/((1+z)(1-z))
        std::complex<double> Li2mzLogzplus(std::complex<double> N);       // Li2(-z)Log(z)/(1+z)
        std::complex<double> Li2mzLogminus2plus(std::complex<double> N);  // Li2(-z)Log[1-z]^2/(1+z)
        std::complex<double> Li2mzLogplusplus(std::complex<double> N);    // Li2(-z)Log[1+z]/(1+z)
        std::complex<double> Logz2Logminusminus(std::complex<double> N);  // Log(z)^2Log(1-z)/(1-z)
        std::complex<double> LogzLogminus2minus(std::complex<double> N);  // Log(z)Log(1-z)^2/(1-z)
        std::complex<double> Logz3minusplus(std::complex<double> N);      //  Log(z)^3 /((1-z)(1+z))
        std::complex<double> Logplusplus(std::complex<double> N);         // Log(1+z)/(1+z)
        std::complex<double> Li2zLogplusplus(std::complex<double> N);     // Li2(z)Log(1+z)/(1+z)
        std::complex<double> Logz2Logplusplus(std::complex<double> N);    // Log(1+z)Log(z)^2/(1+z)
        std::complex<double> LogzLogplus(std::complex<double> N);         // Log(z)Log(1+z)
        std::complex<double> Logz2Logplus(std::complex<double> N);        // Log(z)^2Log(1+z)
        std::complex<double> LogzLogplus2(std::complex<double> N);        // Log(z)Log(1+z)^2
        std::complex<double> LogzLogplus2plus(std::complex<double> N);    // Log(z)Log(1+z)^2/(1+z)

        std::complex<double> Logplus3plus(std::complex<double > N);       // Log(1+z)^3/(1+z)
        std::complex<double> Li3zregminusplus(std::complex<double> N);    // (Li3-Zeta[3])/(1-z)/(1+z)
        std::complex<double> Li3zoverplusplus(std::complex<double> N);    // (Li3(z/(1+z))/(1+z)
        std::complex<double> Logplus3(std::complex<double> N);            // Log[1+z]^3
        std::complex<double> Li2z2(std::complex<double> N);               // Li2(z^2)
        std::complex<double> Li2z2Logz(std::complex<double> N);           // Li2(z^2) Log(z)
        std::complex<double> Li3z2(std::complex<double> N);               // Li3(z^2)
        std::complex<double> Li3overplus(std::complex<double> N);         // Li3(1/(1+z))
        std::complex<double> Li2minus(std::complex<double> N);            // Li2(1-z)

    private:
        HSum H;
        double zeta2;
        double zeta3;
        double Li4;
        double log2;
        double log2q;
        double log2c;
        double zeta2q;
        double EulerGamma;
};
