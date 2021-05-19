#pragma once

#include <cmath>
#include <stdio.h>
#include <complex>
#include <iostream>
#include <gsl/gsl_math.h>

#include "./HarmonicSum.h"
#include "./ComplexDefs.h"

#include "higgs-fo/params.h"

using namespace std;

struct gamma1sums
{
    std::complex<double> S1,S2,S3,S4;
    std::complex<double> N3,N4,N5,N6;
    std::complex<double> SSTR2P,SSTR3P;
    std::complex<double> N2,NM,NMS,N1S,N1T,N2S,N2T;
    std::complex<double> NS,NT,NFO,NFI,NSI,NSE,NE,NN;
    std::complex<double> S11,S12,S13,S14,S15,S16,S1M,S21,S31,S2M;
    std::complex<double> SPMOM,SLC,SLV,SSCHLM,SSTR2M,SSTR3M,SSCHLP;
    std::complex<double> NI,NI2,NI3,NMI, NMI2,N1I,N1I2,N1I3,N2I,N1;
};


class AnomDimensions
{
    public:
        AnomDimensions(void *params);
        virtual ~AnomDimensions();

        void sums(std::complex<double> N);
        void ComputeGamma(std::complex<double> N, int order);

        std::complex<double> plus0,minus0;
        std::complex<double> gg1,qq1,qg1,gq1,WW1,TT1,VV1;
        std::complex<double> gg2,qq2,qg2,gq2,WW2,TT2,VV2;
        std::complex<double> gg0,qq0,qg0,gq0,WW0,TT0,VV0;

    private:
        HSum HAP;
        gamma1sums g1s;

        double NF, NC, CA, CF;
        double EulerGamma=0.5772156649;

        std::complex<double> PNPA,PNSB,PNSC,PNMA,PPSA,PQGA;
        std::complex<double> PQGB,PGQA,PGQB,PGQC,PGGA,PGGB,PGGC;
        std::complex<double> P2PLSN,P2MINN,P2VALN,P2QGN,P2GQN,P2GGN,P2PSN;

        // LO anomalous dimensions
        std::complex<double> gammagg0(std::complex<double> N);
        std::complex<double> gammaSg0(std::complex<double> N);
        std::complex<double> gammagS0(std::complex<double> N);
        std::complex<double> gammansplus0(std::complex<double> N);

        // NLO important definitions
        void DEF1(std::complex<double> N);

        // NLO anomalous dimensions
        std::complex<double> gammagg1(std::complex<double> N);
        std::complex<double> gammaSg1(std::complex<double> N);
        std::complex<double> gammagS1(std::complex<double> N);
        std::complex<double> gammaps1(std::complex<double> N);
        std::complex<double> gammansplus1(std::complex<double> N);
        std::complex<double> gammansminus1(std::complex<double> N);

        // NNLO important definitions
        void DEF2(std::complex<double> N);

        // NNLO anomalous dimensions
        std::complex<double> gammagg2(std::complex<double> N);
        std::complex<double> gammaSg2(std::complex<double> N);
        std::complex<double> gammagS2(std::complex<double> N);
        std::complex<double> gammaps2(std::complex<double> N);
        std::complex<double> gammansplus2(std::complex<double> N);
        std::complex<double> gammansminus2(std::complex<double> N);
        std::complex<double> gammansval2(std::complex<double> N);
};
