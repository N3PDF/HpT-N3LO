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
    std::complex<long double> S1,S2,S3,S4;
    std::complex<long double> N3,N4,N5,N6;
    std::complex<long double> SSTR2P,SSTR3P;
    std::complex<long double> N2,NM,NMS,N1S,N1T,N2S,N2T;
    std::complex<long double> NS,NT,NFO,NFI,NSI,NSE,NE,NN;
    std::complex<long double> S11,S12,S13,S14,S15,S16,S1M,S21,S31,S2M;
    std::complex<long double> SPMOM,SLC,SLV,SSCHLM,SSTR2M,SSTR3M,SSCHLP;
    std::complex<long double> NI,NI2,NI3,NMI, NMI2,N1I,N1I2,N1I3,N2I,N1;
};


class AnomDimensions
{
    public:
        AnomDimensions(void *params);
        virtual ~AnomDimensions();

        void sums(std::complex<long double> N);
        void ComputeGamma(std::complex<long double> N, int order);

        std::complex<long double> plus0,minus0;
        std::complex<long double> gg1,qq1,qg1,gq1,WW1,TT1,VV1;
        std::complex<long double> gg2,qq2,qg2,gq2,WW2,TT2,VV2;
        std::complex<long double> gg0,qq0,qg0,gq0,WW0,TT0,VV0;

    private:
        HSum HAP;
        gamma1sums g1s;

		int NC, NF;
        long double CA, CF;
        long double EulerGamma=0.5772156649;

        std::complex<long double> PNPA,PNSB,PNSC,PNMA,PPSA,PQGA;
        std::complex<long double> PQGB,PGQA,PGQB,PGQC,PGGA,PGGB,PGGC;
        std::complex<long double> P2PLSN,P2MINN,P2VALN,P2QGN,P2GQN,P2GGN,P2PSN;

        // LO anomalous dimensions
        std::complex<long double> gammagg0(std::complex<long double> N);
        std::complex<long double> gammaSg0(std::complex<long double> N);
        std::complex<long double> gammagS0(std::complex<long double> N);
        std::complex<long double> gammansplus0(std::complex<long double> N);

        // NLO important definitions
        void DEF1(std::complex<long double> N);

        // NLO anomalous dimensions
        std::complex<long double> gammagg1(std::complex<long double> N);
        std::complex<long double> gammaSg1(std::complex<long double> N);
        std::complex<long double> gammagS1(std::complex<long double> N);
        std::complex<long double> gammaps1(std::complex<long double> N);
        std::complex<long double> gammansplus1(std::complex<long double> N);
        std::complex<long double> gammansminus1(std::complex<long double> N);

        // NNLO important definitions
        void DEF2(std::complex<long double> N);

        // NNLO anomalous dimensions
        std::complex<long double> gammagg2(std::complex<long double> N);
        std::complex<long double> gammaSg2(std::complex<long double> N);
        std::complex<long double> gammagS2(std::complex<long double> N);
        std::complex<long double> gammaps2(std::complex<long double> N);
        std::complex<long double> gammansplus2(std::complex<long double> N);
        std::complex<long double> gammansminus2(std::complex<long double> N);
        std::complex<long double> gammansval2(std::complex<long double> N);
};
