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


#include <cmath>
#include <higgs-fo/params.h>

#include "../include/SmallptExp.h"


SmallptExp::SmallptExp(int order, int channel, void *params):AD(params)
{
	PhysParams param = *reinterpret_cast<PhysParams*>(params); 

    ORD = order;
    CHANNEL = channel;
	NC = static_cast<long double>(param.nc);
	NF = static_cast<long double>(param.nf);

	CA = NC;
	CF = (NC*NC-1.)/(2.*NC);

    MH2  = static_cast<long double>(std::pow(param.mh, 2));
    MUR2 = static_cast<long double>(std::pow(param.mur, 2));
    MUF2 = static_cast<long double>(std::pow(param.muf, 2));

 	// TODO: Define Q globaly and independent of MH2
	QS2 = MH2;

 	LQ = std::log(MH2/QS2);
 	LF = std::log(MH2/MUF2);
 	LR = std::log(MH2/MUR2);

    aass   = static_cast<long double>(param.alphas/M_PIl); // TODO; long double-check PI normalization
    SROOT  = static_cast<long double>(param.sroot);
    SIGMA0 = static_cast<long double>(param.sigma0);

    // Compute Beta Functions
    Beta0 = (11.*CA-2.*NF)/12;
    Beta1 = (17.*CA*CA-5.*CA*NF-3.*CF*NF)/24.;
    Beta2 = ((2857./54.*CA*CA*CA+(CF*CF-205./18.*CF*CA-1415./54.*CA*CA)
            *NF+(11./9.*CF+79./54.*CA)*NF*NF))/64.;

	// Cusp Anomlaous Dimensions
	// Here, we follow the convention in Bozzi's paper, i.e. the
	// 1/PI is not contained in the definition of the cusp anomalous
	// dimensions but rather in the normalization of as
	Apt1g = CA;
	Bpt1g = -1./6.*(11.*CA-2.*NF);
	Apt2g = CA/2.*((67./18.-std::pow(M_PIl,2)/6.)*CA-5./9.*NF);
    Bpt2g = 0.; // TODO: Implemente correct expression here
}

SmallptExp::~SmallptExp(){}


//==========================================================================================//
//                            Fourier Inverse of ln^k b                                     //
//------------------------------------------------------------------------------------------//

long double SmallptExp::LC1(long double xp)
{
	return -1./xp;
}

long double SmallptExp::LC2(long double xp)
{
	return -2./xp*std::log(1./xp);
}

long double SmallptExp::LC3(long double xp)
{
	return -3./xp*std::pow(std::log(1./xp),2);
}

long double SmallptExp::LC4(long double xp)
{
	return -4./xp*(std::pow(std::log(1./xp),3)-4.*zeta3);
}


//==========================================================================================//
//                                Coefficient functions                                     //
//------------------------------------------------------------------------------------------//

std::complex<long double> SmallptExp::C1GG(std::complex<long double> N)
{
    std::complex<long double> zero(0.,0.);
    return zero;
}

std::complex<long double> SmallptExp::C1GQ(std::complex<long double> N)
{
    return CF/2./(N+1);
}


//==========================================================================================//
//                            Sigma pt-enhanced functions                                   //
//------------------------------------------------------------------------------------------//

std::complex<long double> SmallptExp::SmallptExpExpr(std::complex<long double> N, long double pt)
{
    int pc = 2; // TODO: long double-check this
	AD.ComputeGamma(N,1); // Init. Anomalous Dimensions
    std::complex<long double> zero(0.,0.);
	std::complex<long double> ones(1.,0.);
    std::complex<long double> result(0.,0.);

    // Hard functions
    long double h1gg = (11.+3.*M_PIl*M_PIl)/2.;

 	// TODO: re-check definition MH2 vs. Qs2
 	long double xp = std::pow(pt,2)/std::pow(MH2,2);

	switch (ORD)
	{
        case(0): // order as^0
		{
			if ((CHANNEL==0)||(CHANNEL==4)) {result += ones;}
			if ((CHANNEL==1)||(CHANNEL==4)) {result += zero;}
			if ((CHANNEL==2)||(CHANNEL==4)) {result += zero;}
		}
		break;
        case(1): // order as^1
        {
            if ((CHANNEL==0)||(CHANNEL==4)) // gg-channel or ALL
            {
                // pt-enhanced Sigma terms
                std::complex<long double> Sigma12gg = -Apt1g/2.;
                std::complex<long double> Sigma11gg = -Bpt1g-2*AD.gg0-Apt1g*LQ;

                // constant-terms when pt->0
                std::complex<long double> HH1GG = h1gg+2.*C1GG(N)+2.*(LF-LQ)*AD.gg0 \
                    -LQ*(Bpt1g+(Bpt1g*LQ)/2.)-pc*LR*Beta0;

                result += aass*(Sigma12gg*LC2(xp)+Sigma11gg*LC1(xp)+HH1GG);
            }
			if ((CHANNEL==1)||(CHANNEL==4))
            {
                result += zero;
            }
			if ((CHANNEL==2)||(CHANNEL==4))
            {
                result += zero;
            }
        }
        break;
        case(2):
        {
            if ((CHANNEL==0)||(CHANNEL==4)) // gg-channel or ALL
            {
                // pt-enhanced Sigma terms
                std::complex<long double> Sigma24gg = std::pow(Apt1g,2)/8.;
                std::complex<long double> Sigma23gg = -Apt1g*(Beta0/3.+1./2. \
                    *(-Bpt1g-Apt1g*LQ-2.*AD.gg0));
                std::complex<long double> Sigma22gg = (-Apt2g-(Bpt1g+Apt1g*LQ-Beta0) \
                    *(-Bpt1g-Apt1g-2.*AD.gg0))/2.-Apt1g/2.*(h1gg-LQ*(Bpt1g+(Bpt1g \
                    *LQ)/2.)-(-LQ+LR)*Beta0-pc*LR*Beta0+2.*C1GG(N)+2.*(LF-LQ)*AD.gg0) \
                    +1./2.*(-2.*(-Bpt1g-Apt1g*LQ-2.*AD.gg0)*AD.gg0+2.*AD.gq0*AD.qg0);
                std::complex<long double> Sigma21gg = -Bpt2g-Apt2g*LQ+(LQ-LR)*Beta0 \
                    *(-Bpt1g-Apt1g-2.*AD.gg0)-(Bpt1g+Apt1g*LQ+2.*AD.gg0)*(h1gg-LQ \
                    *(Bpt1g+(Bpt1g*LQ)/2.)-pc*LR*Beta0+2.*C1GG(N)+2.*(LF-LQ)*AD.gg0) \
                    -2.*(C1GQ(N)+(LF-LQ)*AD.gq0)*AD.qg0+Beta0*(2.*C1GG(N)+2.*AD.gg1);

                // constant terms when pt->0
                long double h2gg = 0.;
                std::complex<long double> HH2GG = 2.*h2gg+zero;

                result += aass*aass*(Sigma24gg*LC2(xp)+Sigma23gg*LC3(xp)+Sigma22gg \
                    *LC2(xp)+Sigma21gg*LC1(xp)+HH2GG);
            }
			if ((CHANNEL==1)||(CHANNEL==4))
            {
                result += zero;
            }
			if ((CHANNEL==2)||(CHANNEL==4))
            {
                result += zero;
            }
        }
        break;
	}
	
    return 2.*pt/MH2*SIGMA0*result;
}
