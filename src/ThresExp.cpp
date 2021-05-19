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

#include "../include/ThresExp.h"


ThresExp::ThresExp(int order, int channel, void *params)
{
	PhysParams param = *reinterpret_cast<PhysParams*>(params); 


    ORD = order;
    CHANNEL = channel;
	NC = param.nc;
	NF = param.nf;

	CA = NC;
	CF = (NC*NC-1.)/(2.*NC);

    MH2  = std::pow(param.mh, 2);
    MUR2 = std::pow(param.mur, 2);
    MUF2 = std::pow(param.muf, 2);

 	LF = std::log(MH2/MUF2);
 	LR = std::log(MH2/MUR2);

    aass = param.alphas; // TODO; double-check PI normalization
    SROOT = param.sroot;
    SIGMA0 = param.sigma0;

    // Compute Beta Functions
    Beta0 = (11.*CA-2.*NF)/(12.*M_PI);
    Beta1 = ((17.*CA*CA-5.*CA*NF-3.*CF*NF)*2./3.)/(16.*M_PI*M_PI);
    Beta2 = ((2857./54.*CA*CA*CA+(CF*CF-205./18.*CF*CA-1415./54.*CA*CA)
            *NF+(11./9.*CF+79./54.*CA)*NF*NF))/std::pow(4.*M_PI,3);

    // One & too loop Cusp Anomalous Dimensions
	// In the following, the 1/PI is contained in the definitions
	// of the cusp Anomalous dimensions not in the normalization of as
    Bth1g = -Beta0;
    Ath1g = CA/M_PI;
    Ath1q = CF/M_PI;
    Bth1q = -3./4.*CF/M_PI;
    Ath2g = (CA/2.*(CA*(67./18.-zeta2)-5./9.*NF))/std::pow(M_PI,2);
    Ath2q = (CF/2.*(CA*(67./18.-zeta2)-5./9.*NF))/std::pow(M_PI,2);
}

ThresExp::~ThresExp(){}


//==========================================================================================//
//                                  LO pt-distributions                                     //
//------------------------------------------------------------------------------------------//

// gg->g
std::complex<double> ThresExp::LOgggH(std::complex<double> NN, double xp)
{
    std::complex<double> CLOgggH;
    // TODO: verify correspondence with small-pt
    std::complex<double> N = NN;
    std::complex<double> half(0.5,0.);
    std::complex<double> xprad(std::pow(std::sqrt(1.+xp)-std::sqrt(xp),4),0.);

    CLOgggH = 2. * aass * CA/std::sqrt(M_PI) * 1./xp*std::exp(LogGamma(N)-LogGamma(N+0.5))* \
    (Hyp2F1(half,N,N+0.5,xprad) -2. *(1.+xp)/(std::pow(std::sqrt(1.+xp)+std::sqrt(xp),2.))* \
    N/(N+0.5) *Hyp2F1(half,N+1.,N+1.5,xprad) +((1.+xp)*(3.+xp))/(std::pow(std::sqrt(1.+xp)+ \
    std::sqrt(xp),4.))*N*(N+1.)/((N+0.5)*(N+1.5))*Hyp2F1(half,N+2.,N+2.5,xprad)-2.*(1.+xp)/ \
    (std::pow(std::sqrt(1.+xp)+std::sqrt(xp),6.))*N*(N+1.)*(N+2.)/((N+0.5)*(N+1.5)*(N+2.5)) \
    *Hyp2F1(half,N+3.,N+3.5,xprad)+1./(std::pow(std::sqrt(1.+xp)+std::sqrt(xp),8))*N*(N+1.) \
    *(N+2.)*(N+3.)/((N+0.5)*(N+1.5)*(N+2.5)*(N+3.5))*Hyp2F1(half,N+4.,N+4.5,xprad));

    return CLOgggH;
}

// gq->q
std::complex<double> ThresExp::LOgqqH(std::complex<double> NN, double xp)
{
    std::complex<double> CLOgqqH;
    // TODO: verify correspondence with small-pt
    std::complex<double> N=NN;
    std::complex<double> half(0.5,0.);
    std::complex<double> xprad(std::pow(std::sqrt(1.+xp)-std::sqrt(xp),4.),0.);

    CLOgqqH = aass * CF/std::sqrt(M_PI) * 1./xp * std::exp(LogGamma(N)-LogGamma(N+0.5))* \
    (Hyp2F1(half,N,N+0.5,xprad) - (4.+3.*xp)/(std::pow(std::sqrt(1.+xp)+std::sqrt(xp),2.))* \
    N/(N+0.5) * Hyp2F1(half,N+1.,N+1.5,xprad) + 3. * (1.+xp)/(std::pow(std::sqrt(1.+xp)+ \
    std::sqrt(xp),4.)) * N * (N+1.)/((N+0.5) * (N+1.5)) * Hyp2F1(half,N+2.,N+2.5,xprad)-1./ \
    (std::pow(std::sqrt(1.+xp)+std::sqrt(xp),6.))*N*(N+1.)*(N+2.)/((N+0.5)*(N+1.5)*(N+2.5)) \
    *Hyp2F1(half,N+3.,N+3.5,xprad));

    return CLOgqqH;
}

// qq->g
std::complex<double> ThresExp::LOqqgH(std::complex<double> NN, double xp)
{
    std::complex<double> CLOqqgH;
    // TODO: verify correspondence with small-pt
    std::complex<double> N = NN;
    std::complex<double> half(0.5,0.);
    std::complex<double> xprad(std::pow(std::sqrt(1.+xp)-std::sqrt(xp),4.),0.);

    CLOqqgH = 2.*aass*CF*CF/std::sqrt(M_PI)*1./(std::pow(std::sqrt(1.+xp)+std::sqrt(xp),2.)) \
    * std::exp(LogGamma(N)-LogGamma(N+0.5)) * (Hyp2F1(half,N,N+0.5,xprad)-2.*(1.+xp)/ \
    (std::pow(std::sqrt(1.+xp)+std::sqrt(xp),2.))*N/(N+0.5)*Hyp2F1(half,N+1.,N+1.5,xprad) \
    + 1./(std::pow(std::sqrt(1.+xp)+std::sqrt(xp),4.)) * N * (N+1.)/((N+0.5)*(N+1.5))* \
    Hyp2F1(half,N+2.,N+2.5,xprad));

    return CLOqqgH;
}


//==========================================================================================//
//                  Sigma functions (factor of ln N-enhanced terms)                         //
//------------------------------------------------------------------------------------------//

// gg->g
double ThresExp::Sigma22ggg(double xp)
{
	return 3.*Ath1g/8.;
}

double ThresExp::Sigma21ggg(double xp)
{
	double xQp = std::sqrt(1+xp)+std::sqrt(xp);
	return -Bth1g/2.-Ath1g/2.*std::log(xp/xQp)-2.*Ath1g*LF;
}

double ThresExp::Sigma20ggg(double xp)
{
	return 3./2.*Ath1g*zeta2;
}


//==========================================================================================//
//                             Expanded Resummed Expression                                 //
//------------------------------------------------------------------------------------------//


// TODO: change this into a vector
std::complex<double> ThresExp::ThresExpExpr(std::complex<double> N, double xp)
{
	std::complex<double> zero(0.,0.);
	std::complex<double> result;

	std::complex<double> Nbar = N*std::exp(EulerGamma);
	std::complex<double> LNbar = 2.*aass*Beta0*std::log(Nbar);

	switch (ORD)
	{
		// TODO: Match CASES with HpT-MON
		case(0): // order as^0
		{
			if ((CHANNEL==0)||(CHANNEL==5)) return zero; // gg-channel or ALL
			if ((CHANNEL==1)||(CHANNEL==5)) return zero; // gq-channel or ALL
			if ((CHANNEL==2)||(CHANNEL==5)) return zero; // qq-channel or ALL
		}
		break;
        case(1): // order as^1
		{
			if ((CHANNEL==0)||(CHANNEL==5)) {result += LOgggH(N,xp);}
			if ((CHANNEL==1)||(CHANNEL==5)) {result += LOgqqH(N,xp);}
			if ((CHANNEL==2)||(CHANNEL==5)) {result += LOqqgH(N,xp);}
		}
		break;
        case(2): // order as^2
		{
			if ((CHANNEL==0)||(CHANNEL==5))
			{
				std::complex<double> SIGMAGG = Sigma22ggg(xp)*std::pow(LNbar,2) \
					+Sigma21ggg(xp)*LNbar+Sigma20ggg(xp);
				result += aass*LOgggH(N,xp)*SIGMAGG;
			}
			if ((CHANNEL==1)||(CHANNEL==5))
			{
				return zero;
			}
			if ((CHANNEL==2)||(CHANNEL==5))
			{
				return zero;
			}
		}
		break;
	}
}
