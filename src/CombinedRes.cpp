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

#include "../include/CombinedRes.h"


CombinedRes::CombinedRes(int order, int channel, std::string pdfname, void *params)
{
	PhysParams param = *reinterpret_cast<PhysParams*>(params); 

    MH2  = static_cast<long double>(std::pow(param.mh, 2));
    MUR2 = static_cast<long double>(std::pow(param.mur, 2));
    MUF2 = static_cast<long double>(std::pow(param.muf, 2));

    EXACT_ORD = order-1;

 	SMALLPT = new SmallptExp(order,channel,params);
 	THRESHOLD = new ThresExp(order,channel,params);
 	MELLINPARTONIC = new CrossHiggs(EXACT_ORD,channel,pdfname,params);

}

CombinedRes::~CombinedRes()
{
	delete MELLINPARTONIC;
	delete SMALLPT;
	delete THRESHOLD;
}
        


std::complex<long double> CombinedRes::Matching(std::complex<long double> N, long double pt, int scheme)
{
 	// TODO: re-check definition MH2 vs. Qs2
 	long double xp = std::pow(pt,2)/std::pow(MH2,2);
	if (scheme==0) 			// small-pt only
	{
		return (0.);
	} else if (scheme==1) 	// threshold only
	{
		return (1.);
	} else  				// combined
	{
		return (std::pow(N,3)*std::pow(xp,2)/(1.+std::pow(N,3)*std::pow(xp,2)));
	}
}
    


std::complex<long double> CombinedRes::CombinedResExpr(std::complex<long double> N, long double pt, int scheme)
{
	double pp = static_cast<double>(pt);
	double nn = static_cast<double>(N.real()); // take only real part. Does not work for complex

	std::vector<double> ResultsMellin = MELLINPARTONIC->partonichiggsdpt(pp,nn);
	std::vector<long double> ExactMellin(ResultsMellin.begin(), ResultsMellin.end());

	std::complex<long double> SptMellin = SMALLPT->SmallptExpExpr(N,pt);
	std::complex<long double> ThresMellin = THRESHOLD->ThresExpExpr(N,pt);

	std::complex<long double> ExactMellinCmpx(ExactMellin[0],0.);
	return (1.-Matching(N,pt,scheme))*SptMellin+Matching(N,pt,scheme)*ThresMellin;
}
