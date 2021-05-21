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

    MH2  = std::pow(param.mh, 2);
    MUR2 = std::pow(param.mur, 2);
    MUF2 = std::pow(param.muf, 2);

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
        


std::complex<double> CombinedRes::Matching(std::complex<double> N, double pt, int scheme)
{
 	// TODO: re-check definition MH2 vs. Qs2
 	double xp = std::pow(pt,2)/std::pow(MH2,2);
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
    


std::complex<double> CombinedRes::CombinedResExpr(std::complex<double> N, double pt, int scheme)
{
	double nn = N.real(); // take only real part. Does not work for complex

	std::vector<double> ExactMellin = MELLINPARTONIC->partonichiggsdpt(pt,nn);
	std::complex<double> SptMellin = SMALLPT->SmallptExpExpr(N,pt);
	std::complex<double> ThresMellin = THRESHOLD->ThresExpExpr(N,pt);

	std::complex<double> ExactMellinCmpx(ExactMellin[0],0.);
	return (1.-Matching(N,pt,scheme))*SptMellin+Matching(N,pt,scheme)*ThresMellin;
}