#include "../include/dirc_probability_separation.h"
#include <stdio.h>
#include <math.h>

// DircProbabilitySeparation::DircProbabilitySeparation(\
// 	DircProbabilitySpread* ineg_dens, \
// 	DircProbabilitySpread* ipos_dens)
DircProbabilitySeparation::DircProbabilitySeparation(\
	DircSpreadGaussian* ineg_dens, \
	DircSpreadGaussian* ipos_dens)
{
	neg_dens = ineg_dens;
	pos_dens = ipos_dens;
}

float DircProbabilitySeparation::spread_function(float neg, float pos)
{
	float rval;
	float pmn = pos-neg;
	float power = 1;
	if (pmn > 0)
	{
// 		rval = pmn*pmn;
// 		rval = 1;
// 		rval = sqrt(pmn);
		rval = pow(pmn,power);
	}
	else
	{
// 		rval = -pmn*pmn;
// 		rval = -1;
// 		rval = -sqrt(fabs(pmn));
		rval = -pow(fabs(pmn),power);
	}
// 	rval = pmn;
	return rval;
}
float DircProbabilitySeparation::get_log_likelihood_spread_diff(std::vector<dirc_point> inpoints)
{
	float rval = 0;
	float tpos = 0;
	float tneg = 0;
	for (unsigned int i = 0; i < inpoints.size(); i++)
	{
		tpos = pos_dens->get_single_likelihood(inpoints[i]);
		tneg = neg_dens->get_single_likelihood(inpoints[i]);
		
		if (tpos < 1e-6 || tneg < 1e-6)
		{
			continue;
		}
		tpos = log(tpos);
		tneg = log(tneg);

		rval += spread_function(tneg,tpos);
	}
	return rval;
}
