#include "dirc_point.h"
#include "../include/dirc_probability_spread.h"
#include <math.h>
#include <stdio.h>
//by reference - later for speed
DircProbabilitySpread::DircProbabilitySpread(\
	std::vector<dirc_point> isupport, \
	bool itest_time_dir /*=true*/)
{
	support_points = isupport;
	test_time_dir = itest_time_dir;
}
float DircProbabilitySpread::get_log_likelihood(std::vector<dirc_point> inpoints)
{
	float tprob = 0;
	float rval = 0;
	int eval_count = 0;
	float log_mult = 4;
	float weight = 1;
	for (unsigned int i = 0; i < inpoints.size(); i++)
	{
		tprob = 0;
		weight = get_weight(inpoints[i]);
		for (unsigned int j = 0; j < support_points.size(); j++)
		{
			if (test_time_dir && support_points[j].t*inpoints[i].t < -.00001)
			{
				//Only look at if they went up or down for now
				continue;
			}
			tprob += support_spread_function(support_points[j],inpoints[i]);
			eval_count++;
		}
		tprob /= support_points.size();
		
		rval += weight*log_mult*log(tprob);
	}
// 	printf("eval_count: %d\n",eval_count);
	rval -= log(inpoints.size());
	
	return rval;
}
float DircProbabilitySpread::get_single_log_likelihood(dirc_point inpoint)
{
	float tprob = 0;
	float log_mult = 4;
	float weight = 1;
	weight = get_weight(inpoint);
	for (unsigned int j = 0; j < support_points.size(); j++)
	{
		if (test_time_dir && support_points[j].t*inpoint.t < -.00001)
		{
			//Only look at if they went up or down for now
			continue;
		}
		tprob += support_spread_function(support_points[j],inpoint);
	}
	tprob /= support_points.size();
	tprob = weight*log_mult*log(tprob + 10e-3);
	
	return tprob;
}
	
float DircProbabilitySpread::get_weight(dirc_point inpoint)
{
	return 1;
// 	return 1 + sqrt(sqrt(fabs(inpoint.x)));
}