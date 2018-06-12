#include "../include/dirc_point.h"
#include "../include/dirc_spread_linear_soft.h"
#include <vector>
#include <math.h>

DircSpreadLinearSoft::DircSpreadLinearSoft(\
	float ilin_slope, \
	float ir_trans, \
	float isigma, \
	std::vector<dirc_point> isupport,\
	bool itest_time_dir /*=true*/)\
	: DircSpreadRadius(isupport,itest_time_dir)
{
	lin_slope = ilin_slope;
	r_trans = ir_trans;
	sigma2 = isigma*isigma;
	max_val = r_trans*lin_slope + 1;
}
	
float DircSpreadLinearSoft::radius_spread_function(float r)
{
	if (r < r_trans)
	{
		return max_val - lin_slope*r;
	}
	else
	{
		float dev2 = r - r_trans;
		dev2 *= dev2;
		
		return exp(-dev2/sigma2);
// 		return 0;
		
	}
}
