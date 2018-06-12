#include "dirc_point.h"
#include "dirc_spread_radius.h"
#include <vector>
#ifndef DIRC_SPREAD_LINEAR_SOFT
#define DIRC_SPREAD_LINEAR_SOFT
class DircSpreadLinearSoft : public DircSpreadRadius
{
private:
	float lin_slope, r_trans, sigma2, max_val;
public:
	//by reference - later for speed
	DircSpreadLinearSoft(\
		float ilin_slope, \
		float ir_trans, \
		float isigma, \
		std::vector<dirc_point> isupport,\
		bool itest_time_dir = true);
	
	float radius_spread_function(float r);
};
#endif
