#include "dirc_point.h"
#include "dirc_probability_spread.h"
#include <vector>
#ifndef DIRC_SPREAD_RELATIVE
#define DIRC_SPREAD_RELATIVE
class DircSpreadRelative : public DircProbabilitySpread
{
public:
	//by reference - later for speed
	DircSpreadRelative(\
		std::vector<dirc_point> isupport,\
		bool itest_time_dir = true);
	
	float support_spread_function(dirc_point support, dirc_point test);
	
	virtual float relative_spread_function(dirc_point vec) = 0;
};
#endif
