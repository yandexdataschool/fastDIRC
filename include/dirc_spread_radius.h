#include "dirc_point.h"
#include "dirc_spread_relative.h"
#include <vector>
#ifndef DIRC_SPREAD_RADIUS
#define DIRC_SPREAD_RADIUS
class DircSpreadRadius : public DircSpreadRelative
{
public:
	//by reference - later for speed
	DircSpreadRadius(\
		std::vector<dirc_point> isupport,\
		bool itest_time_dir = true);
	
	float relative_spread_function(dirc_point vec);
	
	virtual float radius_spread_function(float r) = 0;
};
#endif
