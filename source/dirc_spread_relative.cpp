#include "../include/dirc_spread_relative.h"
#include "../include/dirc_probability_spread.h"
#include "../include/dirc_point.h"
#include <vector>

DircSpreadRelative::DircSpreadRelative(\
	std::vector<dirc_point> isupport,\
	bool itest_time_dir /*=true*/)\
	: DircProbabilitySpread(isupport, itest_time_dir) {}
	
float DircSpreadRelative::support_spread_function(dirc_point support, dirc_point test)
{
	dirc_point ovec;
	
	ovec.x = test.x - support.x;
	ovec.y = test.y - support.y;
	
	return relative_spread_function(ovec);
}