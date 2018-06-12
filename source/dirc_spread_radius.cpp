#include "../include/dirc_point.h"
#include "../include/dirc_spread_relative.h"
#include "../include/dirc_spread_radius.h"
#include <vector>
#include <math.h>

DircSpreadRadius::DircSpreadRadius(\
	std::vector<dirc_point> isupport,\
	bool itest_time_dir /*=true*/)\
	: DircSpreadRelative(isupport,itest_time_dir) {}
	
float DircSpreadRadius::relative_spread_function(dirc_point vec)
{
	return radius_spread_function(sqrt(vec.x*vec.x + vec.y*vec.y));
// 	return radius_spread_function(fabs(vec.x) + fabs(vec.y));
}
