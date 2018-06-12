#include "dirc_point.h"
#include <vector>

#ifndef DIRC_PROBABILITY_SPREAD
#define DIRC_PROBABILITY_SPREAD

class DircProbabilitySpread
{
protected:
	std::vector<dirc_point> support_points;
	bool test_time_dir;
	float get_weight(dirc_point inpoint);

	
public:
	//by reference - later for speed
	DircProbabilitySpread(\
		std::vector<dirc_point> isupport,\
		bool itest_time_dir = true);
	
	inline float get_log_likelihood(std::vector<dirc_point> inpoints) __attribute__((always_inline));
	float get_single_log_likelihood(dirc_point inpoint);
	virtual float support_spread_function(dirc_point support, dirc_point test) = 0;
};
#endif
