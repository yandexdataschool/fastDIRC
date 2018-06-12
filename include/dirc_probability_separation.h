#include "dirc_point.h"
#include "dirc_probability_spread.h"
#include "dirc_spread_gaussian.h"
#include <vector>

#ifndef DIRC_PROBABILITY_SPERATION
#define DIRC_PROBABILITY_SPERATION

class DircProbabilitySeparation
{
protected:
// 	DircProbabilitySpread* neg_dens;
// 	DircProbabilitySpread* pos_dens;
	DircSpreadGaussian* neg_dens;
	DircSpreadGaussian* pos_dens;
public:
// 	DircProbabilitySeparation(DircProbabilitySpread* ineg_dens, DircProbabilitySpread* ipos_dens);
	DircProbabilitySeparation(DircSpreadGaussian* ineg_dens, DircSpreadGaussian* ipos_dens);
  
  
	//make virtual for polymorphism later
	virtual float spread_function(float neg, float pos);
	float get_log_likelihood_spread_diff(std::vector<dirc_point> inpoints);
};

#endif