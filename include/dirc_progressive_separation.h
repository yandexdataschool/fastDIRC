#include "dirc_optical_sim.h"
#include "dirc_point.h"
#include "dirc_spread_gaussian.h"

#ifndef DIRC_PROGRESSIVE_SEPARATION
#define DIRC_PROGRESSIVE_SEPARATION

class DircProgressiveSeparation
{
protected:
	int max_sim_phots;
	int step_sim_phots;
	float ll_threshold;
	float E,x,y,phi,theta;
	float mass_1,mass_2;
	float BAR, tracking_unc, ckov_unc;
	float fudge_sigma;
	
	DircOpticalSim* dirc_model;
	DircSpreadGaussian* spread_func;
	
	float ll_diff(\
		std::vector<dirc_point> &hit_points, \
		int num_support,\
		float beta_1,\
		float beta_2);
	float ll_pts(\
		std::vector<dirc_point> &hit_points, \
		int num_support,\
		float beta);
	
	void get_inf_and_ll(\
		int &num_neg_inf,\
		float &ll_sum,\
		std::vector<float> prob_vals);
public:
	DircProgressiveSeparation(\
		DircOpticalSim* imodel,\
		int imax_phots,\
		int istep_phots,\
		float isigma, \
		float x_unc,\
		float y_unc,\
		float t_unc,\
		float im1 /*= .4937*/,\
		float im2 /*= .1396*/,\
		float ithresh /*= 20*/);
	
	void set_masses(float im1,float im2);
	void set_threshold(float ithresh);
	void set_max_step_phots(int im, int is);
	
	
	
	float get_ll_progressive(\
		std::vector<dirc_point> &hit_points,\
		int iBAR,\
		float iE,\
		float ix,\
		float iy,\
		float itheta,\
		float iphi,\
		float itracking_unc,\
		float ickov_unc);
	float get_ll_max(\
		std::vector<dirc_point> &hit_points,\
		int iBAR,\
		float iE,\
		float ix,\
		float iy,\
		float itheta,\
		float iphi,\
		float itracking_unc,\
		float ickov_unc);
};
#endif