#ifndef DIRC_SPREAD_GAUSSIAN
#define DIRC_SPREAD_GAUSSIAN
#include <vector>
#include <TRandom3.h>
#include "dirc_point.h"
#include "dirc_spread_radius.h"

class DircSpreadGaussian {
private:
	float x_sig2inv,y_sig2inv,t_sig2inv;
	float spread_func_norm, spread_func_norm_inv;
	float lin_slope, r_trans, sigma2, sigma2inv,max_val;
	std::unique_ptr<TRandom3> rand_gen;
	std::vector<dirc_point> support_points;
	float get_weight(dirc_point inpoint);
	float min_probability;
public:
	//by reference - later for speed
	DircSpreadGaussian(\
		float isigma, \
		std::vector<dirc_point> isupport,\
		float x_unc,\
		float y_unc,\
		float t_unc,\
		float imin_prob = 1e-6);
	void support_spread(float spread_sig);
	void support_x_weight();
	void set_support(std::vector<dirc_point> isupport);
	void set_gaus_sigma(float isigma);

	const inline float radius_spread_function(const float r2) {
	    if (r2 < 5*sigma2) {
		return exp(-r2*sigma2inv);
	    } else {
		return 0.;
	    }
	};

	const inline float support_spread_function(const dirc_point& support,
						    const dirc_point& test) {
		float dx2,dy2,dt2;
		dx2 = support.x - test.x;
		dx2 *= dx2;
		dy2 = support.y - test.y;
		dy2 *= dy2;
		dt2 = support.t - test.t;
		dt2 *= dt2;
		return radius_spread_function(dx2*x_sig2inv+dy2*y_sig2inv+dt2*t_sig2inv);
	};

	float get_single_likelihood(dirc_point inpoint);
	const float get_log_likelihood(const std::vector<dirc_point>& inpoints);
	float get_log_likelihood_new_support(
                std::vector<dirc_point> &inpoints, std::vector<dirc_point> &t_support);
	void fill_likelihood_new_support(\
		std::vector<float> &likelihood_vals,\
		std::vector<dirc_point> new_support,\
		std::vector<dirc_point> inpoints);
};
#endif
