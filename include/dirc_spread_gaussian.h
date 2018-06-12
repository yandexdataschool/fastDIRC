#include "dirc_point.h"
#include "dirc_spread_radius.h"
#include <vector>
//#include <math.h>
#include <TRandom3.h>
#ifndef DIRC_SPREAD_GAUSSIAN
#define DIRC_SPREAD_GAUSSIAN
// class DircSpreadGaussian : public DircSpreadRadius
// class DircSpreadGaussian : public DircProbabilitySpread
class DircSpreadGaussian
{
private:
	double x_sig2inv,y_sig2inv,t_sig2inv;
	double spread_func_norm, spread_func_norm_inv;
	double lin_slope, r_trans, sigma2, sigma2inv,max_val;
	TRandom3 rand_gen;
	std::vector<dirc_point> support_points;
	double get_weight(dirc_point inpoint);
	double min_probability;
public:
	//by reference - later for speed
	DircSpreadGaussian(\
		double isigma, \
		std::vector<dirc_point> isupport,\
		double x_unc,\
		double y_unc,\
		double t_unc,\
		double imin_prob = 1e-6);
	void support_spread(double spread_sig);
	void support_x_weight();
	void set_support(std::vector<dirc_point> isupport);
	void set_gaus_sigma(double isigma);

	const inline double radius_spread_function(const double r2) {
	    if (r2 < 5*sigma2) {
		return exp(-r2*sigma2inv);
	    } else {
		return 0.;
	    }
	};

	const inline double support_spread_function(const dirc_point& support,
						    const dirc_point& test) {
		double dx2,dy2,dt2;
		dx2 = support.x - test.x;
		dx2 *= dx2;
		dy2 = support.y - test.y;
		dy2 *= dy2;
		dt2 = support.t - test.t;
		dt2 *= dt2;
		return radius_spread_function(dx2*x_sig2inv+dy2*y_sig2inv+dt2*t_sig2inv);
	};

	double get_single_likelihood(dirc_point inpoint);
	const double get_log_likelihood(const std::vector<dirc_point>& inpoints);
	double get_log_likelihood_new_support(
                std::vector<dirc_point> &inpoints, std::vector<dirc_point> &t_support);
	void fill_likelihood_new_support(\
		std::vector<double> &likelihood_vals,\
		std::vector<dirc_point> new_support,\
		std::vector<dirc_point> inpoints);
};
#endif
