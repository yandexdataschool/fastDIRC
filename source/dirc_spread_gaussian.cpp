#include "../include/dirc_point.h"
#include "../include/dirc_spread_gaussian.h"
#include <vector>
#include <math.h>

DircSpreadGaussian::DircSpreadGaussian(
	float isigma,
	std::vector<dirc_point> isupport,
	float x_unc,
	float y_unc,
	float t_unc,
	float imin_prob) {
	sigma2 = isigma*isigma;
	sigma2inv = 1/sigma2;
// 	spread_func_norm=5.56833*isigma;
	spread_func_norm = 1;
	spread_func_norm_inv=1/spread_func_norm;
	
	
	x_sig2inv = 1/(x_unc*x_unc);
	y_sig2inv = 1/(y_unc*y_unc);
	t_sig2inv = 1/(t_unc*t_unc);
	
	min_probability = imin_prob;

	support_points = isupport;
	
}
void DircSpreadGaussian::set_support(std::vector<dirc_point> isupport)
{
	support_points = isupport;
}
void DircSpreadGaussian::support_spread(float spread_sig)
{
	unsigned int start_size = support_points.size();
	int mult_add = 1;
	for (unsigned int i = 0; i < start_size; i++)
	{
		for (int j = 0; j < mult_add; j++)
		{
			dirc_point new_point;
			new_point.x = support_points[i].x + rand_gen->Gaus(0,spread_sig);
			new_point.y = support_points[i].y + rand_gen->Gaus(0,spread_sig);
			new_point.t = support_points[i].t;
			support_points.push_back(new_point);
		}
	}
}
void DircSpreadGaussian::support_x_weight()
{
	float x;
	for (unsigned int i = 0; i < support_points.size(); i++)
	{
		x = support_points[i].x;
		support_points[i].weight = 1 + sqrt(fabs(x));
	}
}

float DircSpreadGaussian::get_single_likelihood(dirc_point inpoint)
{
	float tprob = 0;
//	float log_mult = 1;
//	float weight = 1;
// 	weight = get_weight(inpoint);
	for (unsigned int j = 0; j < support_points.size(); j++)
	{
		tprob += support_spread_function(support_points[j],inpoint);
	}
	tprob /= support_points.size();
	
// 	tprob = std::max(tprob,min_probability);
	
	return tprob+min_probability;
}

const float DircSpreadGaussian::get_log_likelihood(const std::vector<dirc_point>& inpoints) {
        float tprob;
	float rval = 0;
	float log_mult = 1.;
	float weight = 1.;
	for (unsigned int i = 0; i < inpoints.size(); ++i)
	{
		tprob = 0;
		for (unsigned int j = 0; j < support_points.size(); ++j)
		{
			tprob += support_spread_function(support_points[j], inpoints[i]);
		}
		tprob /= support_points.size();
		tprob *= spread_func_norm_inv;
		
		// TODO deal with normalization....
		rval += weight*log_mult*log(tprob+min_probability);
	}
	rval -= log(inpoints.size());
	return rval;
}

void DircSpreadGaussian::fill_likelihood_new_support(\
	std::vector<float> &likelihood_vals,\
	std::vector<dirc_point> new_support,\
	std::vector<dirc_point> inpoints)
{
	// unnormailized to support point size
	likelihood_vals.clear();
	float tprob = 0;
	int eval_count = 0;
	for (unsigned int i = 0; i < inpoints.size(); ++i)
	{
		tprob = 0;
		for (unsigned int j = 0; j < new_support.size(); ++j)
		{
			tprob += support_spread_function(new_support[j],inpoints[i]);
			eval_count++;
		}
		tprob *= spread_func_norm_inv;
		
		likelihood_vals.push_back(tprob);
	}
}
float DircSpreadGaussian::get_log_likelihood_new_support(std::vector<dirc_point> &inpoints, std::vector<dirc_point> &t_support)
{
  //Ideally, I find a way tp do this without replicating code....
  //maybe change above func to call this
	float tprob = 0;
	float rval = 0;
	int eval_count = 0;
	float log_mult = 1;
	float weight = 1;
	
	for (unsigned int i = 0; i < inpoints.size(); i++)
	{
		tprob = 0;
		for (unsigned int j = 0; j < t_support.size(); j++)
		{
			tprob += support_spread_function(t_support[j],inpoints[i]);
			eval_count++;
		}
		tprob /= t_support.size();
		tprob *= spread_func_norm_inv;
		
		rval += weight*log_mult*log(tprob+min_probability);
		
	}
// 	rval -= log(inpoints.size());
	
	return rval;
}
float DircSpreadGaussian::get_weight(dirc_point inpoint)
{
	return 1;
}
void DircSpreadGaussian::set_gaus_sigma(float isigma)
{
	sigma2 = isigma*isigma;
	sigma2inv = 1/sigma2;
// 	spread_func_norm=5.56833*isigma;
	spread_func_norm = 1;
	spread_func_norm_inv=1/spread_func_norm;
}
