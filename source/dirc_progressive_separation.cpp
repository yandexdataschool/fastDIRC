#include "../include/dirc_progressive_separation.h"
#include "../include/dirc_spread_gaussian.h"

DircProgressiveSeparation::DircProgressiveSeparation(\
	DircOpticalSim* imodel,\
	int imax_phots,\
	int istep_phots,\
	float ifudge_sigma, \
	float x_unc,\
	float y_unc,\
	float t_unc,\
	float im1 /*= .4937*/,\
	float im2 /*= .1396*/,\
	float ithresh /*= 20*/)
{
	max_sim_phots = imax_phots;
	step_sim_phots = istep_phots;
  
	mass_1 = im1;
	mass_2 = im2;
	ll_threshold = ithresh;
	
	//be good - try not to abuse side effects here...
	dirc_model = imodel;
	
	std::vector<dirc_point> no_hits;
	
	spread_func = new DircSpreadGaussian(\
		1,\
		no_hits,\
		x_unc,\
		y_unc,\
		t_unc);
	
	//Use rule of thumb times a ocnstant
	fudge_sigma = ifudge_sigma;
}
	
void DircProgressiveSeparation::set_max_step_phots(int im, int is)
{
	max_sim_phots = im;
	step_sim_phots = is;
}

void DircProgressiveSeparation::set_masses(float im1,float im2)
{
	mass_1 = im1;
	mass_2 = im2;
}
void DircProgressiveSeparation::set_threshold(float ithresh)
{
	ll_threshold = ithresh;
}
float DircProgressiveSeparation::ll_pts(\
	std::vector<dirc_point> &hit_points, \
	int num_support,\
	float beta)
{
	std::vector<dirc_point> support; 
	dirc_model->sim_rand_n_photons(\
		support,\
		num_support,\
		0,\
		BAR,\
		x,\
		y,\
		theta,\
		phi,\
		tracking_unc,\
		ckov_unc,\
		beta);
	
	return spread_func->get_log_likelihood_new_support(hit_points, support);
}
float DircProgressiveSeparation::ll_diff(\
	std::vector<dirc_point> &hit_points, \
	int num_support,\
	float beta_1,\
	float beta_2)
{
	float ll_1,ll_2;
	//NOTE: arg 2 is ignored for beta >= 0
	//Also, declare these vectors beforehand preferably - manage memory better for speed if needed
	std::vector<dirc_point> support_1; 
	dirc_model->sim_rand_n_photons(\
		support_1,\
		num_support,\
		0,\
		BAR,\
		x,\
		y,\
		theta,\
		phi,\
		tracking_unc,\
		ckov_unc,\
		beta_1);
	
	std::vector<dirc_point> support_2;
	dirc_model->sim_rand_n_photons(\
		support_2,\
		num_support,\
		0,\
		BAR,\
		x,\
		y,\
		theta,\
		phi,\
		tracking_unc,\
		ckov_unc,\
		beta_2);
	
	ll_1 = spread_func->get_log_likelihood_new_support(hit_points, support_1);
	ll_2 = spread_func->get_log_likelihood_new_support(hit_points, support_2);
	
	return ll_2 - ll_1;
}

float DircProgressiveSeparation::get_ll_progressive(\
	std::vector<dirc_point> &hit_points,\
	int iBAR,\
	float iE,\
	float ix,\
	float iy,\
	float itheta,\
	float iphi,\
	float itracking_unc,\
	float ickov_unc)
{
	float beta_1 = dirc_model->get_beta(iE,mass_1);
	float beta_2 = dirc_model->get_beta(iE,mass_2);
// 	printf("betas: %12.04f %12.04f\n",beta_1,beta_2);
	BAR = iBAR;
	x = ix;
	y = iy;
	theta = itheta;
	phi = iphi;
	tracking_unc = itracking_unc;
	ckov_unc = ickov_unc;
	std::vector<dirc_point> support_1; 
	std::vector<dirc_point> support_2; 
	std::vector<float> accum_likelihood_1;
	std::vector<float> accum_likelihood_2;
	std::vector<float> fill_likelihood_1;
	std::vector<float> fill_likelihood_2;
	int support_total_1 = 0;
	int support_total_2 = 0;
	
	for (unsigned int i = 0; i < hit_points.size(); i++)
	{
		//initialize arrays
		accum_likelihood_1.push_back(0);
		accum_likelihood_2.push_back(0);
	}
	
	float cur_ll = 0;
	float cur_ll_1 = 0;
	float cur_ll_2 = 0;
//	float add_ll_1 = 0;
//	float add_ll_2 = 0;
	
	for (int i = 0; i < max_sim_phots/step_sim_phots; i++)
	{
		//sim another set
		dirc_model->sim_rand_n_photons(\
			support_1,\
			step_sim_phots,\
			0,\
			BAR,\
			x,\
			y,\
			theta,\
			phi,\
			tracking_unc,\
			ckov_unc,\
			beta_1);
		
		support_total_1 += support_1.size();
		
		dirc_model->sim_rand_n_photons(\
			support_2,\
			step_sim_phots,\
			0,\
			BAR,\
			x,\
			y,\
			theta,\
			phi,\
			tracking_unc,\
			ckov_unc,\
			beta_2);
		
		support_total_2 += support_2.size();
		
		//3d rule of thumb - 3.46 ~ sqrt[12] for the uniform uncertainty of the grid
		//TODO fix how this applies to timing
// 		spread_func->set_gaus_sigma(std::max(10*fudge_sigma,fudge_sigma*pow(((float)step_sim_phots*i),-1./7.)));
		spread_func->set_gaus_sigma(fudge_sigma);
		
		spread_func->fill_likelihood_new_support(\
			fill_likelihood_1,\
			support_1,\
			hit_points);
		
		spread_func->fill_likelihood_new_support(\
			fill_likelihood_2,\
			support_2,\
			hit_points);
		
		cur_ll_1 = 0;
		cur_ll_2 = 0;
		
		bool all_hit_points_covered = true;
		for(unsigned int j = 0; j < fill_likelihood_1.size(); j++)
		{
			accum_likelihood_1[j] += fill_likelihood_1[j];
			accum_likelihood_2[j] += fill_likelihood_2[j];
			
			all_hit_points_covered = all_hit_points_covered && (accum_likelihood_1[j] > 1e-6) && (accum_likelihood_2[j] > 1e-6);
			
			cur_ll_1 += log(accum_likelihood_1[j]/support_total_1);
			cur_ll_2 += log(accum_likelihood_2[j]/support_total_2);
			
// 			printf("accum_likelihood_1 accum_likelihood_2: %04d %12.04f %12.04f\n",j,accum_likelihood_1[j],accum_likelihood_2[j]);
// 			printf("fill_likelihood_1 fill_likelihood_2: %04d %12.04f %12.04f\n",j,fill_likelihood_1[j],fill_likelihood_2[j]);
		  
		}
		
		//One or more of the points does not have hits around it - throw more
// 		if (all_hit_points_covered == false) continue;
		
		cur_ll = cur_ll_2 - cur_ll_1;
		
		if (fabs(cur_ll) > ll_threshold)
		{
			//Assume we've IDed it and return;
			return cur_ll;
		}
	}
	return cur_ll;
}
void DircProgressiveSeparation::get_inf_and_ll(\
	int &num_neg_inf,\
	float &ll_sum,\
	std::vector<float> prob_vals)
{
	num_neg_inf = 0;
	ll_sum = 0;
	for(unsigned int i = 0; i < prob_vals.size(); i++)
	{
		if (prob_vals[i] < 1e-5)
		{
			//too small
			num_neg_inf++;
		}
		else
		{
			ll_sum += log(prob_vals[i]);
		}
	}
}
float DircProgressiveSeparation::get_ll_max(\
	std::vector<dirc_point> &hit_points,\
	int iBAR,\
	float iE,\
	float ix,\
	float iy,\
	float itheta,\
	float iphi,\
	float itracking_unc,\
	float ickov_unc)
{
	float beta_1 = dirc_model->get_beta(iE,mass_1);
	float beta_2 = dirc_model->get_beta(iE,mass_2);
	
	x = ix;
	y = iy;
	theta = itheta;
	phi = iphi;
	tracking_unc = itracking_unc;
	ckov_unc = ickov_unc;
	
	return ll_diff(\
		hit_points,\
		max_sim_phots,\
		beta_1,\
		beta_2);
		
}
