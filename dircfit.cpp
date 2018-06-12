#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <utility>

#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "include/dirc_optical_sim.h"
#include "include/dirc_threesegbox_sim.h"
#include "include/dirc_point.h"
#include "include/dirc_probability_spread.h"
#include "include/dirc_probability_separation.h"
#include "include/dirc_spread_radius.h"
#include "include/dirc_spread_relative.h"
#include "include/dirc_spread_linear_soft.h"
#include "include/dirc_spread_gaussian.h"
#include "include/dirc_digitizer.h"
#include "include/dirc_rect_digitizer.h"
#include "include/dirc_babar_digitizer.h"
#include "include/dirc_babar_sim.h"
#include "include/dirc_progressive_separation.h"
#include "include/dirc_gluex_lut_enum.h"
#include "include/dirc_lut_enum.h"
#include "include/dirc_lut.h"
#include <TFile.h>
#include <TTree.h>
#include <TH3.h>
#include <TH2.h>
#include <TH1.h>
#include <TF1.h>
#include <TRandom3.h>
#include <TMinuit.h>


//#define USE_BABAR_BOX

DircBaseSim* global_calib_model;
DircDigitizer* global_digitizer;
float global_particle_theta;
float global_particle_phi;
float global_particle_x;
float global_particle_y;
float global_particle_t;
float global_particle_energy;
float global_particle_flight_distance;
float global_refrac_index;
float global_ckov_unc;
float global_tracking_unc;
int global_n_sim_phots;

std::vector<dirc_point> global_pion_line;
std::vector<dirc_point> global_kaon_line;
std::vector<std::vector<dirc_point> > global_pion_hits;
std::vector<std::vector<dirc_point> > global_kaon_hits;

void init_global_data()
{	
	int num_line_points = 5000;
	int num_particle_hits = 20000;

	float kmass = .4937;
	float pimass = .1396;

	global_pion_line.clear();
	global_kaon_line.clear();
	global_pion_hits.clear();
	global_kaon_hits.clear();

	DircBaseSim* dirc_model = global_calib_model;
	DircDigitizer* digitizer = global_digitizer;

	float particle_theta = global_particle_theta;
	float particle_phi = global_particle_phi;
	float particle_x = global_particle_x;
	float particle_y = global_particle_y;
	float energy = global_particle_energy;
	float particle_flight_distance = global_particle_flight_distance;
//	float particle_t = global_particle_t;
	int n_sim_phots = global_n_sim_phots;

	float tracking_unc = global_tracking_unc;
	float ckov_unc = global_ckov_unc;

	float pion_old_angle = -1;
	float kaon_old_angle = -1;

//`	printf("%12.04f %12.04f %12.04f %12.04f %12.04f %12.04f %12.04f\n",particle_theta,particle_phi,particle_x,particle_y,energy,particle_flight_distance,ckov_unc);

	float refrac_index = global_refrac_index;
//	printf("%12.04f\n",refrac_index);
	
	float refrac_split = -.000;
	float pion_refrac_mod = 1.000;
	float kaon_refrac_mod = 1.000;

	pion_refrac_mod += refrac_split;
	kaon_refrac_mod -= refrac_split;

	float pion_beta = dirc_model->get_beta(energy,pimass);
	float kaon_beta = dirc_model->get_beta(energy,kmass);
	//ns
	float pion_time = particle_flight_distance/(pion_beta*.3);
	float kaon_time = particle_flight_distance/(kaon_beta*.3);

	float pion_emit_angle = 57.3*acos(1/(pion_beta*pion_refrac_mod*refrac_index));
	float kaon_emit_angle = 57.3*acos(1/(kaon_beta*kaon_refrac_mod*refrac_index));

//	printf("angles: %12.04f %12.04f   Beta: %12.04f %12.04f\n",pion_emit_angle,kaon_emit_angle,pion_beta,kaon_beta);

	std::vector<dirc_point> provisional_points_lwn;
	std::vector<dirc_point> provisional_points_lwp;
	dirc_model->track_all_line_photons(\
			provisional_points_lwn,\
			provisional_points_lwp,\
			num_line_points,\
			pion_emit_angle,\
			particle_theta,\
			particle_phi,\
			particle_x,\
			particle_y,\
			-17.25/2,\
			pion_time,\
			1);
	for (unsigned int i = 0; i < provisional_points_lwn.size(); i++)
	{
		global_pion_line.push_back(provisional_points_lwn[i]);
	}
	for (unsigned int i = 0; i < provisional_points_lwp.size(); i++)
	{
		global_pion_line.push_back(provisional_points_lwp[i]);
	}
	provisional_points_lwn.clear();
	provisional_points_lwp.clear();
	dirc_model->track_all_line_photons(\
			provisional_points_lwn,\
			provisional_points_lwp,\
			num_line_points,\
			kaon_emit_angle,\
			particle_theta,\
			particle_phi,\
			particle_x,\
			particle_y,\
			-17.25/2,\
			kaon_time,\
			1);
	for (unsigned int i = 0; i < provisional_points_lwn.size(); i++)
	{
		global_kaon_line.push_back(provisional_points_lwn[i]);
	}
	for (unsigned int i = 0; i < provisional_points_lwp.size(); i++)
	{
		global_kaon_line.push_back(provisional_points_lwp[i]);
	}

	//Lines are drawn, now simulate data
	for (int i = 0; i < num_particle_hits; i++)
	{
		std::vector<dirc_point> pion_push;
		std::vector<dirc_point> kaon_push;
		dirc_model->sim_rand_n_photons(\
				pion_push,\
				n_sim_phots,\
				pion_old_angle,\
				1,\
				particle_x,\
				particle_y,\
				pion_time,\
				particle_theta,\
				particle_phi,\
				tracking_unc,\
				ckov_unc,\
				pion_beta);
		digitizer->digitize_points(pion_push);
	
		dirc_model->sim_rand_n_photons(\
				kaon_push,\
				n_sim_phots,\
				kaon_old_angle,\
				1,\
				particle_x,\
				particle_y,\
				kaon_time,\
				particle_theta,\
				particle_phi,\
				tracking_unc,\
				ckov_unc,\
				kaon_beta);
		digitizer->digitize_points(kaon_push);

		//printf("%5u %5u\n",pion_push.size(),kaon_push.size());

		//printf("model run: %12.04f %12.04f %12.04f %12.04f %12.04f %12.04f %12.04f %12.04f %12.04f\n",kaon_emit_angle,particle_x,particle_y,kaon_time,particle_theta,particle_phi,tracking_unc,ckov_unc,kaon_beta);
		global_pion_hits.push_back(pion_push);
		global_kaon_hits.push_back(kaon_push);
	}
}


void dirc_calib_line(int &npar, float *gin, float &f, float *par, int iflag)
{
	//do not compute gradient, so flag is irrelevant
	//pars: pion_x_adj, pion_y_adj, kaon_x_adj, kaon_y_adj

	TH1F *ll_diff_pion = new TH1F("ll_diff_pion_tmp","Difference of log likelihood real = pion",60000,-600,600);
	TH1F *ll_diff_kaon = new TH1F("ll_diff_kaon_tmp","Difference of log likelihood real = kaon",60000,-600,600);

	float time_spread = 1000;
	//float dist_spread = 40;
	float dist_spread = 50;
	float max_dev_sq = 5;

	max_dev_sq *= max_dev_sq;
	
	float llc = 0;
	float llf = 0;

	float mult_params = 1;

	float x_shift = par[0];
	x_shift = -.6;
	float y_shift = par[0];
	float y_split = par[1];

	float pion_x_adj = x_shift;
	float kaon_x_adj = x_shift;
	float pion_y_adj = mult_params*(y_shift+y_split);
	float kaon_y_adj = mult_params*(y_shift-y_split);


	std::vector<dirc_point> pion_line;
	std::vector<dirc_point> kaon_line;

	for (unsigned int i = 0; i < global_pion_line.size(); i++)
	{
		dirc_point add_point = global_pion_line[i];
		add_point.x += pion_x_adj;
		add_point.y += pion_y_adj;
		pion_line.push_back(add_point);

	//	printf("%12.04f %12.04f %12.04f\n",global_pion_line[i].y,global_kaon_line[i].y,global_pion_line[i].y-global_kaon_line[i].y);
	}
	for (unsigned int i = 0; i < global_kaon_line.size(); i++)
	{
		dirc_point add_point = global_kaon_line[i];
		add_point.x += kaon_x_adj;
		add_point.y += kaon_y_adj;
		kaon_line.push_back(add_point);
	}


	std::vector<dirc_point> sim_points;
	for (unsigned int k = 0; k < global_pion_hits.size(); k++)
	{//global pion hits and global kaon hits should have the same size
		//do pion hit first
		sim_points = global_pion_hits[k];
		
		llc = 0;
		llf = 0;

		
		for (unsigned int i = 0; i < sim_points.size(); i++)
		{
			float min_dist_sq = 10000;
			float cur_dist_sq = -1;
			for (unsigned int j = 0; j < pion_line.size(); j++)
			{
				cur_dist_sq = (sim_points[i].x - pion_line[j].x)*(sim_points[i].x - pion_line[j].x);		
				cur_dist_sq += (sim_points[i].y - pion_line[j].y)*(sim_points[i].y - pion_line[j].y);		
				cur_dist_sq /= dist_spread;
				cur_dist_sq += (sim_points[i].t - pion_line[j].t)*(sim_points[i].t - pion_line[j].t)/time_spread;	
				min_dist_sq = std::min(min_dist_sq,cur_dist_sq);
				//printf("%12.04f\n",cur_dist_sq);
			}

			llc += std::min(min_dist_sq,max_dev_sq);
			//printf("%4d %12.04f %12.04f\n",i,min_dist_sq,llc);
		}
		for (unsigned int i = 0; i < sim_points.size(); i++)
		{
			float min_dist_sq = 10000;
			float cur_dist_sq = -1;
			for (unsigned int j = 0; j < kaon_line.size(); j++)
			{
				cur_dist_sq = (sim_points[i].x - kaon_line[j].x)*(sim_points[i].x - kaon_line[j].x);		
				cur_dist_sq += (sim_points[i].y - kaon_line[j].y)*(sim_points[i].y - kaon_line[j].y);		
				cur_dist_sq /= dist_spread;
				cur_dist_sq += (sim_points[i].t - kaon_line[j].t)*(sim_points[i].t - kaon_line[j].t)/time_spread;	
				min_dist_sq = std::min(min_dist_sq,cur_dist_sq);
			}
			//cur_mean_phi /= total_prov_points;
			//printf("pion/pion: %12.04f %12.04f\n",cur_mean_phi/total_prov_points,atan2(cur_mean_y,cur_mean_x));
			llf += std::min(min_dist_sq,max_dev_sq);
		}
		ll_diff_pion->Fill(1*(llc-llf));

		llc=0;
		llf=0;

		sim_points = global_kaon_hits[k];
		for (unsigned int i = 0; i < sim_points.size(); i++)
		{
			float min_dist_sq = 10000;
			float cur_dist_sq = -1;
			for (unsigned int j = 0; j < pion_line.size(); j++)
			{
				cur_dist_sq = (sim_points[i].x - pion_line[j].x)*(sim_points[i].x - pion_line[j].x);		
				cur_dist_sq += (sim_points[i].y - pion_line[j].y)*(sim_points[i].y - pion_line[j].y);		
				cur_dist_sq /= dist_spread;
				cur_dist_sq += (sim_points[i].t - pion_line[j].t)*(sim_points[i].t - pion_line[j].t)/time_spread;	
				min_dist_sq = std::min(min_dist_sq,cur_dist_sq);
			}
			llc += std::min(min_dist_sq,max_dev_sq);
		}

		for (unsigned int i = 0; i < sim_points.size(); i++)
		{
			float min_dist_sq = 10000;
			float cur_dist_sq = -1;
			for (unsigned int j = 0; j < kaon_line.size(); j++)
			{
				cur_dist_sq = (sim_points[i].x - kaon_line[j].x)*(sim_points[i].x - kaon_line[j].x);		
				cur_dist_sq += (sim_points[i].y - kaon_line[j].y)*(sim_points[i].y - kaon_line[j].y);		
				cur_dist_sq /= dist_spread;
				cur_dist_sq += (sim_points[i].t - kaon_line[j].t)*(sim_points[i].t - kaon_line[j].t)/time_spread;	
				min_dist_sq = std::min(min_dist_sq,cur_dist_sq);
			}
			llf += std::min(min_dist_sq,max_dev_sq);

		}
		ll_diff_kaon->Fill(1*(llc-llf));
	}
	//I feel unclean
	//ll_diff_pion->Write();
	//ll_diff_kaon->Write();

	//ll_diff_pion->Rebin(20);
	//ll_diff_kaon->Rebin(20);

        std::vector<float> xr;
        std::vector<float> yr;
        float ival = 0;
	float pion_integral = 0;
	float kaon_reverse_integral = 1;
	float integral_scale = global_pion_hits.size();

        for (int i = 0; i < ll_diff_pion->GetNbinsX(); i++)
        {
		pion_integral += ll_diff_pion->GetBinContent(i)/integral_scale;
		kaon_reverse_integral -= ll_diff_kaon->GetBinContent(i)/integral_scale;
                xr.push_back(pion_integral);
                yr.push_back(kaon_reverse_integral);
		if (xr[i] < .001 && xr[i] > .00)
		{
			//printf("%8d %12.04f %12.04f\n",i,xr[i],yr[i]);
       		}
		if (ll_diff_pion->GetBinContent(i)/integral_scale > 0)
		{
			//printf("%8d %12.04f\n",i,ll_diff_pion->GetBinContent(i)/integral_scale);
		}
	 }

        float last_x = xr[0];
        float last_y = yr[0];

        for (int i = 0; i < ll_diff_pion->GetNbinsX()-1; i++)
        {
                ival += (yr[i]+last_y)*(xr[i] - last_x)/2;
		//printf("%6d %12.09f %12.04f %12.04f %12.04f %12.04f\n",i,ival,xr[i],yr[i],last_x,last_y);
                last_x = xr[i];
                last_y = yr[i];
        }
	
//	printf("params: %12.04f %12.04f %12.04f %12.04f %12.04f  val: %12.04f\n",pion_x_adj,pion_y_adj,kaon_x_adj,kaon_y_adj,dist_spread,-ival);
	f = -ival; //return negative ROC integral
	delete ll_diff_pion;
	delete ll_diff_kaon;
}


std::vector<dirc_point> fold_x(std::vector<dirc_point> inpoints) {
	std::vector<dirc_point> outvec;
	for (unsigned int i = 0; i < inpoints.size(); i++) {
		dirc_point tpoint = inpoints[i];
		tpoint.x = fabs(tpoint.x);
		outvec.push_back(tpoint);
	}
	return outvec;
}
int main(int nargs, char* argv[])
{  

//	printf("%12.04f\n",atan2(0,10)*57.3);

	clock_t timing_clock;

	const char* in_str;
	bool inputfile = false;
	bool inputrootfile = false;
	bool out_csv = false;
	bool slac_run = false;
	int output_box_angles_n = -1;
	float time_window=-1;//time window for confounded pmt hits, in ns	

	float energy = 5.0;
	float energy_mean = energy;
	float energy_spread = 0;
	float kmass = .4937;
	float pimass = .1396;
	float mumass = .1057;

	float particle_x = 0;
	float particle_y = 0;
	float particle_x_mean = particle_x;
	float particle_y_mean = particle_y;
	float particle_x_spread = 0;
	float particle_y_spread = 0;
	float particle_theta = 4;
	float particle_theta_mean = particle_theta;
	float particle_theta_spread = 0;
	float particle_phi = 40;
	float particle_phi_mean = particle_phi;
	float const_track_off = 0;

	float particle_flight_distance = 0;

	int box_check_n = -1;
	float box_check_theta = 0;
	float box_check_phi = 0;
	float box_check_theta_unc = 0;
	float box_check_phi_unc = 0;
	float box_check_overall_theta = -15;
	int phot_check_n = -1;
	float phot_check_max_theta = 4;


	bool force_kinematics = false;
	bool use_prog_sep = false;
	bool kaleidoscope_plot = false;
	bool monochrome_plot = false;
	bool flatten_time = false;
	bool sep_updown = false;
	bool output_peak_lambda = false;
	bool lut_slac = false;
	bool run_minuit_calibration = false;
	bool perform_chromatic_correction = true;
	bool use_moliere_scattering = false;

	int num_runs = 1000;
	int max_particles = 60000000;
	int phi_phots_reduce = 1;
	int refraction_sim_n = -1;
	int sparse_sim_n = -1;
	int sparse_recon_n = -1;
	int line_recon_n = -1;
	int line_output_n = -1;
	int fill_d_midline_n = -1;
	int lut_sim_n = -1;
	int gaus_ll_n = -1;
	int sim_time_test_n = -1;
	int bounces_phot_n = -1;

	float gaus_ll_spread = 2;

	float mean_n_phot = 40;
	float spread_n_phot = 0;

	float wedge_uncertainty = 0/57.3;
	float refrac_index=1.47;
	float mirror_angle_change = 0;
	float mirror_angle_change_unc = 0;
	float mirror_angle_change_yunc = 0;
	float box_rot = 0;
	float box_rot_unc = 0;
	float bar_box_box_angle = 0/57.3;
	float mirror_r_difference = 400;//1200 - 400 = 800.  Changed 05/09/2016.  Does not affect threeseg mirror reconstruction as far as I can tell - this was known.
	//	float mirror_r_difference = 0;
	float wedge_non_uniformity = 0;
	float pmt_offset = 0;
	float main_mirror_nonuniformity = 0;
	float foc_mirror_size = 288;


	float main_mirror_angle_off = 0;
	float main_mirror_yangle_off = 0;
	float main_mirror_zangle_off = 0;
	float main_mirror_yoff = 0;
	float main_mirror_zoff = 0;

	float bar_box_xoff = 0;
	float bar_box_yoff = 0;
	float bar_box_zoff = 0;

	float pmt_min_z = -1000;
	float pmt_max_z = 1000;
	float large_mirror_min_z = -1000;
	float large_mirror_max_z = 1000;

	pmt_min_z = -559;
	pmt_max_z = -329;
	//pmt_max_z = -338;
	large_mirror_min_z = -559;
	large_mirror_max_z = -130;

	float upper_wedge_yang_spread = 0;
	int rseed = 1337;

	int broaden_events = 0;

	float tracking_unc = .0000*57.3; //mrad
	// 	float ckov_unc = .0077*57.3; //chromatic + optical aberation = 7.7mrad
	float ckov_unc = .003*57.3; //transport = 3mrad

	float resx = 6;
	float resy = 6;
	float rest = 1;
	//	float minx = -8000;
	//	float maxx = -minx;
	//	float miny = -800;
	//	float maxy = -miny;
	float minx = -1500;
	float maxx = 1500;
	float miny = -500;
	float maxy = 500;
	float mint = 0;
	float maxt = 1000;
	float t_unc = .27;
	float t_bin_size = 1;

	float digit_miny = -50;
	float digit_maxy = 300;

	digit_miny = miny;
	digit_maxy = maxy;



	//Sets the side boundarys of the distributions
	float sm_xl = -10000000;
	float sm_xr = -sm_xl;

	//float s_func_x = 6;
	float s_func_x = 6;
	float s_func_y = s_func_x;
	//float s_func_t = 2;
	float s_func_t = 1.0;
	float sfunc_sig = 1;

#ifdef USE_BABAR_BOX
	s_func_x = 29;
	s_func_y = 29;
#endif 

	int n_sim_phots = 40;

//	float pion_x_adj = -.43;
//	float pion_y_adj = -3.6;
//	float kaon_x_adj = -.48;
//	float kaon_y_adj = 2.87;

	float pion_x_adj = 0;
	float pion_y_adj = 0;
	float kaon_x_adj = 0;
	float kaon_y_adj = 0;

	//int n_phi_phots = 900000;
	int n_phi_phots = 150000;
	//n_phi_phots = 75000;
	int n_z_phots = 4;
	// 	n_step_phots = n_z_phots*n_phi_phots;

/*
	sm_xl = -50;
	sm_xr = sm_xl + 440;
*/
//	sm_xl = 0;
//	sm_xr = sm_xl + 1100;

//	float overlap_x = -1;

	bool use_quartz_for_liquid = false;
	bool three_seg_mirror = true;
	bool fill_distributions	= false;
	bool fill_kinematics_yields = false;
	float kinematics_yields_min_theta = 0;
	float kinematics_yields_max_theta = 12;
	float kinematics_yields_step_theta = .25;
	float kinematics_yields_min_phi = 0;
	float kinematics_yields_max_phi = 45;
	float kinematics_yields_step_phi = 1;
	int  kinematics_n_phots = 100000;


	float liquid_absorbtion = 0*-log(.7)/1000;
	float liquid_index = 1.33;

	bool coverage_plot = false;
	int num_cov = 100000;

	char* rootfilename = new char[256];
	char* inputrootfilename = new char[256];
	sprintf(rootfilename,"fitdirc.root");	

	printf("Arguments Passed=%d\n",nargs);

	if(nargs==2){
		in_str = argv[1];
		printf("%s\n",in_str);
		printf("nargs=2\n");
		inputfile = true;
	}
	else{
		for (int i = 1; i < nargs; i++)
		{
			if (strcmp(argv[i], "-if") == 0)
			{
				i++;
				in_str = argv[i];
				inputfile = true;
			}
			else if (strcmp(argv[i], "-of") == 0)
			{
				i++;
				sprintf(rootfilename,"%s",argv[i]);	
			}
			else if (strcmp(argv[i], "-fill_dist") == 0)
			{
				fill_distributions = true;	
			}
			else if (strcmp(argv[i], "-cylindrical_mirror") == 0)
			{
				three_seg_mirror = false;	
			}
			else if (strcmp(argv[i], "-updown") == 0)
			{
				printf("Up-Down Histogram filling only implemented in loop mode, option currently does nothing in other modes\n");
				sep_updown = true;	
			}
			else if (strcmp(argv[i], "-root_input") == 0)
			{
				i++;
				inputrootfile = true;
				sprintf(inputrootfilename,"%s",argv[i]);
			}
			else if (strcmp(argv[i], "-particle_phi") == 0)
			{
				i++;
				particle_phi = atof(argv[i]);
				particle_phi_mean = particle_phi;
			}
			else if (strcmp(argv[i], "-particle_flight_distance") == 0)
			{
				//meters
				i++;
				particle_flight_distance = atof(argv[i]);
			}
			else if (strcmp(argv[i], "-tracking_unc") == 0)
			{
				i++;
				tracking_unc = atof(argv[i]);
			}
			else if (strcmp(argv[i], "-output_peak_lambda") == 0)
			{
				output_peak_lambda = true;
			}
			else if (strcmp(argv[i], "-const_track_off") == 0)
			{
				i++;
				const_track_off = atof(argv[i]);
			}
			else if (strcmp(argv[i], "-ckov_unc") == 0)
			{
				i++;
				ckov_unc = atof(argv[i]);
			}
			else if (strcmp(argv[i], "-refraction_sim_n") == 0)
			{
				i++;
				refraction_sim_n = atoi(argv[i]);
			}
			else if (strcmp(argv[i], "-sparse_sim_n") == 0)
			{
				i++;
				sparse_sim_n = atoi(argv[i]);
			}
			else if (strcmp(argv[i], "-line_output_n") == 0)
			{
				i++;
				line_output_n = atoi(argv[i]);
			}
			else if (strcmp(argv[i], "-sparse_recon_n") == 0)
			{
				i++;
				sparse_recon_n = atoi(argv[i]);
			}
			else if (strcmp(argv[i], "-line_recon_n") == 0)
			{
				i++;
				line_recon_n = atoi(argv[i]);
			}
			else if (strcmp(argv[i], "-side_mirror") == 0)
			{
				i++;
				sm_xl = atof(argv[i]);
				i++;
				sm_xr = atof(argv[i]);
			}
			else if (strcmp(argv[i], "-line_recon_calib") == 0)
			{
				//setting line recon calibrations, expects 4 arguements
				i++;
				pion_x_adj = atof(argv[i]);
				i++;
				pion_y_adj = atof(argv[i]);
				i++;
				kaon_x_adj = atof(argv[i]);
				i++;
				kaon_y_adj = atof(argv[i]);
			}
			else if (strcmp(argv[i], "-lut_sim_n") == 0)
			{
				i++;
				lut_sim_n = atoi(argv[i]);
			}
			else if (strcmp(argv[i], "-lut_slac") == 0)
			{
				lut_slac = true;
			}
			else if (strcmp(argv[i], "-fill_d_midline_n") == 0)
			{
				i++;
				fill_d_midline_n = atoi(argv[i]);
			}
			else if (strcmp(argv[i], "-gaus_ll_n") == 0)
			{
				i++;
				gaus_ll_n = atoi(argv[i]);
			}
			else if (strcmp(argv[i], "-gaus_ll_spread") == 0)
			{
				i++;
				gaus_ll_spread = atof(argv[i]);
			}
			else if (strcmp(argv[i], "-particle_theta") == 0)
			{
				i++;
				particle_theta = atof(argv[i]);
				particle_theta_mean = particle_theta;
			}
			else if (strcmp(argv[i], "-sim_time_test_n") == 0)
			{
				i++;
				sim_time_test_n = atoi(argv[i]);
			}
			else if (strcmp(argv[i], "-bounces_phot_n") == 0)
			{
				i++;
				bounces_phot_n = atoi(argv[i]);
			}
			else if (strcmp(argv[i], "-use_quartz_for_liquid") == 0)
			{
				use_quartz_for_liquid = true;
			}
			else if (strcmp(argv[i], "-use_moliere_scattering") == 0)
			{
				use_moliere_scattering = true;
				printf("Enabling Moliere Scattering - only implemented in loop mode\n");
			}
			else if (strcmp(argv[i], "-force_kinematics") == 0)
			{
				force_kinematics = true;
			}
			else if (strcmp(argv[i], "-fill_kinematics_yields") == 0)
			{
				fill_kinematics_yields = true;
			}
			else if (strcmp(argv[i], "-out_csv") == 0)
			{
				out_csv = true;
			}
			else if (strcmp(argv[i], "-coverage_plot") == 0)
			{
				coverage_plot = true;
			}
			else if (strcmp(argv[i], "-kaleidoscope_plot") == 0)
			{
				kaleidoscope_plot = true;
			}
			else if (strcmp(argv[i], "-no_chromatic_correction") == 0)
			{
				perform_chromatic_correction = false;
			}
			else if (strcmp(argv[i], "-monochrome_light") == 0)
			{
				monochrome_plot = true;
			}
			else if (strcmp(argv[i], "-run_minuit_calibration") == 0)
			{
				run_minuit_calibration = true;
			}
			else if (strcmp(argv[i], "-output_box_angles_n") == 0)
			{
				i++;
				output_box_angles_n = atoi(argv[i]);
			}
			else if (strcmp(argv[i], "-prog_sep") == 0)
			{
				use_prog_sep = true;
			}
			else if (strcmp(argv[i], "-slac_run") == 0)
			{
				slac_run = true;
				three_seg_mirror = false;
				mirror_r_difference = 0;
				//mean_n_phot = 31.1;
				mean_n_phot = 32.4;
				spread_n_phot = 6;
				liquid_index = 1.47;
				phi_phots_reduce = 10;
			}
			else if (strcmp(argv[i], "-slac_geometry") == 0)
			{
				three_seg_mirror = false;
				mirror_r_difference = 0;
				//mean_n_phot = 31.1;
				mean_n_phot = 32.4;
				spread_n_phot = 6;
				liquid_index = 1.47;
				phi_phots_reduce = 10;
				//	sm_xl = -300;
				//	sm_xr = sm_xl + 1000;
				miny = -1500;
				maxy = 1500;
				digit_miny = miny;
				digit_maxy = maxy;
			}
			else if (strcmp(argv[i], "-flatten_time") == 0)
			{
				flatten_time = true;
			}
			else if (strcmp(argv[i], "-open_image_plane") == 0)
			{
				pmt_min_z = -1000;
				pmt_max_z = 1000;
				large_mirror_min_z = -1000;
				large_mirror_max_z = 1000;
			}
			else if (strcmp(argv[i], "-mean_n_phot") == 0)
			{
				i++;
				mean_n_phot = atof(argv[i]);
			}
			else if (strcmp(argv[i], "-s_func_t") == 0)
			{
				i++;
				s_func_t = atof(argv[i]);
			}
			else if (strcmp(argv[i], "-t_unc") == 0)
			{
				i++;
				t_unc = atof(argv[i]);
			}
			else if (strcmp(argv[i], "-t_bin_size") == 0)
			{
				i++;
				t_bin_size = atof(argv[i]);
			}
			else if (strcmp(argv[i], "-n_phi_phots") == 0)
			{
				i++;
				n_phi_phots = atoi(argv[i]);
			}
			else if (strcmp(argv[i], "-n_z_phots") == 0)
			{
				i++;
				n_z_phots = atoi(argv[i]);
			}
			else if (strcmp(argv[i], "-spread_n_phot") == 0)
			{
				i++;
				spread_n_phot = atof(argv[i]);
			}
			else if (strcmp(argv[i], "-t") == 0)
			{
				i++;
				time_window = atof(argv[i]);
			}
			else if (strcmp(argv[i], "-particle_x_spread") == 0)
			{
				i++;
				particle_x_spread = atof(argv[i]);
			}
			else if (strcmp(argv[i], "-particle_y_spread") == 0)
			{
				i++;
				particle_y_spread = atof(argv[i]);
			}
			else if (strcmp(argv[i], "-particle_y") == 0)
			{
				i++;
				particle_y = atof(argv[i]);
			}
			else if (strcmp(argv[i], "-particle_x") == 0)
			{
				i++;
				particle_x = atof(argv[i]);
			}
			else if (strcmp(argv[i], "-particle_theta_spread") == 0)
			{
				i++;
				particle_theta_spread = atof(argv[i]);
			}
			else if (strcmp(argv[i], "-energy_spread") == 0)
			{
				i++;
				energy_spread = atof(argv[i]);
			}
			else if (strcmp(argv[i], "-E") == 0)
			{
				i++;
				energy = atof(argv[i]);
				energy_mean = energy;
			}
			else if (strcmp(argv[i], "-n") == 0)
			{
				i++;
				num_runs = atoi(argv[i]);
			}
			else if (strcmp(argv[i], "-max_particles") == 0)
			{
				i++;
				max_particles = atoi(argv[i]);
			}
			else if (strcmp(argv[i], "-pmt_res") == 0)
			{
				i++;
				resx = atof(argv[i]);
				resy = resx;
			}
			else if (strcmp(argv[i], "-liquid_index") == 0)
			{
				i++;
				liquid_index = atof(argv[i]);
			}
			else if (strcmp(argv[i], "-broaden_events") == 0)
			{
				i++;
				broaden_events = atoi(argv[i]);
			}
			else if (strcmp(argv[i], "-wedge_uncertainty") == 0)
			{
				i++;
				wedge_uncertainty = atof(argv[i]);
			}
			else if (strcmp(argv[i], "-refrac_index") == 0)
			{
				i++;
				refrac_index = atof(argv[i]);
			}
			else if (strcmp(argv[i], "-mirror_angle_change") == 0)
			{
				i++;
				mirror_angle_change = atof(argv[i]);
			}
			else if (strcmp(argv[i], "-mirror_angle_change_unc") == 0)
			{
				i++;
				mirror_angle_change_unc = atof(argv[i]);
			}
			else if (strcmp(argv[i], "-mirror_angle_change_yunc") == 0)
			{
				i++;
				mirror_angle_change_yunc = atof(argv[i]);
			}
			else if (strcmp(argv[i], "-box_rot") == 0)
			{
				i++;
				box_rot = atof(argv[i]);
			}
			else if (strcmp(argv[i], "-box_rot_unc") == 0)
			{
				i++;
				box_rot_unc = atof(argv[i]);
			}
			else if (strcmp(argv[i], "-bar_box_box_angle") == 0)
			{
				i++;
				bar_box_box_angle = atof(argv[i]);
			}
			else if (strcmp(argv[i], "-mirror_r_difference") == 0)
			{
				i++;
				mirror_r_difference = atof(argv[i]);
			}
			else if (strcmp(argv[i], "-wedge_non_uniformity") == 0)
			{
				i++;
				wedge_non_uniformity = atof(argv[i]);
			}
			else if (strcmp(argv[i], "-pmt_offset") == 0)
			{
				i++;
				pmt_offset = atof(argv[i]);
			}
			else if (strcmp(argv[i], "-main_mirror_nonuniformity") == 0)
			{
				i++;
				main_mirror_nonuniformity = atof(argv[i]);
			}
			else if (strcmp(argv[i], "-upper_wedge_yang_spread") == 0)
			{
				i++;
				upper_wedge_yang_spread = atof(argv[i]);
			}
			else if (strcmp(argv[i], "-main_mirror_angle_off") == 0)
			{
				i++;
				main_mirror_angle_off = atof(argv[i]);
			}
			else if (strcmp(argv[i], "-main_mirror_yangle_off") == 0)
			{
				i++;
				main_mirror_yangle_off = atof(argv[i]);
			}
			else if (strcmp(argv[i], "-main_mirror_zangle_off") == 0)
			{
				i++;
				main_mirror_zangle_off = atof(argv[i]);
			}
			else if (strcmp(argv[i], "-main_mirror_yoff") == 0)
			{
				i++;
				main_mirror_yoff = atof(argv[i]);
			}
			else if (strcmp(argv[i], "-main_mirror_zoff") == 0)
			{
				i++;
				main_mirror_zoff = atof(argv[i]);
			}
			else if (strcmp(argv[i], "-rseed") == 0)
			{
				i++;
				rseed = atof(argv[i]);
			}
			else if (strcmp(argv[i], "-upper_wedge_yang_spread") == 0)
			{
				i++;
				upper_wedge_yang_spread = atof(argv[i]);
			}
			else if (strcmp(argv[i], "-box_check_n") == 0)
			{
				i++;
				box_check_n = atoi(argv[i]);
			}
			else if (strcmp(argv[i], "-box_check_theta") == 0)
			{
				i++;
				box_check_theta = atof(argv[i]);
			}
			else if (strcmp(argv[i], "-box_check_phi") == 0)
			{
				i++;
				box_check_phi = atof(argv[i]);
			}
			else if (strcmp(argv[i], "-box_check_theta_unc") == 0)
			{
				i++;
				box_check_theta_unc = atof(argv[i]);
			}
			else if (strcmp(argv[i], "-box_check_phi_unc") == 0)
			{
				i++;
				box_check_phi_unc = atof(argv[i]);
			}
			else if (strcmp(argv[i], "-box_check_overall_theta") == 0)
			{
				i++;
				box_check_overall_theta = atof(argv[i]);
			}
			else if (strcmp(argv[i], "-phot_check_n") == 0)
			{
				i++;
				phot_check_n = atoi(argv[i]);
			}
			else if (strcmp(argv[i], "-phot_check_max_theta") == 0)
			{
				i++;
				phot_check_max_theta = atof(argv[i]);
			}
			else if (strcmp(argv[i], "-bar_box_xoff") == 0)
			{
				i++;
				bar_box_xoff = atof(argv[i]);
			}
			else if (strcmp(argv[i], "-bar_box_yoff") == 0)
			{
				i++;
				bar_box_yoff = atof(argv[i]);
			}
			else if (strcmp(argv[i], "-bar_box_zoff") == 0)
			{
				i++;
				bar_box_zoff = atof(argv[i]);
			}
			else
			{
				printf("Unrecognized argument: %s\n",argv[i]);
			}
		}
	}


	float main_mirror_angle = 74.11+mirror_angle_change;

	//main_mirror_angle += .23;

	float rad_to_deg = 57.2958;

	//float res_enhance = 1/6.;
	float res_enhance = 1;

//	float prog_thresh = 500;

	if (flatten_time == true)
	{
		s_func_t = 100000000;
	}	

	float outcsv_x,outcsv_y;
	//float outcsv_t;
	outcsv_x = 0*35;//bars are 35mm wide
	outcsv_y = 0;//mm
	//outcsv_t = 0;//ns


	if(out_csv){
		n_phi_phots = 4000;
		n_z_phots = 8;
		particle_y += outcsv_y;}
	float pdf_unc_red_fac = 1;


	TRandom3 spread_ang(rseed+3);
/*
	DircOpticalSim *dirc_model_2 = new DircOpticalSim(\
			rseed,\
			-1200 + mirror_r_difference,\
			foc_mirror_size,\
			main_mirror_angle,\
			600,\
			47.87 + box_rot + mirror_angle_change);
*/

//hackish, especially with the uncertainty functions, but easiest
#ifdef USE_BABAR_BOX
	DircBaBarSim *dirc_model = new DircBaBarSim(\
                        rseed);
	minx = -429.0/2;//assum it's babar's size
	minx = -1500;
	maxx = -minx;
	miny = 0;
	maxy = dirc_model->get_sens_r()*dirc_model->get_sens_subtend_angle()/57.3;
	float babar_t_unc = 1;
	float babar_t_bin = 1;
	float babar_pmt_r = 29./2;
#else
	DircThreeSegBoxSim *dirc_model = new DircThreeSegBoxSim(\
			rseed,\
			-1200 + mirror_r_difference,\
			foc_mirror_size,\
			main_mirror_angle,\
			600,\
			47.87 + box_rot + mirror_angle_change);
#endif

	dirc_model->set_store_traveled(false);// uses LOTS of memory if set to true.
	dirc_model->set_liquid_index(liquid_index);
	dirc_model->set_wedge_mirror_rand(wedge_non_uniformity);
	dirc_model->set_three_seg_mirror(three_seg_mirror);
	dirc_model->set_kaleidoscope_plot(kaleidoscope_plot);	
	dirc_model->set_pmt_plane_zs(pmt_min_z,pmt_max_z);
	dirc_model->set_large_mirror_zs(large_mirror_min_z,large_mirror_max_z);
	dirc_model->set_use_quartz_n_for_liquid(use_quartz_for_liquid);


	float muon_beta, pion_beta, kaon_beta/*, electron_beta:=1*/;
	muon_beta=pion_beta=kaon_beta=-1;
	float muon_angle, pion_angle, kaon_angle;
	muon_angle=pion_angle=kaon_angle = -1;


	TFile* tfile = new TFile(rootfilename,"RECREATE");

	TH1F *ll_diff_pion = new TH1F("ll_diff_pion","Difference of log likelihood real = pion",200000,-200,200);
	TH1F *ll_diff_kaon = new TH1F("ll_diff_kaon","Difference of log likelihood real = kaon",200000,-200,200);
	TH1F *ll_diff_pion_up = new TH1F("ll_diff_pion_up","Difference of log likelihood real = pion \"up\" going photons",200000,-200,200);
	TH1F *ll_diff_kaon_up = new TH1F("ll_diff_kaon_up","Difference of log likelihood real = kaon \"up\" going photons",200000,-200,200);
	TH1F *ll_diff_pion_down = new TH1F("ll_diff_pion_down","Difference of log likelihood real = pion \"down\" going photons",200000,-200,200);
	TH1F *ll_diff_kaon_down = new TH1F("ll_diff_kaon_down","Difference of log likelihood real = kaon \"down\" going photons",200000,-200,200);
	TH1F *phot_found_pion = new TH1F("phot_found_pion","number of photons found on pion angle", 1001,-.5,1000.5);
	TH1F *phot_found_kaon = new TH1F("phot_found_kaon","number of photons found on kaon angle", 1001,-.5,1000.5);
	TH1F *phot_found_pion_up = new TH1F("phot_found_pion_up","number of photons found on pion angle \"up\" going photons", 1001,-.5,1000.5);
	TH1F *phot_found_kaon_up = new TH1F("phot_found_kaon_up","number of photons found on kaon angle \"up\" going photons", 1001,-.5,1000.5);
	TH1F *phot_found_pion_down = new TH1F("phot_found_pion_down","number of photons found on pion angle \"down\" going photons", 1001,-.5,1000.5);
	TH1F *phot_found_kaon_down = new TH1F("phot_found_kaon_down","number of photons found on kaon angle \"down\" going photons", 1001,-.5,1000.5);

	TH1F *ref_pion_before = new TH1F("ref_pion_before","Angle of Pion Photons going into interface", 9000,0,90);
	TH1F *ref_kaon_before = new TH1F("ref_kaon_before","Angle of Kaon Photons going into interface", 9000,0,90);
	TH1F *ref_pion_after = new TH1F("ref_pion_after","Angle of Pion Photons going out of interface", 9000,0,90);
	TH1F *ref_kaon_after = new TH1F("ref_kaon_after","Angle of Kaon Photons going out of interface", 9000,0,90);
	TH1F *ref_pion_sens_plane = new TH1F("ref_pion_sens_plane","Angle of Pion Photons going into window to PMTs", 9000,0,90);
	TH1F *ref_kaon_sens_plane = new TH1F("ref_kaon_sens_plane","Angle of Kaon Photons going into window to PMTs", 9000,0,90);

	TH1F *n_bounces_x = new TH1F("n_bounces_x","Number of bounces in x direction",1000,-.5,999.5);
	TH1F *n_bounces_z = new TH1F("n_bounces_z","Number of bounces in z direction",1000,-.5,999.5);
	TH1F *n_bounces_x_direct = new TH1F("n_bounces_x_direct","Number of bounces in x direction (direct photons)",1000,-.5,999.5);
	TH1F *n_bounces_z_direct = new TH1F("n_bounces_z_direct","Number of bounces in z direction (direct photons)",1000,-.5,999.5);
	TH1F *n_bounces_x_indirect = new TH1F("n_bounces_x_indirect","Number of bounces in x direction (indirect photons)",1000,-.5,999.5);
	TH1F *n_bounces_z_indirect = new TH1F("n_bounces_z_indirect","Number of bounces in z direction (indirect photons)",1000,-.5,999.5);



	TH1F *pion_dist_x = new TH1F("pion_dist_x","x val of intercepted points - pion",(maxx-minx)/(res_enhance*resx),minx,maxx);
	TH1F *pion_dist_y = new TH1F("pion_dist_y","y val of intercepted points - pion",(maxy-miny)/(res_enhance*resy),miny,maxy);
	TH1F *pion_dist_t = new TH1F("pion_dist_t","t val of intercepted points - pion",(maxt-mint)/(res_enhance*rest),mint,maxt);
	TH1F *kaon_dist_x = new TH1F("kaon_dist_x","x val of intercepted points - kaon",(maxx-minx)/(res_enhance*resx),minx,maxx);
	TH1F *kaon_dist_y = new TH1F("kaon_dist_y","y val of intercepted points - kaon",(maxy-miny)/(res_enhance*resy),miny,maxy);
	TH1F *kaon_dist_t = new TH1F("kaon_dist_t","t val of intercepted points - kaon",(maxt-mint)/(res_enhance*rest),mint,maxt);

	int d_bins = 2000;
	float d_bin_min = -200;
	float d_bin_max = -d_bin_min;

	TH1F *pion_midline_dx = new TH1F("pion_midline_dx","dx from midline for pions",d_bins,d_bin_min,d_bin_max);
	TH1F *pion_midline_dy = new TH1F("pion_midline_dy","dy from midline for pions",d_bins,d_bin_min,d_bin_max);
	TH1F *pion_midline_dt = new TH1F("pion_midline_dt","dt from midline for pions",d_bins,d_bin_min,d_bin_max);
	TH1F *kaon_midline_dx = new TH1F("kaon_midline_dx","dx from midline for kaons",d_bins,d_bin_min,d_bin_max);
	TH1F *kaon_midline_dy = new TH1F("kaon_midline_dy","dy from midline for kaons",d_bins,d_bin_min,d_bin_max);
	TH1F *kaon_midline_dt = new TH1F("kaon_midline_dt","dt from midline for kaons",d_bins,d_bin_min,d_bin_max);

	TH1F *pion_lut_vals = new TH1F("pion_lut_vals","Pion angle from LUT (deg)",4500,0,90);
	TH1F *kaon_lut_vals = new TH1F("kaon_lut_vals","Kaon angle from LUT (deg)",4500,0,90);
	TH1F *pion_lut_means = new TH1F("pion_lut_means","Pion reconstructed angle from LUT (deg)",5000,45,50);
	TH1F *kaon_lut_means = new TH1F("kaon_lut_means","Kaon reconstructed angle from LUT (deg)",5000,45,50);
	TH1F *pion_lut_angles = new TH1F("pion_lut_angles","Number of angles from the LUT for pion event",5001,-.05,5000.5);
	TH1F *kaon_lut_angles = new TH1F("kaon_lut_angles","Number of angles from the LUT for kaon event",5001,-.05,5000.5);

	TH2F *pion_lut_dt_v_dang = new TH2F("pion_lut_dt_v_dang","Normed Time difference (ns/m) versus angle error (mrad) for Pion LUT",1000,-50,50,800,-2,2);
	TH2F *kaon_lut_dt_v_dang = new TH2F("kaon_lut_dt_v_dang","Normed Time difference (ns/m) versus angle error (mrad) for Kaon LUT",1000,-50,50,800,-2,2);

	TH1F *box_check_x = new TH1F("box_check_x","Check Box X",(maxx-minx)/(.1*res_enhance*resx),minx,maxx);
	TH1F *box_check_y = new TH1F("box_check_y","Check Box Y",(maxy-miny)/(.1*res_enhance*resy),miny,maxy);
	TH1F *box_check_d = new TH1F("box_check_d","Check Box Dist",800000,0,8000);

	TH1F *pion_dist_geant_x = new TH1F("pion_dist_geant_x","x val of intercepted points (from geant) - pion",(maxx-minx)/(res_enhance*resx),minx,maxx);
	TH1F *pion_dist_geant_y = new TH1F("pion_dist_geant_y","y val of intercepted points (from geant) - pion",(maxy-miny)/(res_enhance*resy),miny,maxy);
	TH1F *pion_dist_geant_t = new TH1F("pion_dist_geant_t","t val of intercepted points (from geant) - pion",(maxt-mint)/(res_enhance*rest),mint,maxt);

	TH2F *ref_theta_cphi_pion = new TH2F("ref_theta_cphi_pion","Emitted angle of Pion Photons versus angle into interface", 7200,0,360,1800,0,90);
	TH2F *ref_theta_cphi_kaon = new TH2F("ref_theta_cphi_kaon","Emitted angle of Kaon Photons versus angle into interface", 7200,0,360,1800,0,90);

	TH2F *pion_dist_xy = new TH2F("pion_dist_xy","xy val of intercepted points - pion",(maxx-minx)/(res_enhance*resx),minx,maxx,(maxy-miny)/(res_enhance*resy),miny,maxy);
	TH2F *pion_dist_xy_geantbins_cm = new TH2F("pion_dist_xy_geantbins_cm","xy val of intercepted points with geant bins - pion",500,-150.,150.,67,-5.,35.2);
	TH2F *pion_dist_xy_geantbins = new TH2F("pion_dist_xy_geantbins","xy val of intercepted points with geant bins - pion",500,-1500,1500,67,-50,352);
	//TH2F *pion_dist_xy_geantbins = new TH2F("pion_dist_xy_geantbins","xy val of intercepted points - pion",500,-1500,1500,67,-150,252);
	TH2F *kaon_dist_xy = new TH2F("kaon_dist_xy","xy val of intercepted points - kaon",(maxx-minx)/(res_enhance*resx),minx,maxx,(maxy-miny)/(res_enhance*resy),miny,maxy);

	TH2F *box_check_xy = new TH2F("box_check_xy","xy val of points from top of wedge",(maxx-minx)/(res_enhance*resx),minx,maxx,(maxy-miny)/(res_enhance*resy),miny,maxy);

	TH2F *pion_dist_geant_xy = new TH2F("pion_dist_geant_xy","xy val of intercepted points (from geant) - pion",(maxx-minx)/(res_enhance*resx),minx,maxx,(maxy-miny)/(res_enhance*resy),miny,maxy);

	TH2F *pion_dist_xt = new TH2F("pion_dist_xt","xt val of intercepted points - pion",(maxx-minx)/(res_enhance*resx),minx,maxx,(maxt-mint)/(res_enhance*rest),mint,maxt);
	TH2F *kaon_dist_xt = new TH2F("kaon_dist_xt","xt val of intercepted points - kaon",(maxx-minx)/(res_enhance*resx),minx,maxx,(maxt-mint)/(res_enhance*rest),mint,maxt);

	TH2F *pion_dist_geant_xt = new TH2F("pion_dist_geant_xt","xt val of intercepted points (from geant) - pion",(maxx-minx)/(res_enhance*resx),minx,maxx,(maxt-mint)/(res_enhance*rest),mint,maxt);

	TH2F *pion_dist_yt = new TH2F("pion_dist_yt","yt val of intercepted points - pion",(maxy-miny)/(res_enhance*resy),miny,maxy,(maxt-mint)/(res_enhance*rest),mint,maxt);
	TH2F *kaon_dist_yt = new TH2F("kaon_dist_yt","yt val of intercepted points - kaon",(maxy-miny)/(res_enhance*resy),miny,maxy,(maxt-mint)/(res_enhance*rest),mint,maxt);

	TH2F *pion_dist_geant_yt = new TH2F("pion_dist_geant_yt","yt val of intercepted points (from geant) - pion",(maxy-miny)/(res_enhance*resy),miny,maxy,(maxt-mint)/(res_enhance*rest),mint,maxt);


	maxy *= 5;

	TH2F *pion_coverage_xy = new TH2F("pion_coverage_xy","xy val of generated points - pion",(maxx-minx)/(res_enhance*resx),minx,maxx,(maxy-miny)/(res_enhance*resy),miny,maxy);
	TH2F *kaon_coverage_xy = new TH2F("kaon_coverage_xy","xy val of generated points - kaon",(maxx-minx)/(res_enhance*resx),minx,maxx,(maxy-miny)/(res_enhance*resy),miny,maxy);

	TH2F *kinematics_yields = new TH2F("kinematics_yields","Fractional Yields versus angles",(kinematics_yields_max_theta-kinematics_yields_min_theta)/kinematics_yields_step_theta*1.00001,\
			kinematics_yields_min_theta,\
			kinematics_yields_max_theta,\
			(kinematics_yields_max_phi-kinematics_yields_min_phi)/kinematics_yields_step_phi*1.00001,\
			kinematics_yields_min_phi,\
			kinematics_yields_max_phi);

	TH1F *pion_cerenkov = new TH1F("pion_cerenkov","pion cerenkov angle distribution",300,45,48);
	TH1F *kaon_cerenkov = new TH1F("kaon_cerenkov","kaon cerenkov angle distribution",300,45,48);
	TH1F *pion_lambda = new TH1F("pion_lambda","pion wavelength distribution",450,250,700);
	TH1F *kaon_lambda = new TH1F("kaon_lambda","kaon wavelength distribution",450,250,700);

	TH1F *upper_wedge_reflect_angle = new TH1F("upper_wedge_reflect_angles","angles reflected from the upper wedge",180,0,90);
	TH1F *focus_reflect_angle = new TH1F("focus_reflect_angles","angles reflected from the focusing mirror",180,0,90);
	TH1F *large_flat_reflect_angle = new TH1F("large_flat_reflect_angles","angles reflected from the large flat mirror",180,0,90);
	TH1F *side_reflect_angle = new TH1F("side_reflect_angles","angle photon makes with the side mirror (if it were to reflect)",180,0,90);

	TH1F *pion_unused_photons = new TH1F("pion_unused_photons","pion photonons not used on pion data, pion pdf",100,-.5,99.5);
	TH1F *kaon_unused_photons = new TH1F("kaon_unused_photons","kaon photonons not used on kaon data, kaon pdf",100,-.5,99.5);

	TH1F *liquid_dist = new TH1F("liquid_dist","distance travelled in liquid (mm)",1500,0,1500);

	TH1F *simulation_time = new TH1F("simulation_time","Simulation time per particle",4001,-.5,4000.5);

	TH1F *zero_pion_id = new TH1F("zero_pion_id","Correct pion ID with cut at 0",2,-0.5,1.5);
	TH1F *zero_kaon_id = new TH1F("zero_kaon_id","Correct kaon ID with cut at 0",2,-0.5,1.5);

#ifdef USE_BABAR_BOX
	DircBaBarDigitizer digitizer(\
			minx,\
			maxx,\
			miny,\
			maxy,\
			babar_pmt_r,\
			babar_t_unc,\
			babar_t_bin);
			
#else
	DircRectDigitizer digitizer(\
			minx,\
			maxx,\
			resx,\
			digit_miny,\
			digit_maxy,\
			resy,\
			t_unc,\
			t_bin_size);
#endif
	//  	
	printf("Beginning Run\n");
	float llc, llf, ll_diff;
	llc=llf=ll_diff=0;
	int pion_kaon_diff = 0;
	std::vector<dirc_point> sim_points;
	std::vector<dirc_point> confound_points;
	dirc_model->set_focmirror_nonuniformity(main_mirror_nonuniformity);
	/*   read in file for events and 
	 *   initialize an array 
	 *    */
	if (inputrootfile)
	{
		printf("Opening root file: %s\n",inputrootfilename);
		//Gross... I got root types in my loops...
		Float_t px[10];
		Float_t py[10];
		Float_t pz[10];
		Float_t particle_x[10];
		Float_t particle_y[10];
		Float_t particle_z[10];
		Float_t gen_t[10];
		Float_t gen_E[10];
		Float_t hitw[2000];
		Float_t hity[2000];
		Float_t hit_t[2000];
		int nHit;
		Int_t particle_id[10];

		TFile* inputrootfile = new TFile(inputrootfilename,"READONLY");
		TTree *in_tree = (TTree*) inputrootfile->Get("dircTreePID");

		int nentries = in_tree->GetEntries();

		in_tree->SetBranchAddress("ParticleID",&particle_id);
		in_tree->SetBranchAddress("genTrk_px",&px);
		in_tree->SetBranchAddress("genTrk_py",&py);
		in_tree->SetBranchAddress("genTrk_pz",&pz);
		in_tree->SetBranchAddress("genTrk_x",&particle_x);
		in_tree->SetBranchAddress("genTrk_y",&particle_y);
		in_tree->SetBranchAddress("genTrk_z",&particle_z);
		in_tree->SetBranchAddress("genTrk_t",&gen_t);
		in_tree->SetBranchAddress("genTrk_E",&gen_E);

		in_tree->SetBranchAddress("hitpixel_w",&hitw);
		in_tree->SetBranchAddress("hitpixel_y",&hity);
		in_tree->SetBranchAddress("hit_t",&hit_t);

		in_tree->SetBranchAddress("nHit",&nHit);

		std::vector<std::vector<dirc_point> > root_hits_array;

		unsigned int r = 0;

		int iPID, iBAR;
		//int ievent_index;
		float ix,iy,it,itheta,iphi,iE;
		std::vector<int> PID,BAR;
		std::vector<int> event_index;
		std::vector<float> x,y,t,theta,phi,E;

		float min_t = 1000;
		float avg_x = 0;
		float avg_y = 0;
		float avg_t = 0;
		int tot_nhit = 0;
		for (int i = 0; i < nentries; i++) {

			if (r >= (unsigned int) max_particles) break;

			in_tree->GetEntry(i);
			// 			if (iE > 4) continue;	
			//TODO read more than 1 particle per event
			iPID = particle_id[0];
			iBAR = dirc_model->get_bar_from_x(particle_x[0]*10);

			if (iBAR != 6) {continue;}
			//			if (iPID == 8) {continue;}
			ix = particle_x[0]*10 - dirc_model->get_bar_offset(iBAR);
			iy = particle_y[0]*10 + 700;//Figure out that geometry
			it = gen_t[0] + 2*r;//2ns spacing of events at least.
			iE = gen_E[0];
			float p_trans = sqrt(py[0]*py[0]+px[0]*px[0]);
			itheta = atan(p_trans/pz[0])*180/3.14159;
			iphi = atan2(py[0],px[0])*180/3.14159;

			if (force_kinematics == true)
			{
				itheta = particle_theta;
				iphi = particle_phi;
				iE = energy;
			}
			event_index.push_back(i);
			PID.push_back(iPID);
			BAR.push_back(iBAR);
			x.push_back(ix);
			y.push_back(iy);
			t.push_back(it);
			theta.push_back(itheta);
			phi.push_back(iphi);
			E.push_back(iE);

			if (t[r] < 0) printf("low t: %12.04f\n",t[r]);

			std::vector<dirc_point> root_hits_add;
			if (iPID == 8 || iPID == 9){
				tot_nhit += nHit;
			}
			for (int j = 0; j < nHit; j++)
			{
				dirc_point addpt;
				//these are swapped knowingly
				addpt.x = hity[j]*10 + 151.4 + 22;
				addpt.y = hitw[j]*10 - 28.85 + 25.1;
				//Strange extra hits at 149.15 for some reason
				if (addpt.y > 149.1 && addpt.y < 149.2) continue;
				//addpt.t = t[r] + hit_t[j] - 11;
				addpt.t = hit_t[j] + 2*r;

				if (addpt.t - 2*r < 0)
				{//hmmmm can't be - ignore it
					printf("low addt: %12.04f %d %d\n", addpt.t - 2*r, i, j);
					if (iPID == 8 || iPID == 9)
					{
						tot_nhit--;
					}
					continue;
				}

				root_hits_add.push_back(addpt);



				if (fill_distributions == true)
				{
					addpt.t -= 2*r;
					pion_dist_geant_x->Fill(addpt.x);
					pion_dist_geant_y->Fill(addpt.y);
					pion_dist_geant_t->Fill(addpt.t);
					pion_dist_geant_xy->Fill(addpt.x,addpt.y);
					pion_dist_geant_xt->Fill(addpt.x,addpt.t);
					pion_dist_geant_yt->Fill(addpt.y,addpt.t);
					addpt.t += 2*r;
				}

				if (addpt.t < 0) printf("low addpt: %12.04f\n",addpt.t);
				if ((iPID == 8 || iPID == 9))
				{
					min_t = std::min(min_t,addpt.t - 2*r);
					avg_x += addpt.x;
					avg_y += addpt.y;
					avg_t += addpt.t - 2*r;
				}
			}
			root_hits_array.push_back(root_hits_add);

			r++;
		}
		float high_avg = 0;
		float low_avg = 0;
		float high_min = 10000000;
		float low_min = 10000000;
		int high_count = 0;
		int low_count = 0;
		float tmp_t = -1;
		for (unsigned int j = 0; j < root_hits_array.size(); j++)
		{
			for (unsigned int k = 0; k < root_hits_array[j].size(); k++)
			{
				tmp_t = root_hits_array[j][k].t - 2*j;
				if (tmp_t < 0) printf("low t: %12.04f %d %d\n",tmp_t,j,k);
				if (tmp_t > avg_t/tot_nhit)
				{
					high_count++;
					high_avg += tmp_t;
					high_min = std::min(high_min, tmp_t);
				}
				else
				{
					low_count++;
					low_avg += tmp_t;
					low_min = std::min(low_min, tmp_t);
				}
			}
		}


		printf("Avg hit x y t: %12.04f %12.04f %12.04f\n",avg_x/tot_nhit,avg_y/tot_nhit, avg_t/tot_nhit);
		printf("min input t delta: %12.04f\n",min_t);
		printf("T abv avg : below avg : diff : %12.04f %12.04f %12.04f\n", high_avg/high_count, low_avg/low_count, high_avg/high_count - low_avg/low_count);
		printf("T abv min : below min : diff : %12.04f %12.04f %12.04f\n", high_min, low_min, high_min - low_min);
		inputrootfile->Close("R");

		int mc_tally=0;
//		int confounded_tally=0;
		num_runs=0;
		coverage_plot=false;

		//declare memory addresses for inputs
		//units: mm, ns, deg, GeV 

		//declare for beta and angle based on PID, and constructors for hitpoints and pdfs
		float pion_mc_beta,kaon_mc_beta,pion_mc_angle,kaon_mc_angle;
		std::vector<dirc_point> hits_trk_is_pion;
		std::vector<dirc_point> hits_trk_is_kaon;

		float min_pval = 1e-2;
		DircSpreadGaussian * pdf_as_pion = new DircSpreadGaussian(\
				sfunc_sig,\
				hits_trk_is_pion,\
				s_func_x,\
				s_func_y,\
				s_func_t,\
				min_pval);
		DircSpreadGaussian * pdf_as_kaon = new DircSpreadGaussian(\
				sfunc_sig,\
				hits_trk_is_kaon,\
				s_func_x,\
				s_func_y,\
				s_func_t,\
				min_pval);//could make this an array that is filled...
		//Might be worth making a new copy each time
//		DircProbabilitySeparation * sep_pdfs_mc;
/*
		DircProgressiveSeparation *progressive_separation = \
								    new DircProgressiveSeparation(\
										    dirc_model,\
										    n_phi_phots*n_z_phots,\
										    n_step_phots,\
										    sfunc_sig,\
										    s_func_x,\
										    s_func_y,\
										    s_func_t,\
										    kmass,\
										    pimass,\
										    prog_thresh);
*/
		//		r = 100;

		for(unsigned int n=0;n < r; n++){

			min_t = 1000;

			if(fabs(PID[n])>7 && fabs(PID[n])<13)//8 and 9 are pions, 11 and 12 are kaons
			{

				mc_tally++;

				//
				//and replace acos......
				//cdd replace 1.47 with float refrac_index=1.47;
				pion_mc_beta = dirc_model->get_beta(E[n],pimass);
				kaon_mc_beta = dirc_model->get_beta(E[n],kmass);
				pion_mc_angle = rad_to_deg*acos(1/(1.47*pion_mc_beta));
				kaon_mc_angle = rad_to_deg*acos(1/(1.47*kaon_mc_beta));

				//sep_pdfs_mc = new DircProbabilitySeparation(pdf_as_kaon,pdf_as_pion);

				dirc_model->set_focus_mirror_angle(\
						spread_ang.Gaus(main_mirror_angle,mirror_angle_change_unc),\
						spread_ang.Gaus(0,mirror_angle_change_yunc));
				dirc_model->set_upper_wedge_angle_diff(\
						spread_ang.Gaus(0,wedge_uncertainty),\
						spread_ang.Gaus(0,upper_wedge_yang_spread));
				printf("\r                                                                                                      ");
				printf("\rRunning monte carlo event  %i/%i. Found %i pions/kaons", n,r-1,mc_tally);

				fflush(stdout);

				//simulate pmt hits from primary mc event to compare with preconstructed pdfs
				//


				if(abs(PID[n])==8 || PID[n]==9 ){
					dirc_model->sim_rand_n_photons(\
							sim_points,\
							n_sim_phots,\
							pion_mc_angle,\
							BAR[n],\
							x[n],\
							y[n],\
							0,\
							theta[n],\
							phi[n],\
							tracking_unc,\
							ckov_unc,\
							pion_mc_beta);
				//TODO: Check timing
					pion_kaon_diff = 1;
				}
				else{
					dirc_model->sim_rand_n_photons(\
							sim_points,\
							n_sim_phots,\
							kaon_mc_angle,\
							BAR[n],\
							x[n],\
							y[n],\
							0,\
							theta[n],\
							phi[n],\
							tracking_unc,\
							ckov_unc,\
							kaon_mc_beta);
				//TODO: Check timing
					pion_kaon_diff = -1;
				}


				for (unsigned int i = 0; i < sim_points.size(); i++)
				{
					sim_points[i].t += t[n];
					// 					printf("%12.04f\n",sim_points[i].t);
				}

				//Begin analysis

				//Zero out uncertainties (assume we are starting from calibrated versions

				//				digitizer.digitize_points(sim_points);

				//Ignore all the simpoints - use the root file
				sim_points = root_hits_array[n];

			timing_clock = clock();

			if (use_prog_sep == true)
			{
				dirc_model->set_focus_mirror_angle(\
						spread_ang.Gaus(main_mirror_angle,0),\
						spread_ang.Gaus(0,0));
				dirc_model->set_upper_wedge_angle_diff(\
						spread_ang.Gaus(0,0),\
						spread_ang.Gaus(0,0));
			}
			else
			{
				dirc_model->set_focus_mirror_angle(\
						spread_ang.Gaus(main_mirror_angle,0),\
						spread_ang.Gaus(0,0));
				dirc_model->set_upper_wedge_angle_diff(\
						spread_ang.Gaus(0,0),\
						spread_ang.Gaus(0,0));

				dirc_model->sim_reg_n_photons(\
						hits_trk_is_pion,\
						n_phi_phots,\
						n_z_phots,\
						pion_mc_angle,\
						BAR[n],\
						x[n],\
						y[n],\
						0,
						theta[n],\
						phi[n],\
						0,\
						ckov_unc/pdf_unc_red_fac,\
						pion_mc_beta);
				//TODO: Check timing

				dirc_model->sim_reg_n_photons(\
						hits_trk_is_kaon,\
						n_phi_phots,\
						n_z_phots,\
						kaon_mc_angle,
						BAR[n],\
						x[n],\
						y[n],\
						0,
						theta[n],\
						phi[n],\
						0,\
						ckov_unc/pdf_unc_red_fac,\
						kaon_mc_beta);
				//TODO: Check timing
				float avg_x = 0;	
				float avg_y = 0;	
				float avg_t = 0;	
				for (unsigned int k = 0; k < hits_trk_is_pion.size(); k++)
				{
					hits_trk_is_pion[k].t += t[n];
					min_t = std::min(min_t,hits_trk_is_pion[k].t - 2*n);
					avg_x += hits_trk_is_pion[k].x;
					avg_y += hits_trk_is_pion[k].y;
					avg_t += hits_trk_is_pion[k].t - 2*n;
				}

				float high_avg = 0;
				float low_avg = 0;
				float high_min = 10000000;
				float low_min = 10000000;
				int high_count = 0;
				int low_count = 0;
				float tmp_t = -1;
				float tmp_x = -1;
				float tmp_y = -1;
				for (unsigned int k = 0; k < hits_trk_is_pion.size(); k++)
				{
					tmp_t = hits_trk_is_pion[k].t - 2*n;
					tmp_x = hits_trk_is_pion[k].x;
					tmp_y = hits_trk_is_pion[k].y;
					if (n == 0 && fill_distributions == true)
					{
						pion_dist_x->Fill(tmp_x);
						pion_dist_y->Fill(tmp_y);
						pion_dist_t->Fill(tmp_t);
						pion_dist_xy->Fill(tmp_x, tmp_y);
						pion_dist_xt->Fill(tmp_x, tmp_t);
						pion_dist_yt->Fill(tmp_t, tmp_t);
					}

					if (tmp_t > avg_t/hits_trk_is_pion.size())
					{
						high_count++;
						high_avg += tmp_t;
						high_min = std::min(high_min,tmp_t);
					}
					else
					{
						low_count++;
						low_avg += tmp_t;
						low_min = std::min(low_min,tmp_t);
					}
				}
				printf("\nAvg t High : low : diff : %12.04f : %12.04f : % 12.04f\n",high_avg/high_count,low_avg/low_count,high_avg/high_count-low_avg/low_count);
				printf("Min t High : low : diff : %12.04f : %12.04f : % 12.04f\n",high_min,low_min,high_min - low_min);

				for (unsigned int k = 0; k < hits_trk_is_kaon.size(); k++)
				{
					hits_trk_is_kaon[k].t += t[n];
				}
				for (unsigned int k = 0; k < hits_trk_is_kaon.size(); k++)
				{
					tmp_t = hits_trk_is_kaon[k].t - 2*n;
					tmp_x = hits_trk_is_kaon[k].x;
					tmp_y = hits_trk_is_kaon[k].y;
					if (n == 0 && fill_distributions == true)
					{
						kaon_dist_x->Fill(tmp_x);
						kaon_dist_y->Fill(tmp_y);
						kaon_dist_t->Fill(tmp_t);
						kaon_dist_xy->Fill(tmp_x, tmp_y);
						kaon_dist_xt->Fill(tmp_x, tmp_t);
						kaon_dist_yt->Fill(tmp_t, tmp_t);
					}

					if (tmp_t > avg_t/hits_trk_is_kaon.size())
					{
						high_count++;
						high_avg += tmp_t;
						high_min = std::min(high_min,tmp_t);
					}
					else
					{
						low_count++;
						low_avg += tmp_t;
						low_min = std::min(low_min,tmp_t);
					}
				}
				pdf_as_pion->set_support(hits_trk_is_pion);
				pdf_as_kaon->set_support(hits_trk_is_kaon);


				llc = pdf_as_pion->get_log_likelihood(sim_points);
				llf = pdf_as_kaon->get_log_likelihood(sim_points);

				ll_diff = llc - llf;

			}

			timing_clock = clock() - timing_clock;

			simulation_time->Fill(((float)timing_clock)/(CLOCKS_PER_SEC)*1000);


			if(abs(PID[n])==8||PID[n]==9){
				ll_diff_pion->Fill(ll_diff);

				if (ll_diff > 0)
				{
					zero_pion_id->Fill(1);
				}
				else
				{
					zero_pion_id->Fill(0);
				}

				phot_found_pion->Fill(sim_points.size());
			}
			else{

				ll_diff_kaon->Fill(ll_diff);
				if (ll_diff < 0)
				{
					zero_kaon_id->Fill(1);
				}
				else
				{
					zero_kaon_id->Fill(0);
				}
				phot_found_kaon->Fill(sim_points.size());
			}
			}
		}
		printf("\n");//ensure we're on a new line at the end
		printf("Kaon missID: %12.04f%% \n",(float) 100.0*zero_kaon_id->GetBinContent(1)/(zero_kaon_id->GetBinContent(1)+zero_kaon_id->GetBinContent(2)));
		printf("Pion missID: %12.04f%% \n",(float) 100.0*zero_pion_id->GetBinContent(1)/(zero_pion_id->GetBinContent(1)+zero_pion_id->GetBinContent(2)));
	}
	else if(inputfile){

		printf("Csv file with name %s selected.\n",in_str);
		printf("Opening %s with time window of %fns\n",in_str,time_window);

		int mc_tally=0;
		int confounded_tally=0;
		num_runs=0;
		coverage_plot=false;
		std::ifstream f(in_str);
		if(!f.is_open()){printf("%s not found!",in_str);return -1;}	

		int line_buf_size = 1000; 
		char  s[line_buf_size];
		f.getline(s,line_buf_size); //skips first line
		//declare memory addresses for inputs
		//units: mm, ns, deg, GeV 
		int iPID, iBAR,ievent_index;
		float ix,iy,it,itheta,iphi,iE;
		std::vector<int> PID,BAR,event_index;
		std::vector<float> x,y,t,theta,phi,E;

		//declare for beta and angle based on PID, and constructors for hitpoints and pdfs
		float pion_mc_beta,kaon_mc_beta,pion_mc_angle,kaon_mc_angle;
		std::vector<dirc_point> hits_trk_is_pion;
		std::vector<dirc_point> hits_trk_is_kaon;

		float min_pval = 1e-2;
		DircSpreadGaussian * pdf_as_pion = new DircSpreadGaussian(\
				sfunc_sig,\
				hits_trk_is_pion,\
				s_func_x,\
				s_func_y,\
				s_func_t,\
				min_pval);
		DircSpreadGaussian * pdf_as_kaon = new DircSpreadGaussian(\
				sfunc_sig,\
				hits_trk_is_kaon,\
				s_func_x,\
				s_func_y,\
				s_func_t,\
				min_pval);//could make this an array that is filled...
		unsigned int r=0;
		while(f>>ievent_index>>iPID>>iBAR>>ix>>iy>>it>>itheta>>iphi>>iE)
		{
			// 			if (iE > 4) continue;
			if (r >= (unsigned int) max_particles) break;

			if (force_kinematics == true)
			{
				itheta = particle_theta;
				iphi = particle_phi;
				iE = energy;
			}

			event_index.push_back(ievent_index);
			PID.push_back(iPID);
			BAR.push_back(iBAR);
			x.push_back(ix);
			y.push_back(iy);
			t.push_back(it);
			theta.push_back(itheta);
			phi.push_back(iphi);
			E.push_back(iE);
			r++;
		}
		f.close();
		for(unsigned int n=0;n<r;n++){
			if(fabs(PID[n])>7 && fabs(PID[n])<13)//8 and 9 are pions, 11 and 12 are kaons
			{

				mc_tally++;

				//
				//and replace acos......
				//cdd replace 1.47 with float refrac_index=1.47;
				pion_mc_beta = dirc_model->get_beta(E[n],pimass);
				kaon_mc_beta = dirc_model->get_beta(E[n],kmass);
				pion_mc_angle = rad_to_deg*acos(1/(1.47*pion_mc_beta));
				kaon_mc_angle = rad_to_deg*acos(1/(1.47*kaon_mc_beta));

				//sep_pdfs_mc = new DircProbabilitySeparation(pdf_as_kaon,pdf_as_pion);

				dirc_model->set_focus_mirror_angle(\
						spread_ang.Gaus(main_mirror_angle,mirror_angle_change_unc),\
						spread_ang.Gaus(0,mirror_angle_change_yunc));
				dirc_model->set_upper_wedge_angle_diff(\
						spread_ang.Gaus(0,wedge_uncertainty),\
						spread_ang.Gaus(0,upper_wedge_yang_spread));
				printf("\r                                                                                                      ");
				printf("\rRunning monte carlo event  %i/%i. Found %i pions/kaons", n,r-1,mc_tally);

				fflush(stdout);

				//simulate pmt hits from primary mc event to compare with preconstructed pdfs
				//


				if(abs(PID[n])==8 || PID[n]==9 ){
					dirc_model->sim_rand_n_photons(\
							sim_points,\
							n_sim_phots,\
							pion_mc_angle,\
							BAR[n],\
							x[n],\
							y[n],\
							0,\
							theta[n],\
							phi[n],\
							0,\
							ckov_unc,\
							pion_mc_beta);
				//TODO: Check timing
					pion_kaon_diff = 1;
				}
				else{
					dirc_model->sim_rand_n_photons(\
							sim_points,\
							n_sim_phots,\
							kaon_mc_angle,\
							BAR[n],\
							x[n],\
							y[n],\
							0,\
							theta[n],\
							phi[n],\
							0,\
							ckov_unc,\
							kaon_mc_beta);
				//TODO: Check timing
					pion_kaon_diff = -1;
				}
				for (unsigned int i = 0; i < sim_points.size(); i++)
				{
					sim_points[i].t += t[n];
					// 					sim_points[i].x += 150;
					// 					printf("%12.04f\n",sim_points[i].t);
				}
				//simulate pmt hits from background events that occur within t ns of the primary event
				confounded_tally=0;
				for (unsigned int j =0; j < r;j++){

					if (abs(event_index[n] - event_index[j]) > broaden_events)
					{
						continue;
					}

					if(fabs(t[n]-t[j])<time_window && j!=n){
						confounded_tally++; 
						if(abs(PID[j])==2 ||PID[j]==3){

							dirc_model->sim_rand_n_photons(\
									confound_points,\
									n_sim_phots,\
									47.135,\
									BAR[j],\
									x[j],\
									y[j],\
									0,\
									theta[j],\
									phi[j],\
									tracking_unc,\
									ckov_unc,\
									1);
						//TODO: Check timing
						}
						else if(abs(PID[j])==8||PID[j]==9){
							//pion_beta = dirc_model->get_beta(E[j],pimass);
							pion_angle = rad_to_deg*acos(1/(refrac_index*pion_beta));


						}
						else if(abs(PID[j])==11 || PID[j]==12){
							kaon_beta = dirc_model->get_beta(E[j],kmass);
							kaon_angle = rad_to_deg*acos(1/(refrac_index*kaon_beta));

						}
						else if(abs(PID[j])==5 || PID[j]==6){
							muon_beta = dirc_model->get_beta(E[j],mumass);
							muon_angle = rad_to_deg*acos(1/(refrac_index*muon_beta));

							dirc_model->sim_rand_n_photons(\
									confound_points,\
									n_sim_phots,\
									muon_angle,\
									BAR[j],\
									x[j],\
									y[j],\
									0,\
									theta[j],\
									phi[j],\
									tracking_unc,\
									ckov_unc,\
									muon_beta);
							//TODO: Check timing
						}
					}
					for (unsigned int i = 0; i < confound_points.size(); i++)
					{
						confound_points[i].t += t[j];
						// 						confound_points[i].x += 150;
						sim_points.push_back(confound_points[i]);
					}
					confound_points.clear();
				}//end confounded point generation

				//Begin analysis

				digitizer.digitize_points(sim_points);



				timing_clock = clock();

				if (use_prog_sep == true)
				{
					dirc_model->set_focus_mirror_angle(\
							spread_ang.Gaus(main_mirror_angle,0),\
							spread_ang.Gaus(0,0));
					dirc_model->set_upper_wedge_angle_diff(\
							spread_ang.Gaus(0,0),\
							spread_ang.Gaus(0,0));
					printf("Progressive separation removed pending optical sim upgrade\n");
				}
				else
				{
					dirc_model->set_focus_mirror_angle(\
							spread_ang.Gaus(main_mirror_angle,0),\
							spread_ang.Gaus(0,0));
					dirc_model->set_upper_wedge_angle_diff(\
							spread_ang.Gaus(0,0),\
							spread_ang.Gaus(0,0));

					dirc_model->sim_reg_n_photons(\
							hits_trk_is_pion,\
							n_phi_phots,\
							n_z_phots,\
							pion_mc_angle,\
							BAR[n],\
							x[n],\
							y[n],\
							0,
							theta[n],\
							phi[n],\
							0,\
							ckov_unc/pdf_unc_red_fac,\
							pion_mc_beta);
						//TODO: Check timing

					dirc_model->sim_reg_n_photons(\
							hits_trk_is_kaon,\
							n_phi_phots,\
							n_z_phots,\
							kaon_mc_angle,
							BAR[n],\
							x[n],\
							y[n],\
							0,
							theta[n],\
							phi[n],\
							0,\
							ckov_unc/pdf_unc_red_fac,\
							kaon_mc_beta);
						//TODO: Check timing

					for (unsigned int k = 0; k < hits_trk_is_pion.size(); k++)
					{
						hits_trk_is_pion[k].t += t[n];
					}
					for (unsigned int k = 0; k < hits_trk_is_kaon.size(); k++)
					{
						hits_trk_is_kaon[k].t += t[n];
					}
					pdf_as_pion->set_support(hits_trk_is_pion);
					pdf_as_kaon->set_support(hits_trk_is_kaon);


					// 					printf("betas: %12.04f %12.04f\n",pion_mc_beta,kaon_mc_beta);

					llc = pdf_as_pion->get_log_likelihood(sim_points);
					llf = pdf_as_kaon->get_log_likelihood(sim_points);

					// 					printf("%06d\n",sim_points.size());

					ll_diff = llc - llf;

					// 					printf("\n%12.04f %12.04f %12.04f\n",llc,llf,ll_diff);
					//int correct_id = 1;
					if (ll_diff*pion_kaon_diff < 0)
					{
						//correct_id = 0;

						printf("wrong: %12.04f %d\n",ll_diff,pion_kaon_diff);
						for (unsigned int k = 0; k < sim_points.size(); k++)
						{
							printf("%12.04f %12.04f %12.04f ll: %12.04f %12.04f, %12.04f\n",sim_points[k].x,sim_points[k].y,sim_points[k].t,\
									log(pdf_as_pion->get_single_likelihood(sim_points[k])),\
									log(pdf_as_kaon->get_single_likelihood(sim_points[k])),\
									log(pdf_as_pion->get_single_likelihood(sim_points[k]))-log(pdf_as_kaon->get_single_likelihood(sim_points[k])));
						}
					}
				}

				timing_clock = clock() - timing_clock;

				simulation_time->Fill(((float)timing_clock)/(CLOCKS_PER_SEC)*1000);


				// 				printf("\nPID[n]=%i loglikehood difference(pion-kaon)=%f\n",PID[n],ll_diff);
				if(abs(PID[n])==8||PID[n]==9){
					ll_diff_pion->Fill(ll_diff);

					if (ll_diff > 0)
					{
						zero_pion_id->Fill(1);
					}
					else
					{
						zero_pion_id->Fill(0);
					}

					phot_found_pion->Fill(sim_points.size());
				}
				else{

					ll_diff_kaon->Fill(ll_diff);
					if (ll_diff < 0)
					{
						zero_kaon_id->Fill(1);
					}
					else
					{
						zero_kaon_id->Fill(0);
					}
					phot_found_kaon->Fill(sim_points.size());
				}
			}
		}
		printf("\n");//ensure we're on a new line at the end
		printf("Kaon missID: %12.04f%%\n",100.0*zero_kaon_id->GetBinContent(1)/(zero_kaon_id->GetBinContent(1)+zero_kaon_id->GetBinContent(2)));
		printf("Pion missID: %12.04f%%\n",100.0*zero_pion_id->GetBinContent(1)/(zero_pion_id->GetBinContent(1)+zero_pion_id->GetBinContent(2)));
		//end mc reading script

	}
	else if (slac_run == true)
	{ 
		printf("Testing the SLAC geometry/run\n");

		pion_beta = dirc_model->get_beta(energy,pimass);
		kaon_beta = dirc_model->get_beta(energy,kmass);

		std::vector<dirc_point> hit_points_pion;
		std::vector<dirc_point> hit_points_kaon;

		dirc_model->set_liquid_absorbtion(liquid_absorbtion);
		dirc_model->set_liquid_index(liquid_index);
		dirc_model->set_three_seg_mirror(three_seg_mirror);
		dirc_model->set_sidemirror(sm_xr,sm_xl);


		dirc_model->set_pmt_offset(pmt_offset);
		dirc_model->set_upper_wedge_angle_diff(wedge_uncertainty);
		dirc_model->set_bar_box_angle(bar_box_box_angle);

		//ns
		float pion_time = particle_flight_distance/(pion_beta*.3);
		float kaon_time = particle_flight_distance/(kaon_beta*.3);

		dirc_model->sim_reg_n_photons(\
				hit_points_pion,\
				n_phi_phots/10,\
				n_z_phots,\
				pion_angle,\
				1,\
				particle_x,\
				particle_y,\
				pion_time,\
				particle_theta,\
				particle_phi,\
				0,\
				ckov_unc/pdf_unc_red_fac,\
				pion_beta);

		dirc_model->sim_reg_n_photons(\
				hit_points_kaon,\
				n_phi_phots/10,\
				n_z_phots,\
				pion_angle,\
				1,\
				particle_x,\
				particle_y,\
				kaon_time,\
				particle_theta,\
				particle_phi,\
				0,\
				ckov_unc/pdf_unc_red_fac,\
				kaon_beta);


		DircSpreadGaussian* pdf_pion = new DircSpreadGaussian(\
				sfunc_sig,\
				hit_points_pion,\
				s_func_x,\
				s_func_y,\
				s_func_t);
		DircSpreadGaussian* pdf_kaon = new DircSpreadGaussian(\
				sfunc_sig,\
				hit_points_kaon,\
				s_func_x,\
				s_func_y,\
				s_func_t);

		for (int i = 0; i < num_runs; i++)
		{
			dirc_model->set_focus_mirror_angle(\
					spread_ang.Gaus(main_mirror_angle,mirror_angle_change_unc),\
					spread_ang.Gaus(0,mirror_angle_change_yunc));
			dirc_model->set_upper_wedge_angle_diff(\
					spread_ang.Gaus(0,wedge_uncertainty),\
					spread_ang.Gaus(0,upper_wedge_yang_spread));
			dirc_model->set_bar_box_angle(spread_ang.Gaus(0,box_rot_unc));

			printf("\r                                                    ");
			printf("\rrunning iter %8d/%d  ",i+1,num_runs);


			fflush(stdout);

			particle_theta = spread_ang.Uniform(0,box_rot_unc);
			particle_phi = spread_ang.Uniform(0,360);
			//Sim new photons
			
			//ns
			float pion_time = particle_flight_distance/(pion_beta*.3);
			float kaon_time = particle_flight_distance/(kaon_beta*.3);

			dirc_model->sim_reg_n_photons(\
					hit_points_pion,\
					n_phi_phots/phi_phots_reduce,\
					n_z_phots,\
					pion_angle,\
					1,\
					particle_x,\
					particle_y,\
					pion_time,\
					particle_theta,\
					particle_phi,\
					0,\
					ckov_unc/pdf_unc_red_fac,\
					pion_beta);

			dirc_model->sim_reg_n_photons(\
					hit_points_kaon,\
					n_phi_phots/phi_phots_reduce,\
					n_z_phots,\
					pion_angle,\
					1,\
					particle_x,\
					particle_y,\
					kaon_time,\
					particle_theta,\
					particle_phi,\
					0,\
					ckov_unc/pdf_unc_red_fac,\
					kaon_beta);

			pdf_pion->set_support(hit_points_pion);
			pdf_kaon->set_support(hit_points_kaon);


			n_sim_phots = spread_ang.Gaus(mean_n_phot,spread_n_phot);

			//assume its a middle bar
			dirc_model->sim_rand_n_photons(\
					sim_points,\
					n_sim_phots,\
					pion_angle,\
					1,\
					particle_x,\
					particle_y,\
					pion_time,\
					particle_theta,\
					particle_phi,\
					tracking_unc,\
					ckov_unc,\
					pion_beta);
			digitizer.digitize_points(sim_points);


			llc = pdf_pion->get_log_likelihood(sim_points);
			llf = pdf_kaon->get_log_likelihood(sim_points);


			ll_diff_pion->Fill(1*(llc-llf));

			phot_found_pion->Fill(sim_points.size());

			dirc_model->sim_rand_n_photons(\
					sim_points,\
					n_sim_phots,\
					kaon_angle,\
					1,\
					particle_x,\
					particle_y,\
					kaon_time,\
					particle_theta,\
					particle_phi,\
					tracking_unc,\
					ckov_unc,\
					kaon_beta);


			digitizer.digitize_points(sim_points);

			llc = pdf_pion->get_log_likelihood(sim_points);
			llf = pdf_kaon->get_log_likelihood(sim_points);

			ll_diff_kaon->Fill(1*(llc-llf));

			// 		ll_diff_kaon->Fill(sep_pdfs->get_log_likelihood_spread_diff(sim_points));

			phot_found_kaon->Fill(sim_points.size());
		}

		printf("\nRun Completed\n");
		if (coverage_plot == true)
		{
			printf("starting coverage ouput\n");

			for (int i = 0; i < num_cov; i++)
			{
				if ((i+1) % 1000 == 0)
				{
					//printf("\r                                                                                                                           ");
					printf("\rcoverage iter %d/%d",i+1,num_cov);
				}
				fflush(stdout);

				particle_theta = spread_ang.Uniform(0,14);
				particle_phi = spread_ang.Uniform(0,360);
				energy = spread_ang.Uniform(2,4.5);
				pion_beta = dirc_model->get_beta(energy,pimass);
				kaon_beta = dirc_model->get_beta(energy,kmass);
				dirc_model->sim_rand_n_photons(\
						sim_points,\
						n_sim_phots,\
						pion_angle,\
						1,\
						particle_x,\
						particle_y,\
						0,\
						particle_theta,\
						particle_phi,\
						tracking_unc,\
						ckov_unc,\
						pion_beta);

				for (unsigned int j = 0; j < sim_points.size(); j++)
				{
					//printf("%12.04f %12.04f\n",sim_points[j].x,sim_points[j].y);
					pion_coverage_xy->Fill(sim_points[j].x+outcsv_x,sim_points[j].y);
				}

				dirc_model->sim_rand_n_photons(\
						sim_points,\
						n_sim_phots,\
						kaon_angle,\
						1,\
						particle_x,\
						particle_y,\
						0,\
						particle_theta,\
						particle_phi,\
						tracking_unc,\
						ckov_unc,\
						kaon_beta);

				for (unsigned int j = 0; j < sim_points.size(); j++)
				{
					kaon_coverage_xy->Fill(sim_points[j].x+outcsv_x,sim_points[j].y);
				}

			}
			printf("\nfinished output of coverage plots\n");
		}
		else
		{
			pion_coverage_xy->Fill(0.0,0.0);
			kaon_coverage_xy->Fill(0.0,0.0);

		}

		//	printf("Found %d pion points on the target\n", (int) hit_points_pion.size());
		float x,y,t_ns;
		for (unsigned int i = 0; i < hit_points_pion.size(); i++)
		{
			x = hit_points_pion[i].x;
			y = hit_points_pion[i].y;
			t_ns = hit_points_pion[i].t;
			// 		if (hit_points_pion[i].t > 0) continue;
			pion_dist_xy->Fill(x,y);
			pion_dist_xt->Fill(x,t_ns);
			pion_dist_yt->Fill(y,t_ns);
			pion_dist_t->Fill(t_ns);
		}

		//	printf("Found %d kaon points on the target\n", (int) hit_points_kaon.size());
		for (unsigned int i = 0; i < hit_points_kaon.size(); i++)
		{
			x = hit_points_kaon[i].x;
			y = hit_points_kaon[i].y;
			t_ns = hit_points_kaon[i].t;
			// 		if (hit_points_pion[i].t > 0) continue;
			kaon_dist_xy->Fill(x,y);
			kaon_dist_xt->Fill(x,t_ns);
			kaon_dist_yt->Fill(y,t_ns);
			kaon_dist_t->Fill(t_ns);
		}
	}
	else if (sparse_recon_n > 0)
	{ 
		printf("Trying sparse reconstruction loop mode\n");
		
		if (flatten_time==true || sep_updown == true || monochrome_plot == true)
		{
			printf("flatten_time, sep_updown, and monochrome_plot  not currently implemented for the sparse reconstruction.\n");
		}

		int num_provisional_points = 25000;
		float points_dist_sq = 100;
		float rethrown_phi_width = .1;
		int num_rethrown_phots = 1200000 * rethrown_phi_width/(3.14159);

		pion_beta = dirc_model->get_beta(energy,pimass);
		kaon_beta = dirc_model->get_beta(energy,kmass);

		std::vector<dirc_point> hit_points_pion;
		std::vector<dirc_point> hit_points_kaon;

		dirc_model->set_liquid_absorbtion(liquid_absorbtion);
		dirc_model->set_liquid_index(liquid_index);
		dirc_model->set_three_seg_mirror(three_seg_mirror);
		dirc_model->set_sidemirror(sm_xr,sm_xl);


		dirc_model->set_pmt_offset(pmt_offset);
		dirc_model->set_upper_wedge_angle_diff(wedge_uncertainty);
		dirc_model->set_bar_box_angle(bar_box_box_angle);

		//ns
		float pion_time = particle_flight_distance/(pion_beta*.3);
		float kaon_time = particle_flight_distance/(kaon_beta*.3);


		DircSpreadGaussian* pdf_pion = new DircSpreadGaussian(\
				sfunc_sig,\
				hit_points_pion,\
				s_func_x,\
				s_func_y,\
				s_func_t);
		DircSpreadGaussian* pdf_kaon = new DircSpreadGaussian(\
				sfunc_sig,\
				hit_points_kaon,\
				s_func_x,\
				s_func_y,\
				s_func_t);

		for (int i = 0; i < sparse_recon_n; i++)
		{
			dirc_model->set_focus_mirror_angle(\
					spread_ang.Gaus(main_mirror_angle,mirror_angle_change_unc),\
					spread_ang.Gaus(0,mirror_angle_change_yunc));
			dirc_model->set_upper_wedge_angle_diff(\
					spread_ang.Gaus(0,wedge_uncertainty),\
					spread_ang.Gaus(0,upper_wedge_yang_spread));
			dirc_model->set_bar_box_angle(spread_ang.Gaus(0,box_rot_unc));

			if (particle_theta_mean < .01)
			{
				particle_phi = spread_ang.Uniform(0,360);
			}

			particle_theta = spread_ang.Gaus(particle_theta_mean, particle_theta_spread);
			particle_x = spread_ang.Gaus(particle_x_mean, particle_x_spread);
			particle_y = spread_ang.Gaus(particle_y_mean, particle_y_spread);
			energy = spread_ang.Gaus(energy_mean,energy_spread);
			pion_beta = dirc_model->get_beta(energy,pimass);
			kaon_beta = dirc_model->get_beta(energy,kmass);

			printf("\r                                                    ");
			printf("\rrunning iter %8d/%d  ",i+1,sparse_recon_n);


			fflush(stdout);

			n_sim_phots = spread_ang.Gaus(mean_n_phot,spread_n_phot);

			//assume its a middle bar
			dirc_model->sim_rand_n_photons(\
					sim_points,\
					n_sim_phots,\
					pion_angle,\
					1,\
					particle_x,\
					particle_y,\
					pion_time,\
					particle_theta+const_track_off,\
					particle_phi,\
					tracking_unc,\
					ckov_unc,\
					pion_beta);
			digitizer.digitize_points(sim_points);

			llc = 0;
			llf = 0;
			
			std::vector<float> start_phi;
			std::vector<dirc_point> provisional_points;


			for (int i = 0; i < num_provisional_points; i++)
			{
				dirc_point tmp_pro_point;
				float cur_phi = i*2*3.14159/num_provisional_points;
				bool hit_pmt;
				hit_pmt = dirc_model->track_single_photon_beta(\
				        tmp_pro_point,\
				        pion_beta,\
       	 				cur_phi,\
       	 				particle_theta,\
        				particle_phi,\
        				particle_x,\
        				particle_y,\
        				-17/2,\
        				0,\
        				1);
					//last 3 are z, t, and bar;

				if (hit_pmt == true)
				{
					//I could have done these as a pair, but I'll just keep them synchronized;
					start_phi.push_back(cur_phi);
					provisional_points.push_back(tmp_pro_point);
				}
			}

			for (unsigned int i = 0; i < sim_points.size(); i++)
			{
				float cur_mean_phi = 0;
				float cur_mean_x = 0;
				float cur_mean_y = 0;
				
				int total_prov_points = 0;
				float cur_dist_sq = 0;
				for (unsigned int j = 0; j < provisional_points.size(); j++)
				{
					cur_dist_sq = 0;
					cur_dist_sq += (sim_points[i].x-provisional_points[j].x)*(sim_points[i].x-provisional_points[j].x);
					cur_dist_sq += (sim_points[i].y-provisional_points[j].y)*(sim_points[i].y-provisional_points[j].y);
					cur_dist_sq += (sim_points[i].t-provisional_points[j].t)*(sim_points[i].t-provisional_points[j].t);
					
					if (cur_dist_sq < points_dist_sq)
					{
						cur_mean_phi += start_phi[j];
						cur_mean_x += cos(start_phi[j]);
						cur_mean_y += sin(start_phi[j]);
						total_prov_points++;
					}
				}
				//cur_mean_phi /= total_prov_points;
				//printf("pion/pion: %12.04f %12.04f\n",cur_mean_phi/total_prov_points,atan2(cur_mean_y,cur_mean_x));
				cur_mean_phi = atan2(cur_mean_y,cur_mean_x);
				hit_points_pion.clear();
				for (int j = 0; j < num_rethrown_phots; j++)
				{
					float cur_phi = spread_ang.Uniform(cur_mean_phi-rethrown_phi_width,cur_mean_phi+rethrown_phi_width);
					float cur_z = spread_ang.Uniform(-17,0);
					dirc_point tmp_point;
					float hit_pmt = dirc_model->track_single_photon_beta(\
					        tmp_point,\
					        pion_beta,\
       	 					cur_phi,\
       	 					particle_theta,\
        					particle_phi,\
        					particle_x,\
        					particle_y,\
        					cur_z,\
        					0,\
        					1);
					if (hit_pmt == true)
					{
						hit_points_pion.push_back(tmp_point);
					}
				}	
				pdf_pion->set_support(hit_points_pion);
				if (hit_points_pion.size() > 0)
				{
					std::vector<dirc_point> single_point;
					single_point.push_back(sim_points[i]);
					float add_ll = pdf_pion->get_log_likelihood_new_support(single_point,hit_points_pion);
//					printf("pion/pion: %d/%d: %12.04f\n",i,sim_points.size(),add_ll);
					llc += add_ll;
					//printf("pion llc: %12.04f\n",llc);
					if (add_ll != add_ll)
					{
						for (unsigned int j = 0; j < hit_points_pion.size(); j++)
						{
					//		printf("pion/pion: %12.04f %12.04f %12.04f : %12.04f %12.04f %12.04f\n",sim_points[i].x,sim_points[i].y,sim_points[i].t,hit_points_pion[j].x,hit_points_pion[j].y,hit_points_pion[j].t);
						}
					}
				}
	
			}
			start_phi.clear();
			provisional_points.clear();	
			for (int i = 0; i < num_provisional_points; i++)
			{
				dirc_point tmp_pro_point;
				float cur_phi = i*2*3.14159/num_provisional_points;
				bool hit_pmt = false;
				hit_pmt = dirc_model->track_single_photon_beta(\
				        tmp_pro_point,\
				        kaon_beta,\
       	 				cur_phi,\
       	 				particle_theta,\
        				particle_phi,\
        				particle_x,\
        				particle_y,\
        				-17/2,\
        				0,\
        				1);
					//last 3 are z, t, and bar;

				if (hit_pmt == true)
				{
					//I could have done these as a pair, but I'll just keep them synchronized;
					start_phi.push_back(cur_phi);
					provisional_points.push_back(tmp_pro_point);
				}
			}

			for (unsigned int i = 0; i < sim_points.size(); i++)
			{
				float cur_mean_phi = 0;
				float cur_mean_x = 0;
				float cur_mean_y = 0;
				int total_prov_points = 0;
				float cur_dist_sq = 0;
				for (unsigned int j = 0; j < provisional_points.size(); j++)
				{
					cur_dist_sq = 0;
					cur_dist_sq += (sim_points[i].x-provisional_points[j].x)*(sim_points[i].x-provisional_points[j].x);
					cur_dist_sq += (sim_points[i].y-provisional_points[j].y)*(sim_points[i].y-provisional_points[j].y);
					cur_dist_sq += (sim_points[i].t-provisional_points[j].t)*(sim_points[i].t-provisional_points[j].t);
					
					if (cur_dist_sq < points_dist_sq)
					{
						cur_mean_phi += start_phi[j];
						cur_mean_x += cos(start_phi[j]);
						cur_mean_y += sin(start_phi[j]);
						total_prov_points++;
					}
				}
				//cur_mean_phi /= total_prov_points;
				//printf("pion/kaon: %12.04f %12.04f\n",cur_mean_phi/total_prov_points,atan2(cur_mean_y,cur_mean_x));
				cur_mean_phi = atan2(cur_mean_y,cur_mean_x);
				hit_points_kaon.clear();
				for (int j = 0; j < num_rethrown_phots; j++)
				{
					float cur_phi = spread_ang.Uniform(cur_mean_phi-rethrown_phi_width,cur_mean_phi+rethrown_phi_width);
					float cur_z = spread_ang.Uniform(-17,0);
					dirc_point tmp_point;
					float hit_pmt = dirc_model->track_single_photon_beta(\
					        tmp_point,\
					        kaon_beta,\
       	 					cur_phi,\
       	 					particle_theta,\
        					particle_phi,\
        					particle_x,\
        					particle_y,\
        					cur_z,\
        					0,\
        					1);
					if (hit_pmt == true)
					{
						hit_points_kaon.push_back(tmp_point);
					}
				}	
				pdf_kaon->set_support(hit_points_kaon);
				if (hit_points_kaon.size() > 0)
				{
					std::vector<dirc_point> single_point;
					single_point.push_back(sim_points[i]);
					float add_ll = pdf_kaon->get_log_likelihood_new_support(single_point,hit_points_kaon);
//					printf("pion/kaon: %d/%d: %12.04f\n",i,sim_points.size(),add_ll);
					llf += add_ll;
					//printf("kaon llf: %12.04f\n",llf);
					if (add_ll != add_ll)
					{
						for (unsigned int j = 0; j < hit_points_kaon.size(); j++)
						{
					//		printf("pion/kaon: %12.04f %12.04f %12.04f : %12.04f %12.04f %12.04f\n",sim_points[i].x,sim_points[i].y,sim_points[i].t,hit_points_kaon[j].x,hit_points_kaon[j].y,hit_points_kaon[j].t);
						}
					}
				}
	
			}


//			llc = pdf_pion->get_log_likelihood(sim_points);
//			llf = pdf_kaon->get_log_likelihood(sim_points);

			ll_diff_pion->Fill(1*(llc-llf));
			phot_found_pion->Fill(sim_points.size());
			printf("\nll_diff pion: %12.04f\n",llc-llf);


			dirc_model->sim_rand_n_photons(\
					sim_points,\
					n_sim_phots,\
					kaon_angle,\
					1,\
					particle_x,\
					particle_y,\
					kaon_time,\
					particle_theta+const_track_off,\
					particle_phi,\
					tracking_unc,\
					ckov_unc,\
					kaon_beta);

			digitizer.digitize_points(sim_points);

			start_phi.clear();
			provisional_points.clear();	
			llc=0;
			llf=0;
			for (int i = 0; i < num_provisional_points; i++)
			{
				dirc_point tmp_pro_point;
				float cur_phi = i*2*3.14159/num_provisional_points;
				bool hit_pmt;
				hit_pmt = dirc_model->track_single_photon_beta(\
				        tmp_pro_point,\
				        pion_beta,\
       	 				cur_phi,\
       	 				particle_theta,\
        				particle_phi,\
        				particle_x,\
        				particle_y,\
        				-17/2,\
        				0,\
        				1);
					//last 3 are z, t, and bar;

				if (hit_pmt == true)
				{
					//I could have done these as a pair, but I'll just keep them synchronized;
					start_phi.push_back(cur_phi);
					provisional_points.push_back(tmp_pro_point);
				}
			}

			for (unsigned int i = 0; i < sim_points.size(); i++)
			{
				float cur_mean_phi = 0;
				float cur_mean_x = 0;
				float cur_mean_y = 0;
				int total_prov_points = 0;
				float cur_dist_sq = 0;
				for (unsigned int j = 0; j < provisional_points.size(); j++)
				{
					cur_dist_sq = 0;
					cur_dist_sq += (sim_points[i].x-provisional_points[j].x)*(sim_points[i].x-provisional_points[j].x);
					cur_dist_sq += (sim_points[i].y-provisional_points[j].y)*(sim_points[i].y-provisional_points[j].y);
					cur_dist_sq += (sim_points[i].t-provisional_points[j].t)*(sim_points[i].t-provisional_points[j].t);
					
					if (cur_dist_sq < points_dist_sq)
					{
						cur_mean_phi += start_phi[j];
						cur_mean_x += cos(start_phi[j]);
						cur_mean_y += sin(start_phi[j]);
						total_prov_points++;
					}
				}
				//cur_mean_phi /= total_prov_points;
				cur_mean_phi = atan2(cur_mean_y,cur_mean_x);
				hit_points_pion.clear();
				for (int j = 0; j < num_rethrown_phots; j++)
				{
					float cur_phi = spread_ang.Uniform(cur_mean_phi-rethrown_phi_width,cur_mean_phi+rethrown_phi_width);
					float cur_z = spread_ang.Uniform(-17,0);
					dirc_point tmp_point;
					float hit_pmt = dirc_model->track_single_photon_beta(\
					        tmp_point,\
					        pion_beta,\
       	 					cur_phi,\
       	 					particle_theta,\
        					particle_phi,\
        					particle_x,\
        					particle_y,\
        					cur_z,\
        					pion_time,\
        					1);
					if (hit_pmt == true)
					{
						hit_points_pion.push_back(tmp_point);
					}
				}	
				pdf_pion->set_support(hit_points_pion);
				if (hit_points_pion.size() > 0)
				{
					std::vector<dirc_point> single_point;
					single_point.push_back(sim_points[i]);
					float add_ll = pdf_pion->get_log_likelihood_new_support(single_point,hit_points_pion);
				//	printf("pion/kaon: %12.04f\n",add_ll);
					llc += add_ll;
				}
	
			}
			start_phi.clear();
			provisional_points.clear();	
			for (int i = 0; i < num_provisional_points; i++)
			{
				dirc_point tmp_pro_point;
				float cur_phi = i*2*3.14159/num_provisional_points;
				bool hit_pmt;
				hit_pmt = dirc_model->track_single_photon_beta(\
				        tmp_pro_point,\
				        kaon_beta,\
       	 				cur_phi,\
       	 				particle_theta,\
        				particle_phi,\
        				particle_x,\
        				particle_y,\
        				-17.25/2,\
        				kaon_time,\
        				1);
					//last 3 are z, t, and bar;

				if (hit_pmt == true)
				{
					//I could have done these as a pair, but I'll just keep them synchronized;
					start_phi.push_back(cur_phi);
					provisional_points.push_back(tmp_pro_point);
				}
			}

			for (unsigned int i = 0; i < sim_points.size(); i++)
			{
				float cur_mean_phi = 0;
				float cur_mean_x = 0;
				float cur_mean_y = 0;
				int total_prov_points = 0;
				float cur_dist_sq = 0;
				for (unsigned int j = 0; j < provisional_points.size(); j++)
				{
					cur_dist_sq = 0;
					cur_dist_sq += (sim_points[i].x-provisional_points[j].x)*(sim_points[i].x-provisional_points[j].x);
					cur_dist_sq += (sim_points[i].y-provisional_points[j].y)*(sim_points[i].y-provisional_points[j].y);
					cur_dist_sq += (sim_points[i].t-provisional_points[j].t)*(sim_points[i].t-provisional_points[j].t);
					
					if (cur_dist_sq < points_dist_sq)
					{
						cur_mean_phi += start_phi[j];
						cur_mean_x += cos(start_phi[j]);
						cur_mean_y += sin(start_phi[j]);
						total_prov_points++;
					}
				}
				//cur_mean_phi /= total_prov_points;
				cur_mean_phi = atan2(cur_mean_y,cur_mean_x);
				hit_points_kaon.clear();
				for (int j = 0; j < num_rethrown_phots; j++)
				{
					float cur_phi = spread_ang.Uniform(cur_mean_phi-rethrown_phi_width,cur_mean_phi+rethrown_phi_width);
					float cur_z = spread_ang.Uniform(-17,0);
					dirc_point tmp_point;
					float hit_pmt = dirc_model->track_single_photon_beta(\
					        tmp_point,\
					        kaon_beta,\
       	 					cur_phi,\
       	 					particle_theta,\
        					particle_phi,\
        					particle_x,\
        					particle_y,\
        					cur_z,\
        					kaon_time,\
        					1);
					if (hit_pmt == true)
					{
						hit_points_kaon.push_back(tmp_point);
					}
				}	
				pdf_kaon->set_support(hit_points_kaon);
				if (hit_points_kaon.size() > 0)
				{
					std::vector<dirc_point> single_point;
					single_point.push_back(sim_points[i]);
					float add_ll = pdf_kaon->get_log_likelihood(single_point);
				//	printf("kaon/kaon: %12.04f\n",add_ll);
					llf += add_ll;
				}
	
			}
//			llc = pdf_pion->get_log_likelihood(sim_points);
//			llf = pdf_kaon->get_log_likelihood(sim_points);

			ll_diff_kaon->Fill(1*(llc-llf));
			phot_found_kaon->Fill(sim_points.size());
			printf("ll_diff kaon: %12.04f\n",llc-llf);

		}

		printf("\nSparse Run Completed\n");
	}
	else if (run_minuit_calibration == true)
	{
		dirc_model->set_liquid_absorbtion(liquid_absorbtion);
		dirc_model->set_liquid_index(liquid_index);
		dirc_model->set_three_seg_mirror(three_seg_mirror);
		dirc_model->set_sidemirror(sm_xr,sm_xl);
	
	
		dirc_model->set_pmt_offset(pmt_offset);
		dirc_model->set_upper_wedge_angle_diff(wedge_uncertainty);
		dirc_model->set_bar_box_angle(bar_box_box_angle);
		
		global_calib_model = dynamic_cast<DircBaseSim*>(dirc_model);
		global_digitizer = &digitizer; //check this address math
		global_particle_theta = particle_theta;
		global_particle_phi = particle_phi;
		global_particle_x = particle_x;
		global_particle_y = particle_y;
		global_particle_t = 0; //implement flight distance later
		global_particle_energy = energy;
		global_particle_flight_distance = particle_flight_distance;
		global_refrac_index = refrac_index;
		global_n_sim_phots = n_sim_phots;
		global_ckov_unc = ckov_unc;
		global_tracking_unc = tracking_unc;


		init_global_data();

		int n_minuit_pars = 2;	
		float *minuit_pars = new float[n_minuit_pars];
/*
		float *minuit_grad;//not_computed
		float minuit_val;
		float minuit_flag = 0;
*/
		minuit_pars[0] = pion_x_adj;
		minuit_pars[1] = pion_y_adj;
		minuit_pars[2] = kaon_x_adj;
//		minuit_pars[3] = kaon_y_adj;

		//dirc_calib_line(n_minuit_pars, minuit_grad, minuit_val, minuit_pars, minuit_flag);
		//printf("Minuit result: %12.04f\n",minuit_val);

		TMinuit minuit(n_minuit_pars);
		minuit.SetFCN(dirc_calib_line);
		// Vector of step, initial min and max value
		float arglist[10];
		float vstrt[n_minuit_pars];
		float stp[n_minuit_pars];
 		float bmin[n_minuit_pars];
 		float bmax[n_minuit_pars];

		vstrt[0] = 0; 
		vstrt[1] = 0;
//		vstrt[2] = 0;
//		vstrt[3] = 0;
//		vstrt[4] = 50;
 		stp[0] = .5; 
		stp[1] = .5;
//		stp[2] = 1;
//		stp[3] = 1;
//		stp[4] = 10;
		bmin[0] = -20; 
		bmin[1] = bmin[0];
//		bmin[2] = bmin[0];
//		bmin[3] = bmin[0];
//		bmin[4] = 5;
 		bmax[0] = 20; 
		bmax[1] = bmax[0];
//		bmax[2] = bmax[0];
//		bmax[3] = bmax[0];
//		bmax[4] = 200;

		char param_title[64];
		int ierflg = 0;
		for (int i = 0; i < n_minuit_pars; i++)
		{
			sprintf(param_title,"minuit_param_%03d",i);
			minuit.mnparm(i, param_title, vstrt[i],stp[i],bmin[i],bmax[i],ierflg);
		}

		arglist[0] = 100;
		arglist[1] = .01;
		minuit.SetPrintLevel(-1);	
		minuit.mnexcm("SIM", arglist ,1,ierflg);
//		minuit.mnexcm("IMP", arglist ,1,ierflg);
		
		float x_shift = -.6;
		float y_shift, y_split;
		float cur_err;
		float mult_params = 1;

		minuit.GetParameter(0,y_shift,cur_err);
		minuit.GetParameter(1,y_split,cur_err);

		pion_x_adj = x_shift;
		kaon_x_adj = x_shift;
		pion_y_adj = mult_params*(y_shift+y_split);
		kaon_y_adj = mult_params*(y_shift-y_split);
		printf("minuit_calib_string:\n");
		printf("%12.04f %12.04f %12.04f %12.04f\n",pion_x_adj,pion_y_adj,kaon_x_adj,kaon_y_adj);

	}
	else if (line_recon_n > 0)
	{ 
		printf("Trying line reconstruction loop mode\n");
		printf("Caliibrations: pion_x_adj pion_y_adj kaon_x_adj kaon_y_adj\n");
		printf("%12.04f %12.04f %12.04f %12.04f\n",pion_x_adj, pion_y_adj, kaon_x_adj, kaon_y_adj);
		
		if (flatten_time==true || sep_updown == true || monochrome_plot == true)
		{
			printf("flatten_time, sep_updown, and monochrome_plot  not currently implemented for the line reconstruction.\n");
		}

		int num_line_points = 5000;
//		float points_dist_sq = 6000;
		float time_spread = 2;
		//float time_spread = 1000;
		//float dist_spread = 40;
		float dist_spread = 50;
		float max_dev_sq = 5;


		max_dev_sq *= max_dev_sq;

		pion_beta = dirc_model->get_beta(energy,pimass);
		kaon_beta = dirc_model->get_beta(energy,kmass);

		std::vector<dirc_point> hit_points_pion;
		std::vector<dirc_point> hit_points_kaon;

		dirc_model->set_liquid_absorbtion(liquid_absorbtion);
		dirc_model->set_liquid_index(liquid_index);
		dirc_model->set_three_seg_mirror(three_seg_mirror);
		dirc_model->set_sidemirror(sm_xr,sm_xl);


		dirc_model->set_pmt_offset(pmt_offset);
		dirc_model->set_upper_wedge_angle_diff(wedge_uncertainty);
		dirc_model->set_bar_box_angle(bar_box_box_angle);

		//ns
		float pion_time = particle_flight_distance/(pion_beta*.3);
		float kaon_time = particle_flight_distance/(kaon_beta*.3);
		float unused_photons = 0;

		float pion_t_adj = .34;
		pion_t_adj = 0;

		for (int i = 0; i < line_recon_n; i++)
		{
			dirc_model->set_focus_mirror_angle(\
					spread_ang.Gaus(main_mirror_angle,mirror_angle_change_unc),\
					spread_ang.Gaus(0,mirror_angle_change_yunc));
			dirc_model->set_upper_wedge_angle_diff(\
					spread_ang.Gaus(0,wedge_uncertainty),\
					spread_ang.Gaus(0,upper_wedge_yang_spread));
			dirc_model->set_bar_box_angle(spread_ang.Gaus(0,box_rot_unc));

			particle_theta = spread_ang.Gaus(particle_theta_mean, particle_theta_spread);
			particle_x = spread_ang.Gaus(particle_x_mean, particle_x_spread);
			particle_y = spread_ang.Gaus(particle_y_mean, particle_y_spread);
			energy = spread_ang.Gaus(energy_mean,energy_spread);
			pion_beta = dirc_model->get_beta(energy,pimass);
			kaon_beta = dirc_model->get_beta(energy,kmass);

			printf("\r                                                    ");
			printf("\rrunning iter %8d/%d  ",i+1,line_recon_n);


			fflush(stdout);

			n_sim_phots = spread_ang.Gaus(mean_n_phot,spread_n_phot);

			//assume its a middle bar
			dirc_model->sim_rand_n_photons(\
				sim_points,\
				n_sim_phots,\
				pion_angle,\
				1,\
				particle_x,\
				particle_y,\
				pion_time,\
				particle_theta+const_track_off,\
				particle_phi,\
				tracking_unc,\
				ckov_unc,\
				pion_beta);
			digitizer.digitize_points(sim_points);

			llc = 0;
			llf = 0;
			
			std::vector<float> phi_last_wall_neg;
			std::vector<float> phi_last_wall_pos;
			std::vector<dirc_point> provisional_points_lwn;
			std::vector<dirc_point> provisional_points_lwp;

			float refrac_split = -.000;
			float pion_refrac_mod = 1.000;
			float kaon_refrac_mod = 1.000;

			pion_refrac_mod += refrac_split;
			kaon_refrac_mod -= refrac_split;

		//	pion_refrac_mod = 1;
		//	kaon_refrac_mod = 1;


			float pion_emit_angle = 57.3*acos(1/(pion_beta*pion_refrac_mod*refrac_index));
			float kaon_emit_angle = 57.3*acos(1/(kaon_beta*kaon_refrac_mod*refrac_index));

			dirc_model->track_all_line_photons(\
			        provisional_points_lwn,\
			        provisional_points_lwp,\
				num_line_points,\
			        pion_emit_angle,\
       				particle_theta,\
       				particle_phi,\
       				particle_x,\
       				particle_y,\
       				-17.25/2,\
       				pion_time,\
       				1);

			unused_photons = 0;

			for (unsigned int i = 0; i < provisional_points_lwn.size(); i++)
			{
				provisional_points_lwn[i].x += pion_x_adj;
				provisional_points_lwn[i].y += pion_y_adj;
				provisional_points_lwn[i].t += pion_t_adj;
			}
			for (unsigned int i = 0; i < provisional_points_lwp.size(); i++)
			{
				provisional_points_lwp[i].x += pion_x_adj;
				provisional_points_lwp[i].y += pion_y_adj;
				provisional_points_lwp[i].t += pion_t_adj;
			}


			for (unsigned int i = 0; i < sim_points.size(); i++)
			{
				float min_dist_sq = 10000;
				float cur_dist_sq = -1;
				for (unsigned int j = 0; j < provisional_points_lwn.size(); j++)
				{
					cur_dist_sq = (sim_points[i].x - provisional_points_lwn[j].x)*(sim_points[i].x - provisional_points_lwn[j].x);		
					cur_dist_sq += (sim_points[i].y - provisional_points_lwn[j].y)*(sim_points[i].y - provisional_points_lwn[j].y);		
					cur_dist_sq /= dist_spread;
					cur_dist_sq += (sim_points[i].t - provisional_points_lwn[j].t)*(sim_points[i].t - provisional_points_lwn[j].t)/time_spread;	
					min_dist_sq = std::min(min_dist_sq,cur_dist_sq);
				}
				for (unsigned int j = 0; j < provisional_points_lwp.size(); j++)
				{
					cur_dist_sq = (sim_points[i].x - provisional_points_lwp[j].x)*(sim_points[i].x - provisional_points_lwp[j].x);		
					cur_dist_sq += (sim_points[i].y - provisional_points_lwp[j].y)*(sim_points[i].y - provisional_points_lwp[j].y);		
					cur_dist_sq /= dist_spread;
					cur_dist_sq += (sim_points[i].t - provisional_points_lwp[j].t)*(sim_points[i].t - provisional_points_lwp[j].t)/time_spread;	
					min_dist_sq = std::min(min_dist_sq,cur_dist_sq);
				}
				//cur_mean_phi /= total_prov_points;
				//printf("pion/pion: %12.04f %12.04f\n",cur_mean_phi/total_prov_points,atan2(cur_mean_y,cur_mean_x));
				
				if (min_dist_sq > max_dev_sq)
				{
					unused_photons++;
				}
				llc += std::min(min_dist_sq,max_dev_sq);
			//	printf("%4d %12.04f %12.04f\n",i,min_dist_sq,llc);

			}
			pion_unused_photons->Fill(unused_photons);
			phi_last_wall_neg.clear();
			phi_last_wall_pos.clear();
			provisional_points_lwn.clear();	
			provisional_points_lwp.clear();	

			dirc_model->track_all_line_photons(\
			        provisional_points_lwn,\
			        provisional_points_lwp,\
				num_line_points,\
			        kaon_emit_angle,\
       				particle_theta,\
       				particle_phi,\
       				particle_x,\
       				particle_y,\
       				-17.25/2,\
       				kaon_time,\
       				1);

			for (unsigned int i = 0; i < provisional_points_lwn.size(); i++)
			{
				provisional_points_lwn[i].x += kaon_x_adj;
				provisional_points_lwn[i].y += kaon_y_adj;
			}
			for (unsigned int i = 0; i < provisional_points_lwp.size(); i++)
			{
				provisional_points_lwp[i].x += kaon_x_adj;
				provisional_points_lwp[i].y += kaon_y_adj;
			}
			for (unsigned int i = 0; i < sim_points.size(); i++)
			{
				float min_dist_sq = 10000;
				float cur_dist_sq = -1;
				for (unsigned int j = 0; j < provisional_points_lwn.size(); j++)
				{
					cur_dist_sq = (sim_points[i].x - provisional_points_lwn[j].x)*(sim_points[i].x - provisional_points_lwn[j].x);		
					cur_dist_sq += (sim_points[i].y - provisional_points_lwn[j].y)*(sim_points[i].y - provisional_points_lwn[j].y);		
					cur_dist_sq /= dist_spread;
					cur_dist_sq += (sim_points[i].t - provisional_points_lwn[j].t)*(sim_points[i].t - provisional_points_lwn[j].t)/time_spread;	
					min_dist_sq = std::min(min_dist_sq,cur_dist_sq);
				}
				for (unsigned int j = 0; j < provisional_points_lwp.size(); j++)
				{
					cur_dist_sq = (sim_points[i].x - provisional_points_lwp[j].x)*(sim_points[i].x - provisional_points_lwp[j].x);		
					cur_dist_sq += (sim_points[i].y - provisional_points_lwp[j].y)*(sim_points[i].y - provisional_points_lwp[j].y);		
					cur_dist_sq /= dist_spread;
					cur_dist_sq += (sim_points[i].t - provisional_points_lwp[j].t)*(sim_points[i].t - provisional_points_lwp[j].t)/time_spread;	
					min_dist_sq = std::min(min_dist_sq,cur_dist_sq);
				}
				//cur_mean_phi /= total_prov_points;
				//printf("pion/pion: %12.04f %12.04f\n",cur_mean_phi/total_prov_points,atan2(cur_mean_y,cur_mean_x));
				llf += std::min(min_dist_sq,max_dev_sq);

			}

			//Flip to make sign agree with ll
			ll_diff_pion->Fill(-1*(llc-llf));
			phot_found_pion->Fill(sim_points.size());
			//printf("\nll_diff pion: %12.04f\n",llc-llf);


			dirc_model->sim_rand_n_photons(\
					sim_points,\
					n_sim_phots,\
					kaon_angle,\
					1,\
					particle_x,\
					particle_y,\
					kaon_time,\
					particle_theta+const_track_off,\
					particle_phi,\
					tracking_unc,\
					ckov_unc,\
					kaon_beta);
			//printf("model run: %12.04f %12.04f %12.04f %12.04f %12.04f %12.04f %12.04f %12.04f %12.04f\n",kaon_angle,particle_x,particle_y,kaon_time,particle_theta+const_track_off,particle_phi,tracking_unc,ckov_unc,kaon_beta);


			digitizer.digitize_points(sim_points);

			phi_last_wall_neg.clear();
			phi_last_wall_pos.clear();
			provisional_points_lwn.clear();	
			provisional_points_lwp.clear();	
			llc=0;
			llf=0;

			dirc_model->track_all_line_photons(\
			        provisional_points_lwn,\
			        provisional_points_lwp,\
				num_line_points,\
			        pion_emit_angle,\
       				particle_theta,\
       				particle_phi,\
       				particle_x,\
       				particle_y,\
       				-17.25/2,\
       				pion_time,\
       				1);

			for (unsigned int i = 0; i < provisional_points_lwn.size(); i++)
			{
				provisional_points_lwn[i].x += pion_x_adj;
				provisional_points_lwn[i].y += pion_y_adj;
				provisional_points_lwn[i].t += pion_t_adj;
			}
			for (unsigned int i = 0; i < provisional_points_lwp.size(); i++)
			{
				provisional_points_lwp[i].x += pion_x_adj;
				provisional_points_lwp[i].y += pion_y_adj;
				provisional_points_lwp[i].t += pion_t_adj;
			}
			for (unsigned int i = 0; i < sim_points.size(); i++)
			{
				float min_dist_sq = 10000;
				float cur_dist_sq = -1;
				for (unsigned int j = 0; j < provisional_points_lwn.size(); j++)
				{
					cur_dist_sq = (sim_points[i].x - provisional_points_lwn[j].x)*(sim_points[i].x - provisional_points_lwn[j].x);		
					cur_dist_sq += (sim_points[i].y - provisional_points_lwn[j].y)*(sim_points[i].y - provisional_points_lwn[j].y);		
					cur_dist_sq /= dist_spread;
					cur_dist_sq += (sim_points[i].t - provisional_points_lwn[j].t)*(sim_points[i].t - provisional_points_lwn[j].t)/time_spread;	
					min_dist_sq = std::min(min_dist_sq,cur_dist_sq);
				}
				for (unsigned int j = 0; j < provisional_points_lwp.size(); j++)
				{
					cur_dist_sq = (sim_points[i].x - provisional_points_lwp[j].x)*(sim_points[i].x - provisional_points_lwp[j].x);		
					cur_dist_sq += (sim_points[i].y - provisional_points_lwp[j].y)*(sim_points[i].y - provisional_points_lwp[j].y);		
					cur_dist_sq /= dist_spread;
					cur_dist_sq += (sim_points[i].t - provisional_points_lwp[j].t)*(sim_points[i].t - provisional_points_lwp[j].t)/time_spread;	
					min_dist_sq = std::min(min_dist_sq,cur_dist_sq);
				}
				//cur_mean_phi /= total_prov_points;
				//printf("pion/pion: %12.04f %12.04f\n",cur_mean_phi/total_prov_points,atan2(cur_mean_y,cur_mean_x));
				llc += std::min(min_dist_sq,max_dev_sq);

			}
			phi_last_wall_neg.clear();
			phi_last_wall_pos.clear();
			provisional_points_lwn.clear();	
			provisional_points_lwp.clear();	

			dirc_model->track_all_line_photons(\
			        provisional_points_lwn,\
			        provisional_points_lwp,\
				num_line_points,\
			        kaon_emit_angle,\
       				particle_theta,\
       				particle_phi,\
       				particle_x,\
       				particle_y,\
       				-17.25/2,\
       				kaon_time,\
       				1);
			unused_photons = 0;
			for (unsigned int i = 0; i < provisional_points_lwn.size(); i++)
			{
				provisional_points_lwn[i].x += kaon_x_adj;
				provisional_points_lwn[i].y += kaon_y_adj;
			}
			for (unsigned int i = 0; i < provisional_points_lwp.size(); i++)
			{
				provisional_points_lwp[i].x += kaon_x_adj;
				provisional_points_lwp[i].y += kaon_y_adj;
			}
			for (unsigned int i = 0; i < sim_points.size(); i++)
			{
				float min_dist_sq = 10000;
				float cur_dist_sq = -1;
				for (unsigned int j = 0; j < provisional_points_lwn.size(); j++)
				{
					cur_dist_sq = (sim_points[i].x - provisional_points_lwn[j].x)*(sim_points[i].x - provisional_points_lwn[j].x);		
					cur_dist_sq += (sim_points[i].y - provisional_points_lwn[j].y)*(sim_points[i].y - provisional_points_lwn[j].y);		
					cur_dist_sq /= dist_spread;
					cur_dist_sq += (sim_points[i].t - provisional_points_lwn[j].t)*(sim_points[i].t - provisional_points_lwn[j].t)/time_spread;	
					min_dist_sq = std::min(min_dist_sq,cur_dist_sq);
				}
				for (unsigned int j = 0; j < provisional_points_lwp.size(); j++)
				{
					cur_dist_sq = (sim_points[i].x - provisional_points_lwp[j].x)*(sim_points[i].x - provisional_points_lwp[j].x);		
					cur_dist_sq += (sim_points[i].y - provisional_points_lwp[j].y)*(sim_points[i].y - provisional_points_lwp[j].y);		
					cur_dist_sq /= dist_spread;
					cur_dist_sq += (sim_points[i].t - provisional_points_lwp[j].t)*(sim_points[i].t - provisional_points_lwp[j].t)/time_spread;	
					min_dist_sq = std::min(min_dist_sq,cur_dist_sq);
				}
				//cur_mean_phi /= total_prov_points;
				//printf("pion/pion: %12.04f %12.04f\n",cur_mean_phi/total_prov_points,atan2(cur_mean_y,cur_mean_x));
				if (min_dist_sq > max_dev_sq)
				{
					unused_photons++;
				}
				llf += std::min(min_dist_sq,max_dev_sq);

			}
			kaon_unused_photons->Fill(unused_photons);

			//Flip to make sign agree with ll
			ll_diff_kaon->Fill(-1*(llc-llf));
			phot_found_kaon->Fill(sim_points.size());
			//printf("ll_diff kaon: %12.04f\n",llc-llf);

		}
		printf("\nLine Recon Run Completed\n");
		//printf("%12.04f\n",ival);
	}
	else if (num_runs > 0)
	{
		 
		printf("no input file specified.  Running in loop mode\n");


		pion_beta = dirc_model->get_beta(energy,pimass);
		kaon_beta = dirc_model->get_beta(energy,kmass);

		if (monochrome_plot == true)
		{
			pion_angle = rad_to_deg*acos(1/(refrac_index*pion_beta));
			kaon_angle = rad_to_deg*acos(1/(refrac_index*kaon_beta));
			pion_beta = -1;
			kaon_beta = -1;
			ckov_unc = 0;
		}

		std::vector<dirc_point> hit_points_pion;
		std::vector<dirc_point> hit_points_kaon;



		// 	dirc_model->set_upper_wedge_angle_diff(0);

		dirc_model->set_use_moliere(use_moliere_scattering);
		dirc_model->set_moliere_p(energy*1000);//assume momentum is the same for both for now - high energy;

		dirc_model->set_liquid_absorbtion(liquid_absorbtion);
		dirc_model->set_liquid_index(liquid_index);
		dirc_model->set_three_seg_mirror(three_seg_mirror);
		dirc_model->set_sidemirror(sm_xr,sm_xl);


		dirc_model->set_pmt_offset(pmt_offset);
		dirc_model->set_upper_wedge_angle_diff(wedge_uncertainty);
		dirc_model->set_bar_box_angle(bar_box_box_angle);

		//ns
		float pion_time = particle_flight_distance/(pion_beta*.3);
		float kaon_time = particle_flight_distance/(kaon_beta*.3);

		dirc_model->sim_reg_n_photons(\
				hit_points_pion,\
				n_phi_phots,\
				n_z_phots,\
				-1,\
				1,\
				particle_x,\
				particle_y,\
				pion_time,\
				particle_theta,\
				particle_phi,\
				0,\
				ckov_unc/pdf_unc_red_fac,\
				pion_beta);
//TODO RETURN

		dirc_model->sim_reg_n_photons(\
				hit_points_kaon,\
				n_phi_phots,\
				n_z_phots,\
				-1,\
				1,\
				particle_x,\
				particle_y,\
				kaon_time,\
				particle_theta,\
				particle_phi,\
				0,\
				ckov_unc/pdf_unc_red_fac,\
				kaon_beta);

/*
		dirc_model->sim_rand_n_photons(\
				hit_points_pion,\
				n_phi_phots*n_z_phots,\
				-1,\
				1,\
				particle_x,\
				particle_y,\
				pion_time,\
				particle_theta,\
				particle_phi,\
				0,\
				ckov_unc/pdf_unc_red_fac,\
				pion_beta);

		dirc_model->sim_reg_n_photons(\
				hit_points_kaon,\
				n_phi_phots*n_z_phots,\
				-1,\
				1,\
				particle_x,\
				particle_y,\
				kaon_time,\
				particle_theta,\
				particle_phi,\
				0,\
				ckov_unc/pdf_unc_red_fac,\
				kaon_beta);
*/
		if (flatten_time == true)
		{
			for (unsigned int i = 0; i < hit_points_pion.size(); i++)
			{
				hit_points_pion[i].t = 0;
			}
			for (unsigned int i = 0; i < hit_points_kaon.size(); i++)
			{
				hit_points_kaon[i].t = 0;
			}
		}


		DircSpreadGaussian* pdf_pion = new DircSpreadGaussian(\
				sfunc_sig,\
				hit_points_pion,\
				s_func_x,\
				s_func_y,\
				s_func_t);
		DircSpreadGaussian* pdf_kaon = new DircSpreadGaussian(\
				sfunc_sig,\
				hit_points_kaon,\
				s_func_x,\
				s_func_y,\
				s_func_t);

		for (int i = 0; i < num_runs; i++)
		{
			dirc_model->set_focus_mirror_angle(\
					spread_ang.Gaus(main_mirror_angle,mirror_angle_change_unc)+main_mirror_angle_off,\
					spread_ang.Gaus(0,mirror_angle_change_yunc)+main_mirror_yangle_off,\
					main_mirror_zangle_off);
			dirc_model->set_mirror_plane_offsets(\
					main_mirror_yoff,\
					main_mirror_zoff);
			dirc_model->set_bar_box_offsets(\
					bar_box_xoff,\
					bar_box_yoff,\
					bar_box_zoff);


			dirc_model->set_upper_wedge_angle_diff(\
					spread_ang.Gaus(0,wedge_uncertainty),\
					spread_ang.Gaus(0,upper_wedge_yang_spread));
			dirc_model->set_bar_box_angle(spread_ang.Gaus(0,box_rot_unc));

			if (particle_theta_mean < .01)
			{
				particle_phi = spread_ang.Uniform(0,360);
			}

			particle_theta = spread_ang.Gaus(particle_theta_mean, particle_theta_spread);
			particle_x = spread_ang.Gaus(particle_x_mean, particle_x_spread);
			particle_y = spread_ang.Gaus(particle_y_mean, particle_y_spread);
			energy = spread_ang.Gaus(energy_mean,energy_spread);
			pion_beta = dirc_model->get_beta(energy,pimass);
			kaon_beta = dirc_model->get_beta(energy,kmass);

			if (monochrome_plot == true)
			{
				pion_angle = rad_to_deg*acos(1/(refrac_index*pion_beta));
				kaon_angle = rad_to_deg*acos(1/(refrac_index*kaon_beta));
				pion_beta = -1;
				kaon_beta = -1;
				ckov_unc = 0;
			}
			//			printf("%12.04f %12.04f %12.04f %12.04f\n",particle_x, particle_y, particle_theta, energy);

			printf("\r                                                    ");
			printf("\rrunning iter %8d/%d  ",i+1,num_runs);


			fflush(stdout);

			n_sim_phots = spread_ang.Gaus(mean_n_phot,spread_n_phot);

			//assume its a middle bar
			dirc_model->sim_rand_n_photons(\
					sim_points,\
					n_sim_phots,\
					pion_angle,\
					1,\
					particle_x,\
					particle_y,\
					pion_time,\
					particle_theta+const_track_off,\
					particle_phi,\
					tracking_unc,\
					ckov_unc,\
					pion_beta);
			digitizer.digitize_points(sim_points);

			if (flatten_time == true)
			{
				for (unsigned int i = 0; i < sim_points.size(); i++)
				{
					sim_points[i].t = 0;
				}
			}
			llc = pdf_pion->get_log_likelihood(sim_points);
			llf = pdf_kaon->get_log_likelihood(sim_points);
			
			ll_diff_pion->Fill(1*(llc-llf));
			phot_found_pion->Fill(sim_points.size());

			if (sep_updown == true)
			{
				std::vector<dirc_point> sim_points_up;
				std::vector<dirc_point> sim_points_down;
				for (unsigned int j = 0; j < sim_points.size(); j++)
				{
					if (sim_points[j].updown == 1)
					{
						sim_points_up.push_back(sim_points[j]);
					}
					else if (sim_points[j].updown == -1)
					{
						sim_points_down.push_back(sim_points[j]);
					}
					else
					{
						printf("Found a sim point without an updown flag.  This shouldn't happen.\n");
					}
				}
				//just with the up
				llc = pdf_pion->get_log_likelihood(sim_points_up);
				llf = pdf_kaon->get_log_likelihood(sim_points_up);


				ll_diff_pion_up->Fill(1*(llc-llf));

				//just with the down
				llc = pdf_pion->get_log_likelihood(sim_points_down);
				llf = pdf_kaon->get_log_likelihood(sim_points_down);

				ll_diff_pion_down->Fill(1*(llc-llf));

				phot_found_pion_up->Fill(sim_points_up.size());
				phot_found_pion_down->Fill(sim_points_down.size());
			}




			dirc_model->sim_rand_n_photons(\
					sim_points,\
					n_sim_phots,\
					kaon_angle,\
					1,\
					particle_x,\
					particle_y,\
					kaon_time,\
					particle_theta+const_track_off,\
					particle_phi,\
					tracking_unc,\
					ckov_unc,\
					kaon_beta);

			if (flatten_time == true)
			{
				for (unsigned int i = 0; i < sim_points.size(); i++)
				{
					sim_points[i].t = 0;
				}
			}

			digitizer.digitize_points(sim_points);

			llc = pdf_pion->get_log_likelihood(sim_points);
			llf = pdf_kaon->get_log_likelihood(sim_points);

			ll_diff_kaon->Fill(1*(llc-llf));
			phot_found_kaon->Fill(sim_points.size());

			if (sep_updown == true)
			{
				std::vector<dirc_point> sim_points_up;
				std::vector<dirc_point> sim_points_down;
				for (unsigned int j = 0; j < sim_points.size(); j++)
				{
					if (sim_points[j].updown == 1)
					{
						sim_points_up.push_back(sim_points[j]);
					}
					else if (sim_points[j].updown == -1)
					{
						sim_points_down.push_back(sim_points[j]);
					}
					else
					{
						printf("Found a sim point without an updown flag.  This shouldn't happen.\n");
					}
				}
				//just with the up
				llc = pdf_pion->get_log_likelihood(sim_points_up);
				llf = pdf_kaon->get_log_likelihood(sim_points_up);

				ll_diff_kaon_up->Fill(1*(llc-llf));

				//just with the down
				llc = pdf_pion->get_log_likelihood(sim_points_down);
				llf = pdf_kaon->get_log_likelihood(sim_points_down);

				ll_diff_kaon_down->Fill(1*(llc-llf));
				phot_found_kaon_up->Fill(sim_points_up.size());
				phot_found_kaon_down->Fill(sim_points_down.size());
			}
			// 		ll_diff_kaon->Fill(sep_pdfs->get_log_likelihood_spread_diff(sim_points));

		}

		printf("\nRun Completed\n");
		if (coverage_plot == true)
		{
			printf("starting coverage ouput\n");

			for (int i = 0; i < num_cov; i++)
			{
				if ((i+1) % 1000 == 0)
				{
					//printf("\r                                                                                                                           ");
					printf("\rcoverage iter %d/%d",i+1,num_cov);
				}
				fflush(stdout);

				particle_theta = spread_ang.Uniform(0,14);
				particle_phi = spread_ang.Uniform(0,360);
				energy = spread_ang.Uniform(2,4.5);
				pion_beta = dirc_model->get_beta(energy,pimass);
				kaon_beta = dirc_model->get_beta(energy,kmass);
				dirc_model->sim_rand_n_photons(\
						sim_points,\
						n_sim_phots,\
						pion_angle,\
						1,\
						particle_x,\
						particle_y,\
						0,\
						particle_theta,\
						particle_phi,\
						tracking_unc,\
						ckov_unc,\
						pion_beta);

				for (unsigned int j = 0; j < sim_points.size(); j++)
				{
					//printf("%12.04f %12.04f\n",sim_points[j].x,sim_points[j].y);
					pion_coverage_xy->Fill(sim_points[j].x+outcsv_x,sim_points[j].y);
				}

				dirc_model->sim_rand_n_photons(\
						sim_points,\
						n_sim_phots,\
						kaon_angle,\
						1,\
						particle_x,\
						particle_y,\
						0,\
						particle_theta,\
						particle_phi,\
						tracking_unc,\
						ckov_unc,\
						kaon_beta);

				for (unsigned int j = 0; j < sim_points.size(); j++)
				{
					kaon_coverage_xy->Fill(sim_points[j].x+outcsv_x,sim_points[j].y);
				}

			}
			printf("\nfinished output of coverage plots\n");
		}
		else
		{
			pion_coverage_xy->Fill(0.0,0.0);
			kaon_coverage_xy->Fill(0.0,0.0);

		}
		if (fill_kinematics_yields == true)
		{
			int itheta_bin = 0;
			int iphi_bin = 0;
			printf("starting kinematics yields ouput\n");
			for (float itheta = kinematics_yields_min_theta + kinematics_yields_step_theta/2; itheta < kinematics_yields_max_theta; itheta += kinematics_yields_step_theta)
			{
				itheta_bin++;//I know there's a better way to do this, but I'm lazy
				iphi_bin = 0;
				for (float iphi = kinematics_yields_min_phi + kinematics_yields_step_phi/2; iphi < kinematics_yields_max_phi; iphi += kinematics_yields_step_phi)
				{
					iphi_bin++;
					printf("\r                                                                             ");
					printf("\rPhi: %12.04f    Theta: %12.04f",iphi,itheta);
					fflush(stdout);	
					pion_beta = dirc_model->get_beta(5,pimass); //High energy
					dirc_model->sim_rand_n_photons(
							sim_points,
							kinematics_n_phots,
							pion_angle,
							1,
							0, //particle_x
							0, //particle_y
							0,\
							itheta,
							iphi,
							0, //ignore tracking
							ckov_unc,
							pion_beta);

					kinematics_yields->SetBinContent(itheta_bin,iphi_bin,(1.0*sim_points.size())/kinematics_n_phots);
				}
			}
			printf("\n");
		}


		//Make the distribution - use random	
		pion_beta = dirc_model->get_beta(energy,pimass);
		kaon_beta = dirc_model->get_beta(energy,kmass);
		if (monochrome_plot == true)
		{
			pion_angle = rad_to_deg*acos(1/(refrac_index*pion_beta));
			kaon_angle = rad_to_deg*acos(1/(refrac_index*kaon_beta));
			pion_beta = -1;
			kaon_beta = -1;
			ckov_unc = 0;
		}
		//	printf("Found %d pion points on the target\n", (int) hit_points_pion.size());
		float x,y,t_ns;
		for (unsigned int i = 0; i < hit_points_pion.size(); i++)
		{
			x = hit_points_pion[i].x;
			y = hit_points_pion[i].y;
			t_ns = hit_points_pion[i].t;
			// 		if (hit_points_pion[i].t > 0) continue;
			pion_dist_x->Fill(x);
			pion_dist_y->Fill(y);
			pion_dist_t->Fill(t_ns);
			pion_dist_xy->Fill(x,y);
			pion_dist_xy_geantbins->Fill(x,y+50);
			//pion_dist_xy_geantbins->Fill(x,y);
			pion_dist_xy_geantbins_cm->Fill(x/10,(y+50)/10);
			pion_dist_xt->Fill(x,t_ns);
			pion_dist_yt->Fill(y,t_ns);
			pion_dist_t->Fill(t_ns);
			//printf("pion_dist_xyt: %12.04f %12.04f %12.04f\n",x,y,t_ns);
		}

		//	printf("Found %d kaon points on the target\n", (int) hit_points_kaon.size());
		for (unsigned int i = 0; i < hit_points_kaon.size(); i++)
		{
			x = hit_points_kaon[i].x;
			y = hit_points_kaon[i].y;
			t_ns = hit_points_kaon[i].t;
			// 		if (hit_points_pion[i].t > 0) continue;
			kaon_dist_x->Fill(x);
			kaon_dist_y->Fill(y);
			kaon_dist_t->Fill(t_ns);
			kaon_dist_xy->Fill(x,y);
			kaon_dist_xt->Fill(x,t_ns);
			kaon_dist_yt->Fill(y,t_ns);
			kaon_dist_t->Fill(t_ns);
		}
	}
	float dist_mean = 0;
	if (box_check_n > 0)
	{	
		printf("Box Checking Histograms...\n");
		std::vector<dirc_point> box_check_points;

		dirc_model->test_from_wedge_top(\
				box_check_points,\
				box_check_n,\
				1 /*Bar*/,\
				particle_x,\
				box_check_theta,\
				box_check_phi,\
				box_check_theta_unc,\
				box_check_phi_unc,\
				box_check_overall_theta);

		dist_mean = 0;

		for (unsigned int i = 0; i < box_check_points.size(); i++)
		{
			box_check_x->Fill(box_check_points[i].x);	
			box_check_y->Fill(box_check_points[i].y);	
			box_check_d->Fill(box_check_points[i].t);
			dist_mean += box_check_points[i].t;
			box_check_xy->Fill(box_check_points[i].x,box_check_points[i].y);
		}
		dist_mean /= box_check_points.size();

		printf("Box Check Results\n");
		printf("n-phots XSpread YSpread Dist X/Dist Y/Dist\n");
		printf("%07d %12.04f %12.04f %12.04f %12.04e %12.04e\n",(int) box_check_points.size(),box_check_x->GetRMS(),box_check_y->GetRMS(),dist_mean,box_check_x->GetRMS()/dist_mean,box_check_y->GetRMS()/dist_mean);


	}
	if (phot_check_n > 0)
	{	
		printf("Checking Number of Photons Found\n");
		std::vector<dirc_point> phot_check_points;

		int num_phot_points = 100;
		float phot_check_theta = 0;
		float phot_check_phi = 0;
		dist_mean = 0;	

		for (int i = 0; i < phot_check_n; i++)
		{
			pion_beta = dirc_model->get_beta(energy,pimass);
			kaon_beta = dirc_model->get_beta(energy,kmass);
			if (monochrome_plot == true)
			{
				pion_angle = rad_to_deg*acos(1/(refrac_index*pion_beta));
				kaon_angle = rad_to_deg*acos(1/(refrac_index*kaon_beta));
				pion_beta = -1;
				kaon_beta = -1;
				ckov_unc = 0;
			}
			phot_check_theta = spread_ang.Uniform(0,phot_check_max_theta);
			phot_check_phi = spread_ang.Uniform(0,360);
			dirc_model->sim_rand_n_photons(\
					phot_check_points,\
					num_phot_points,\
					pion_angle,\
					1,\
					particle_x,\
					particle_y,\
					0,\
					phot_check_theta,\
					phot_check_phi,\
					tracking_unc,\
					ckov_unc,\
					pion_beta);

			dist_mean += phot_check_points.size();
			//printf("%d\n",phot_check_points.size());
			phot_check_points.clear();
		}

		printf("Relative photons found:\n");
		printf("%12.04f\n",dist_mean/(num_phot_points*phot_check_n));


	}
	if (gaus_ll_n > 0)
	{
		printf("Running Filling LL hists with points drawn from gaussian, Do not perform any other fills\n");
		pion_beta = dirc_model->get_beta(energy,pimass);
		kaon_beta = dirc_model->get_beta(energy,kmass);
	
		float quartz_index = 1.473;
		float pion_cerenkov = acos(1/(quartz_index*pion_beta));
		float kaon_cerenkov = acos(1/(quartz_index*kaon_beta));
		float gaus_mean_a = 0;
		float gaus_mean_b = 1000*(pion_cerenkov-kaon_cerenkov);

		printf("Filling Histos with %d points\n",gaus_ll_n);
		printf("Mean separation:  %12.04f \"mrad\"\n",gaus_mean_b);
		printf("Gaussian Sigma:   %12.04f \"mrad\"\n",gaus_ll_spread);

		float filla,fillb;
		for (int i = 0; i < gaus_ll_n; i++)
		{
			filla = spread_ang.Gaus(gaus_mean_a,gaus_ll_spread);
			fillb = spread_ang.Gaus(gaus_mean_b,gaus_ll_spread);

			//Keeps ordering sane;
			ll_diff_kaon->Fill(filla);
			ll_diff_pion->Fill(fillb);
		}
		printf("Gaussian fill to loglikelihood histograms completed.\n");
	}
	if (lut_sim_n > 0)
	{
		int lut_size = 5000000;
		printf("Performing LUT analysis with: %d points\n",lut_size);

		DircLUTEnum* lut_enum = new DircGluexLUTEnum(
			minx,\
	                maxx,\
	                resx,\
	                digit_miny,\
	                digit_maxy,\
	                resy);
		DircLUT* dirc_lut = new DircLUT(lut_enum);

		pion_beta = dirc_model->get_beta(energy,pimass);
		kaon_beta = dirc_model->get_beta(energy,kmass);
	
		float quartz_index = 1.473;
		float pion_cerenkov = acos(1/(quartz_index*pion_beta));
		float kaon_cerenkov = acos(1/(quartz_index*kaon_beta));

		std::vector<dirc_point> hit_points_pion;
		std::vector<dirc_point> hit_points_kaon;
		std::vector<float> pion_ckov;
		std::vector<float> kaon_ckov;
		std::vector<float> pion_dts;
		std::vector<float> kaon_dts;

		std::vector<float> pion_angle_means;
		std::vector<float> kaon_angle_means;

		std::vector<dirc_point> lut_points;
		std::vector<float> lut_phis;
		std::vector<float> lut_thetas;
		dirc_model->sim_lut_points(\
	                lut_points,\
       	        	lut_phis,\
                	lut_thetas,\
                	lut_size, \
                	1); //bar

		for (unsigned int i = 0; i < lut_points.size(); i++)
		{
			dirc_lut->add_table_pt(lut_points[i],lut_phis[i],lut_thetas[i]);
		}
		printf("LUT Table filled with %d points\n", (int) lut_points.size());
		float fmin = 57.3*(pion_cerenkov+kaon_cerenkov)/2 - 1;
		float fmax = 57.3*(pion_cerenkov+kaon_cerenkov)/2 + 1;

		float oval_cut_angle_center = (pion_cerenkov+kaon_cerenkov)/2;
		float oval_cut_angle_spread_sq = .018;//radians
		//float oval_cut_angle_spread_sq = .1;//radians
		oval_cut_angle_spread_sq *= oval_cut_angle_spread_sq;
		float oval_cut_time_spread_sq = .5;//ns/m - scaled by pathlength
		//float oval_cut_time_spread_sq = 4;//ns/m - scaled by pathlength
		oval_cut_time_spread_sq *= oval_cut_time_spread_sq;//ns
		//float oval_cut_val = -1;

		float ll_mean_adjust = 47.1;
		float ll_mean_enhance = 17.45;

		//all in radians
		float iter_ang_spread_sq = oval_cut_angle_spread_sq/2;
		float iter_time_spread_sq = oval_cut_time_spread_sq/2;
		int max_lut_iter = 20;
		float iter_mean_change_stop = .00002;
		float iter_last_mean = -10;
		float iter_mean = oval_cut_angle_center;


		//estimate chromatic correction, use a pion at our kinematics and repost vals
		float m_cc_direct = 0;
		float b_cc_direct = 0;
		float m_cc_indirect = 0;
		float b_cc_indirect = 0;
		int n_cc_particles = 10000;

		if (perform_chromatic_correction == true)
		{
			printf("Simulating for chromatic corrections\n");

			dirc_model->sim_rand_n_photons(\
	                       	hit_points_pion,\
        	                n_sim_phots*n_cc_particles,\
                	        0, //ckov theta - doesn't matter
                        	1, //bar
                        	particle_x,\
                        	particle_y,\
                        	0, //t
                        	particle_theta,\
                        	particle_phi,\
                        	tracking_unc,\
                        	ckov_unc, //ckov_unc
                        	pion_beta);

			//Do an interative mean finding for ceter of oval?
			dirc_lut->get_chromcorr_m_b_single_oval_cut(\
       				hit_points_pion, \
       				particle_phi, \
       				particle_theta, \
       				particle_y, \
       				pion_cerenkov, \
       				oval_cut_angle_spread_sq,\
       				oval_cut_time_spread_sq,\
				pion_cerenkov,\
				m_cc_direct,\
				b_cc_direct,\
				m_cc_indirect,\
				b_cc_indirect);
		
			printf("Direct m and b: %12.04f %12.04f\n",m_cc_direct,b_cc_direct);
			printf("Indirect m and b: %12.04f %12.04f\n",m_cc_indirect,b_cc_indirect);
		}


		for (int i = 0; i < lut_sim_n; i++)
		{
			if (lut_slac == true)
			{
				//SLAC acceptance
				float x_bound = tan(17/57.3);
				float y_bound = tan(27/57.3);
				float tmp_x = spread_ang.Uniform(-x_bound,x_bound);
				float tmp_y = spread_ang.Uniform(-y_bound,y_bound);

				particle_theta = 57.3*atan(sqrt(tmp_x*tmp_x+tmp_y*tmp_y));
				particle_phi = 57.3*atan2(tmp_x,tmp_y); //This is the correct otdering - I want 0 to be along the bar.
			}

			int n_sim_phots = spread_ang.Gaus(mean_n_phot,spread_n_phot);
			dirc_model->sim_rand_n_photons(\
                        	hit_points_pion,\
                                n_sim_phots,\
                                0, //ckov theta - doesn't matter
                                1, //bar
                                particle_x,\
                                particle_y,\
                                0, //t
                                particle_theta,\
                                particle_phi,\
                                tracking_unc,\
                                ckov_unc, //ckov_unc
                                pion_beta);
			
			dirc_model->sim_rand_n_photons(\
                        	hit_points_kaon,\
                                n_sim_phots,\
                                0, //ckov theta - doesn't matter
                                1, //bar
                                particle_x,\
                                particle_y,\
                                0, //t
                                particle_theta,\
                                particle_phi,\
                                tracking_unc,\
                                ckov_unc, //ckov_unc
                                kaon_beta);

			//pass by reference?  Probably pretty small and doesn't matter
			dirc_lut->get_ckov_theta_all(pion_ckov,pion_dts,hit_points_pion, particle_phi, particle_theta, particle_y);
			dirc_lut->get_ckov_theta_all(kaon_ckov,kaon_dts,hit_points_kaon, particle_phi, particle_theta, particle_y);
			float pion_lut_mean = 0;
			float pion_lut_count = 0;
			float kaon_lut_mean = 0;
			float kaon_lut_count = 0;

			iter_last_mean = -10;
			iter_mean = oval_cut_angle_center;
			for (int j = 0; j < max_lut_iter; j++)
			{
				if (fabs(iter_mean - iter_last_mean) < iter_mean_change_stop) break; //converged
				pion_lut_mean = 0;
				pion_lut_count = 0;
				for (unsigned int k = 0; k < pion_ckov.size(); k++)
				{
					if ((pion_ckov[k]-iter_mean)*(pion_ckov[k]-iter_mean)/iter_ang_spread_sq + pion_dts[k]*pion_dts[k]/iter_time_spread_sq > 1)
					{
						continue;
					}
					pion_lut_mean += pion_ckov[k];
					pion_lut_count += 1.0;
				}
				iter_last_mean = iter_mean;
				iter_mean = pion_lut_mean/pion_lut_count;
				//printf("%3d %12.04f %12.04f\n",j,iter_mean, pion_lut_count);	
			}
			//printf("Pion Peak: %12.04f\n",57.3*iter_mean);
			//rerun with oval cut averaging per point
			/*
			dirc_lut->get_ckov_theta_single_oval_cut(pion_ckov, \
			       	pion_dts, \
        			hit_points_pion, \
        			particle_phi, \
        			particle_theta, \
        			particle_y, \
        			iter_mean, \
        			oval_cut_angle_spread_sq,\
        			oval_cut_time_spread_sq);
			*/
				
			dirc_lut->get_ckov_theta_single_oval_cut(pion_ckov, \
			       	pion_dts, \
        			hit_points_pion, \
        			particle_phi, \
        			particle_theta, \
        			particle_y, \
        			iter_mean, \
        			oval_cut_angle_spread_sq,\
        			oval_cut_time_spread_sq,\
				m_cc_direct,\
				b_cc_direct,\
				m_cc_indirect,\
				b_cc_indirect);

			pion_lut_mean = 0;
			pion_lut_count = 0;

			for (unsigned int j = 0; j < pion_ckov.size(); j++)
			{
				if ((pion_ckov[j]-iter_mean)*(pion_ckov[j]-iter_mean)/oval_cut_angle_spread_sq + pion_dts[j]*pion_dts[j]/oval_cut_time_spread_sq > 1)
				{
					continue;
				}
				pion_lut_vals->Fill(57.3*pion_ckov[j]);
				pion_lut_mean += 57.3*pion_ckov[j];
				pion_lut_count += 1.0;
				pion_lut_dt_v_dang->Fill(1000*(pion_ckov[j]-pion_cerenkov),pion_dts[j]);
			}
			iter_last_mean = -10;
			iter_mean = oval_cut_angle_center;
			//Found mean, rerun to get one per initial point

			for (int j = 0; j < max_lut_iter; j++)
			{
				if (fabs(iter_mean - iter_last_mean) < iter_mean_change_stop) break; //converged
				kaon_lut_mean = 0;
				kaon_lut_count = 0;
				for (unsigned int k = 0; k < kaon_ckov.size(); k++)
				{
					if ((kaon_ckov[k]-iter_mean)*(kaon_ckov[k]-iter_mean)/iter_ang_spread_sq + kaon_dts[k]*kaon_dts[k]/iter_time_spread_sq > 1)
					{
						continue;
					}
					kaon_lut_mean += kaon_ckov[k];
					kaon_lut_count += 1.0;
				}
				iter_last_mean = iter_mean;
				iter_mean = kaon_lut_mean/kaon_lut_count;
				//printf("%3d %12.04f %12.04f\n",j,iter_mean, kaon_lut_count);	
			}
			for (int j = 0; j < max_lut_iter; j++)
			{
				if (fabs(iter_mean - iter_last_mean) < iter_mean_change_stop) break; //converged
				kaon_lut_mean = 0;
				kaon_lut_count = 0;
				for (unsigned int k = 0; k < kaon_ckov.size(); k++)
				{
					if ((kaon_ckov[k]-iter_mean)*(kaon_ckov[k]-iter_mean)/iter_ang_spread_sq + kaon_dts[k]*kaon_dts[k]/iter_time_spread_sq > 1)
					{
						continue;
					}
					kaon_lut_mean += kaon_ckov[k];
					kaon_lut_count += 1.0;
				}
				iter_last_mean = iter_mean;
				iter_mean = kaon_lut_mean/kaon_lut_count;
				//printf("%3d %12.04f %12.04f\n",j,iter_mean, kaon_lut_count);	
			}
			//rerun with oval cut averaging per point
			/*
			dirc_lut->get_ckov_theta_single_oval_cut(kaon_ckov, \
			       	kaon_dts, \
        			hit_points_kaon, \
        			particle_phi, \
        			particle_theta, \
        			particle_y, \
        			iter_mean, \
        			oval_cut_angle_spread_sq,\
        			oval_cut_time_spread_sq);
			*/
			dirc_lut->get_ckov_theta_single_oval_cut(kaon_ckov, \
			       	kaon_dts, \
        			hit_points_kaon, \
        			particle_phi, \
        			particle_theta, \
        			particle_y, \
        			iter_mean, \
        			oval_cut_angle_spread_sq,\
        			oval_cut_time_spread_sq,\
				m_cc_direct,\
				b_cc_direct,\
				m_cc_indirect,\
				b_cc_indirect);
			kaon_lut_mean = 0;
			kaon_lut_count = 0;
			for (unsigned int j = 0; j < kaon_ckov.size(); j++)
			{
				if ((kaon_ckov[j]-iter_mean)*(kaon_ckov[j]-iter_mean)/oval_cut_angle_spread_sq + kaon_dts[j]*kaon_dts[j]/oval_cut_time_spread_sq > 1)
				{
					continue;
				}
				kaon_lut_vals->Fill(57.3*kaon_ckov[j]);
				kaon_lut_mean += 57.3*kaon_ckov[j];
				kaon_lut_count += 1.0;
				kaon_lut_dt_v_dang->Fill(1000*(kaon_ckov[j]-kaon_cerenkov),kaon_dts[j]);
			}


			//printf("%12.04f\n",pion_lut_mean/pion_lut_count);
			//printf("%12.04f\n", pion_ckov[j]);
			pion_lut_means->Fill(pion_lut_mean/pion_lut_count);
			kaon_lut_means->Fill(kaon_lut_mean/kaon_lut_count);
			pion_lut_angles->Fill(pion_lut_count);
			kaon_lut_angles->Fill(kaon_lut_count);

			pion_angle_means.push_back(pion_lut_mean/pion_lut_count);
			kaon_angle_means.push_back(kaon_lut_mean/kaon_lut_count);

			ll_diff_pion->Fill(ll_mean_enhance*(pion_lut_mean/pion_lut_count-ll_mean_adjust));
			ll_diff_kaon->Fill(ll_mean_enhance*(kaon_lut_mean/kaon_lut_count-ll_mean_adjust));

		}
		//ROOT Made me do this, I am so sorry.

		TH1F* tmp_pion_lut = new TH1F(*pion_lut_vals);
		TH1F* tmp_kaon_lut = new TH1F(*kaon_lut_vals);
		
		float pion_mean = 0;//Degrees
		float pion_sigma = 0;//mrad
		float kaon_mean = 0;//Degrees
		float kaon_sigma = 0;//mrad

		
		tmp_pion_lut->Fit("gaus","RQL","",fmin,fmax);
		pion_mean = tmp_pion_lut->GetFunction("gaus")->GetParameter(1);
		pion_sigma = 3141.59/180*tmp_pion_lut->GetFunction("gaus")->GetParameter(2);
		tmp_kaon_lut->Fit("gaus","RQL","",fmin,fmax);
		kaon_mean = tmp_kaon_lut->GetFunction("gaus")->GetParameter(1);
		kaon_sigma = 3141.59/180*tmp_kaon_lut->GetFunction("gaus")->GetParameter(2);
		printf("Per photon resolutions of LUT\n");
		
		printf("Pion - Mean: %12.04f deg         Sigma: %12.04f mrad         %12.04f deg\n",pion_mean,pion_sigma,.057*pion_sigma);
		printf("Kaon - Mean: %12.04f deg         Sigma: %12.04f mrad         %12.04f deg\n",kaon_mean,kaon_sigma,.057*kaon_sigma);
		
		pion_mean = 0;
		pion_sigma = 0;//mrad
		kaon_mean = 0;
		kaon_sigma = 0;//mrad
		for (unsigned int i = 0; i < pion_angle_means.size(); i++)
		{
			pion_mean += pion_angle_means[i];
			kaon_mean += kaon_angle_means[i];
		}
		pion_mean /= pion_angle_means.size();
		kaon_mean /= kaon_angle_means.size();
		for (unsigned int i = 0; i < pion_angle_means.size(); i++)
		{
			pion_sigma += (pion_angle_means[i]-pion_mean)*(pion_angle_means[i]-pion_mean);
			kaon_sigma += (kaon_angle_means[i]-kaon_mean)*(kaon_angle_means[i]-kaon_mean);
		}
		pion_sigma = sqrt(pion_sigma/pion_angle_means.size());
		kaon_sigma = sqrt(kaon_sigma/kaon_angle_means.size());

		//pion_mean *= 17.45;//to mrad
		//kaon_mean *= 17.45;//to mrad
		pion_sigma *= 17.45;//to mrad
		kaon_sigma *= 17.45;//to mrad
		/*	
		tmp_pion_means->Fit("gaus","RQL","",mean_fit_min,mean_fit_max);
		pion_mean = tmp_pion_means->GetFunction("gaus")->GetParameter(1);
		pion_sigma = 3141.59/180*tmp_pion_means->GetFunction("gaus")->GetParameter(2);
		tmp_kaon_means->Fit("gaus","RQL","",mean_fit_min,mean_fit_max);
		kaon_mean = tmp_kaon_means->GetFunction("gaus")->GetParameter(1);
		kaon_sigma = 3141.59/180*tmp_kaon_means->GetFunction("gaus")->GetParameter(2);
		*/
		printf("Mean photon resolutions of LUT\n");
		printf("Pion - Mean: %12.04f deg         Sigma: %12.04f mrad         %12.04f deg        Entries: %12.04f\n",pion_mean,pion_sigma,.057*pion_sigma,pion_lut_angles->GetMean());
		printf("Kaon - Mean: %12.04f deg         Sigma: %12.04f mrad         %12.04f deg        Entries: %12.04f\n",kaon_mean,kaon_sigma,.057*kaon_sigma,kaon_lut_angles->GetMean());
	}
	if (fill_d_midline_n > 0)
	{ 
		printf("Check Delta Midline\n");

	
		if (flatten_time==true || sep_updown == true || monochrome_plot == true)
		{
			printf("flatten_time, sep_updown, and monochrome_plot  not currently implemented for the line reconstruction.\n");
		}

		int num_line_points = 5000;
		//float time_spread = 4;
		float time_spread = 4;
		float dist_spread = 40;
		float max_dev_sq = 5;

		float callibration_y_shift = 0;

		max_dev_sq *= max_dev_sq;

		pion_beta = dirc_model->get_beta(energy,pimass);
		kaon_beta = dirc_model->get_beta(energy,kmass);

		std::vector<dirc_point> hit_points_pion;
		std::vector<dirc_point> hit_points_kaon;

		dirc_model->set_liquid_absorbtion(liquid_absorbtion);
		dirc_model->set_liquid_index(liquid_index);
		dirc_model->set_three_seg_mirror(three_seg_mirror);
		dirc_model->set_sidemirror(sm_xr,sm_xl);


		dirc_model->set_pmt_offset(pmt_offset);
		dirc_model->set_upper_wedge_angle_diff(wedge_uncertainty);
		dirc_model->set_bar_box_angle(bar_box_box_angle);

		//ns
		float pion_time = particle_flight_distance/(pion_beta*.3);
		float kaon_time = particle_flight_distance/(kaon_beta*.3);

		for (int i = 0; i < fill_d_midline_n; i++)
		{
			printf("\r                                                    ");
			printf("\rrunning iter %8d/%d  ",i+1,fill_d_midline_n);


			fflush(stdout);

			n_sim_phots = spread_ang.Gaus(mean_n_phot,spread_n_phot);
			//assume its a middle bar
			dirc_model->sim_rand_n_photons(\
					sim_points,\
					n_sim_phots,\
					pion_angle,\
					1,\
					particle_x,\
					particle_y,\
					pion_time,\
					particle_theta+const_track_off,\
					particle_phi,\
					tracking_unc,\
					ckov_unc,\
					pion_beta);
			//		digitizer.digitize_points(sim_points);

			std::vector<float> phi_last_wall_neg;
			std::vector<float> phi_last_wall_pos;
			std::vector<dirc_point> provisional_points_lwn;
			std::vector<dirc_point> provisional_points_lwp;
			float pion_refrac_mod = 1.000;
			float kaon_refrac_mod = 1.000;

			float pion_emit_angle = 57.3*acos(1/(pion_beta*pion_refrac_mod*refrac_index));
			float kaon_emit_angle = 57.3*acos(1/(kaon_beta*kaon_refrac_mod*refrac_index));

			dirc_model->track_all_line_photons(\
					provisional_points_lwn,\
					provisional_points_lwp,\
					num_line_points,\
					pion_emit_angle,\
					particle_theta,\
					particle_phi,\
					particle_x,\
					particle_y,\
					-17.25/2,\
					pion_time,\
					1);

			float t_dx = -1;
			float t_dy = -1;
			float t_dt = -1;

			for (unsigned int i = 0; i < sim_points.size(); i++)
			{
				float min_dist_sq = 10000;
				float cur_dist_sq = -1;
				for (unsigned int j = 0; j < provisional_points_lwn.size(); j++)
				{
					cur_dist_sq = (sim_points[i].x - provisional_points_lwn[j].x)*(sim_points[i].x - provisional_points_lwn[j].x);		
					cur_dist_sq += (sim_points[i].y - provisional_points_lwn[j].y)*(sim_points[i].y - provisional_points_lwn[j].y);		
					cur_dist_sq /= dist_spread;
					cur_dist_sq += (sim_points[i].t - provisional_points_lwn[j].t)*(sim_points[i].t - provisional_points_lwn[j].t)/time_spread;
					if (cur_dist_sq < min_dist_sq)
					{
						min_dist_sq = cur_dist_sq;
						t_dx = sim_points[i].x - provisional_points_lwn[j].x;
						t_dy = sim_points[i].y - provisional_points_lwn[j].y;
						t_dt = sim_points[i].t - provisional_points_lwn[j].t;
					}	
				}
				for (unsigned int j = 0; j < provisional_points_lwp.size(); j++)
				{
					cur_dist_sq = (sim_points[i].x - provisional_points_lwp[j].x)*(sim_points[i].x - provisional_points_lwp[j].x);		
					cur_dist_sq += (sim_points[i].y - provisional_points_lwp[j].y)*(sim_points[i].y - provisional_points_lwp[j].y);		
					cur_dist_sq /= dist_spread;
					cur_dist_sq += (sim_points[i].t - provisional_points_lwp[j].t)*(sim_points[i].t - provisional_points_lwp[j].t)/time_spread;	
					if (cur_dist_sq < min_dist_sq)
					{
						min_dist_sq = cur_dist_sq;
						t_dx = sim_points[i].x - provisional_points_lwp[j].x;
						t_dy = sim_points[i].y - provisional_points_lwp[j].y;
						t_dt = sim_points[i].t - provisional_points_lwp[j].t;
					}	
				}
				//cur_mean_phi /= total_prov_points;
				//printf("pion/pion: %12.04f %12.04f\n",cur_mean_phi/total_prov_points,atan2(cur_mean_y,cur_mean_x));
				pion_midline_dx->Fill(t_dx);
				pion_midline_dy->Fill(t_dy);
				pion_midline_dt->Fill(t_dt);
			}


			phi_last_wall_neg.clear();
			phi_last_wall_pos.clear();
			provisional_points_lwn.clear();	
			provisional_points_lwp.clear();	

			sim_points.clear();

			dirc_model->sim_rand_n_photons(\
					sim_points,\
					n_sim_phots,\
					kaon_angle,\
					1,\
					particle_x,\
					particle_y,\
					kaon_time,\
					particle_theta+const_track_off,\
					particle_phi,\
					tracking_unc,\
					ckov_unc,\
					kaon_beta);


			phi_last_wall_neg.clear();
			phi_last_wall_pos.clear();
			provisional_points_lwn.clear();	
			provisional_points_lwp.clear();	

			dirc_model->track_all_line_photons(\
					provisional_points_lwn,\
					provisional_points_lwp,\
					num_line_points,\
					kaon_emit_angle,\
					particle_theta,\
					particle_phi,\
					particle_x,\
					particle_y,\
					-17.25/2,\
					kaon_time,\
					1);
			for (unsigned int i = 0; i < sim_points.size(); i++)
			{
				float min_dist_sq = 10000;
				float cur_dist_sq = -1;
				for (unsigned int j = 0; j < provisional_points_lwn.size(); j++)
				{
					cur_dist_sq = (sim_points[i].x - provisional_points_lwn[j].x)*(sim_points[i].x - provisional_points_lwn[j].x);		
					cur_dist_sq += (sim_points[i].y - provisional_points_lwn[j].y)*(sim_points[i].y - provisional_points_lwn[j].y);		
					cur_dist_sq /= dist_spread;
					cur_dist_sq += (sim_points[i].t - provisional_points_lwn[j].t)*(sim_points[i].t - provisional_points_lwn[j].t)/time_spread;	
					if (cur_dist_sq < min_dist_sq)
					{
						min_dist_sq = cur_dist_sq;
						t_dx = sim_points[i].x - provisional_points_lwn[j].x;
						t_dy = sim_points[i].y - provisional_points_lwn[j].y;
						t_dt = sim_points[i].t - provisional_points_lwn[j].t;
					}	
				}
				for (unsigned int j = 0; j < provisional_points_lwp.size(); j++)
				{
					cur_dist_sq = (sim_points[i].x - provisional_points_lwp[j].x)*(sim_points[i].x - provisional_points_lwp[j].x);		
					cur_dist_sq += (sim_points[i].y - provisional_points_lwp[j].y)*(sim_points[i].y - provisional_points_lwp[j].y);		
					cur_dist_sq /= dist_spread;
					cur_dist_sq += (sim_points[i].t - provisional_points_lwp[j].t)*(sim_points[i].t - provisional_points_lwp[j].t)/time_spread;	
					if (cur_dist_sq < min_dist_sq)
					{
						min_dist_sq = cur_dist_sq;
						t_dx = sim_points[i].x - provisional_points_lwp[j].x;
						t_dy = sim_points[i].y - provisional_points_lwp[j].y;
						t_dt = sim_points[i].t - provisional_points_lwp[j].t;
					}	
				}
				kaon_midline_dx->Fill(t_dx);
				kaon_midline_dy->Fill(t_dy);
				kaon_midline_dt->Fill(t_dt);
			}
		}
		printf("\nMidline Delta Run Completed, Calibration String:\n");
		printf("\npion_x_adj pion_y_adj kaon_x_adj kaon_y_adj\n");
		printf("%12.04f %12.04f %12.04f %12.04f\n",-pion_midline_dx->GetMean(),-pion_midline_dy->GetMean() - callibration_y_shift,-kaon_midline_dx->GetMean(),-kaon_midline_dy->GetMean() - callibration_y_shift);
	}
	if (line_output_n > 0)
	{
		printf("Line Output for %d Phi steps\n",sparse_sim_n);
		dirc_point out_val;
		//assume speed of light particle and "thin" cone.
		float emit_angle;
		float pion_beta = dirc_model->get_beta(energy,pimass);
		float wavelength = 0;
		std::vector<dirc_point> left_points;
		std::vector<dirc_point> right_points;
		emit_angle = 57.3*acos(pion_beta/(.996*refrac_index)) + spread_ang.Gaus(0,particle_theta_spread);
		emit_angle = dirc_model->get_cerenkov_angle_rand(pion_beta,0,wavelength);
		float particle_time = particle_flight_distance/(pion_beta*.3);

		//emit_z = -17.25/2;

		dirc_model->track_all_line_photons(\
				left_points,\
				right_points,\
				line_output_n,\
				emit_angle,\
				particle_theta,\
				particle_phi,\
				particle_x,\
				particle_y,\
				-17.25/2,\
				particle_time,\
				1);
		//last 3 are z, t, and bar;
		float x,y,t,phi;
		float last_x,last_y,last_t,last_phi;
		float x_lim,y_lim,t_lim,phi_lim;

		x = 0;
		y = 0;
		t = 0;
		phi = 0;
		last_x = 0;
		last_y = 0;
		last_t = 0;
		last_phi = 0;
		x_lim = 10;
		y_lim = 10;
		t_lim = 10;
		phi_lim = .1;
		for (unsigned int i = 0; i < left_points.size(); i++)
		{
			out_val = left_points[i];
			x = out_val.x; 
			y = out_val.y; 
			t = out_val.t; 
			phi = out_val.init_phi; 
			if (fabs(x-last_x) > x_lim || fabs(y - last_y) > y_lim || fabs(t - last_t) > t_lim || fabs(phi - last_phi) > phi_lim)
			{
				//print a break in the data if there are any jumps
				printf("\n");
			}
			printf("%05d %12.04f %12.04f %12.04f %12.04f %12.04f\n",-1,emit_angle,out_val.init_phi,out_val.x,out_val.y,out_val.t);
			last_x = x; 
			last_y = y; 
			last_t = t; 
			last_phi = phi; 
		}
		for (unsigned int i = 0; i < right_points.size(); i++)
		{
			out_val = right_points[i];
			x = out_val.x; 
			y = out_val.y; 
			t = out_val.t; 
			phi = out_val.init_phi; 
			if (fabs(x-last_x) > x_lim || fabs(y - last_y) > y_lim || fabs(t - last_t) > t_lim || fabs(phi - last_phi) > phi_lim)
			{
				//print a break in the data;
				//printf("\n");
			}
			//printf("%05d %12.04f %12.04f %12.04f %12.04f %12.04f\n",-1,emit_angle,out_val.init_phi,out_val.x,out_val.y,out_val.t);
			last_x = x; 
			last_y = y; 
			last_t = t; 
			last_phi = phi; 
		}
	}
	if (bounces_phot_n > 0)
	{
		printf("Outputing N bounces histograms (With %d initial photons)\n",bounces_phot_n);
		dirc_model->set_store_bounces(true);
		std::vector<dirc_point> sim_points;
		
		dirc_model->sim_rand_n_photons(\
				sim_points,\
				bounces_phot_n,\
				0,\
				1,\
				particle_x,\
				particle_y,\
				0,\
				particle_theta,\
				particle_phi,\
				tracking_unc,\
				0,\
				pion_beta);
		
		std::vector<int> nbouncesx;
		std::vector<int> nbouncesz;
		std::vector<int> nbouncesx_dir;
		std::vector<int> nbouncesz_dir;
		std::vector<int> nbouncesx_indir;
		std::vector<int> nbouncesz_indir;

		dirc_model->fill_bounces_vecs(\
			nbouncesx,\
			nbouncesz,\
			nbouncesx_dir,\
			nbouncesz_dir,\
			nbouncesx_indir,\
			nbouncesz_indir);

		for (unsigned int i = 0; i < nbouncesx.size(); i++)
		{
			n_bounces_x->Fill(nbouncesx[i]);
		}
		for (unsigned int i = 0; i < nbouncesz.size(); i++)
		{
			n_bounces_z->Fill(nbouncesz[i]);
		}
		for (unsigned int i = 0; i < nbouncesx_dir.size(); i++)
		{
			n_bounces_x_direct->Fill(nbouncesx_dir[i]);
		}
		for (unsigned int i = 0; i < nbouncesz_dir.size(); i++)
		{
			n_bounces_z_direct->Fill(nbouncesz_dir[i]);
		}
		for (unsigned int i = 0; i < nbouncesx_indir.size(); i++)
		{
			n_bounces_x_indirect->Fill(nbouncesx_indir[i]);
		}
		for (unsigned int i = 0; i < nbouncesz_indir.size(); i++)
		{
			n_bounces_z_indirect->Fill(nbouncesz_indir[i]);
		}
	
		dirc_model->set_store_bounces(false);
	}
	if (sim_time_test_n > 0)
	{
		printf("Timing %d single particle simulations\n",sim_time_test_n);
		clock_t tmp_clock = clock();
		std::vector<dirc_point> sim_points;
		for (int i = 0; i < sim_time_test_n; i++)
		{
			pion_beta = dirc_model->get_beta(energy,pimass);
			kaon_beta = dirc_model->get_beta(energy,kmass);

			dirc_model->sim_rand_n_photons(\
					sim_points,\
					n_sim_phots,\
					0,\
					1,\
					particle_x,\
					particle_y,\
					0,\
					particle_theta+const_track_off,\
					particle_phi,\
					tracking_unc,\
					ckov_unc,\
					pion_beta);
		}
		tmp_clock = clock() - tmp_clock ;
		float time_taken = ((float)tmp_clock)/(CLOCKS_PER_SEC);
		float per_particle_time = time_taken/sim_time_test_n;
		float sim_rate = 1/per_particle_time;

		printf("Total Time Taken: %12.02es\n",time_taken);
		printf("Time per event:   %12.02ems\n",per_particle_time * 1000);
		printf("Rate:             %12.02eHz\n",sim_rate);
	}
	if (out_csv == true)
	{
		int out_csv_n = 40000;

		
		std::vector<dirc_point> store_points;
		pion_beta = dirc_model->get_beta(energy,pimass);
		kaon_beta = dirc_model->get_beta(energy,kmass);
		dirc_model->sim_rand_n_photons(\
                        store_points,\
                        out_csv_n,\
                        pion_angle,\
                        1,\
                        particle_x,\
                        particle_y,\
                        0,\
                        particle_theta,\
                        particle_phi,\
                        tracking_unc,\
                        ckov_unc,\
                        pion_beta);


		std::ofstream csv_out;
		csv_out.open("pion_loc.csv");

		char outline[256];
		
		for (unsigned int i = 0; i < store_points.size(); i++)
		{
			sprintf(outline,"%12.04f	%12.04f	%12.04f\n",store_points[i].x,store_points[i].y,store_points[i].t);
			csv_out << outline;
		}
		

		csv_out.close();
	}

	if (output_box_angles_n > 0)
	{
		std::vector<dirc_point> store_points;
		std::vector<float> upper_wedge_angle;
		std::vector<float> focus_angle;
		std::vector<float> large_flat_angle;
		std::vector<float> side_angle;

		dirc_model->set_upper_wedge_angle_store(true);
		dirc_model->set_store_optical_angles(true);

		pion_beta = dirc_model->get_beta(energy,pimass);
		kaon_beta = dirc_model->get_beta(energy,kmass);
		if (monochrome_plot == true)
		{
			pion_angle = rad_to_deg*acos(1/(refrac_index*pion_beta));
			kaon_angle = rad_to_deg*acos(1/(refrac_index*kaon_beta));
			pion_beta = -1;
			kaon_beta = -1;
			ckov_unc = 0;
		}
		dirc_model->sim_rand_n_photons(\
			store_points,\
			output_box_angles_n,\
			pion_angle,\
			1,\
			particle_x,\
			particle_y,\
			0,\
			particle_theta,\
			particle_phi,\
			tracking_unc,\
			ckov_unc,\
			pion_beta);

		upper_wedge_angle = dirc_model->get_upper_wedge_incident();
		focus_angle = dirc_model->get_focus_photon_angles();
		side_angle = dirc_model->get_side_photon_angles();
		large_flat_angle = dirc_model->get_large_flat_photon_angles();

		for (unsigned int i = 0; i < upper_wedge_angle.size(); i++)
		{
			upper_wedge_reflect_angle->Fill(upper_wedge_angle[i]);
		}
		for (unsigned int i = 0; i < focus_angle.size(); i++)
		{
			focus_reflect_angle->Fill(focus_angle[i]);
		}
		for (unsigned int i = 0; i < side_angle.size(); i++)
		{
			side_reflect_angle->Fill(side_angle[i]);
		}
		for (unsigned int i = 0; i < large_flat_angle.size(); i++)
		{
			large_flat_reflect_angle->Fill(large_flat_angle[i]);
		}

		dirc_model->set_upper_wedge_angle_store(false);
	}


			
	if (refraction_sim_n > 0)
	{
		std::vector<std::pair<float, float> > pion_theta_cphi;
		std::vector<std::pair<float, float> > kaon_theta_cphi;

		particle_theta=particle_theta_mean;	
		particle_phi=particle_phi_mean;	

		printf("Creating Refraction histograms...\n");
		std::vector<float> before_fill;
		std::vector<float> after_fill;
		std::vector<float> pmt_incidence;

		pion_theta_cphi = dirc_model->get_refraction_rand_phi(\
				before_fill,\
				after_fill,\
				pmt_incidence,\
				refraction_sim_n,\
				pion_angle ,
				particle_x,\
				particle_y,\
				particle_theta,\
				particle_phi,\
				0,\
				0,\
				pion_beta);
		//Since we're running this all at once, ignore tracking and other uncertainties

		for (unsigned int i = 0; i < before_fill.size(); i++)
		{
			ref_pion_before->Fill(before_fill[i]*rad_to_deg);
		}
		for (unsigned int i = 0; i < after_fill.size(); i++)
		{
			ref_pion_after->Fill(after_fill[i]*rad_to_deg);
		}
		for (unsigned int i = 0; i < pmt_incidence.size(); i++)
		{
			ref_pion_sens_plane->Fill(pmt_incidence[i]*rad_to_deg);
		}
		for (unsigned int i = 0; i < pion_theta_cphi.size(); i++)
		{
			ref_theta_cphi_pion->Fill(pion_theta_cphi[i].first*rad_to_deg,pion_theta_cphi[i].second*rad_to_deg);
		}

		kaon_theta_cphi = dirc_model->get_refraction_rand_phi(\
				before_fill,\
				after_fill,\
				pmt_incidence,\
				refraction_sim_n,\
				kaon_angle ,
				particle_x,\
				particle_y,\
				particle_theta,\
				particle_phi,\
				0,\
				0,\
				kaon_beta);


		for (unsigned int i = 0; i < before_fill.size(); i++)
		{
			ref_kaon_before->Fill(before_fill[i]*rad_to_deg);
		}
		for (unsigned int i = 0; i < after_fill.size(); i++)
		{
			ref_kaon_after->Fill(after_fill[i]*rad_to_deg);
		}
		for (unsigned int i = 0; i < pmt_incidence.size(); i++)
		{
			ref_kaon_sens_plane->Fill(pmt_incidence[i]*rad_to_deg);
		}
		for (unsigned int i = 0; i < kaon_theta_cphi.size(); i++)
		{
			ref_theta_cphi_kaon->Fill(kaon_theta_cphi[i].first*rad_to_deg,kaon_theta_cphi[i].second*rad_to_deg);
		}
		printf("Refraction histograms created\n");
	}


	float t_lambda = 0;
	for (int i = 0; i < 1000000; i++)
	{
		float pion_beta = dirc_model->get_beta(energy,pimass);
		float kaon_beta = dirc_model->get_beta(energy,kmass);
		
		//printf("%12.04f\n",dirc_model->get_cerenkov_angle_rand(pion_beta,0,t_lambda));

		pion_cerenkov->Fill(dirc_model->get_cerenkov_angle_rand(pion_beta,0,t_lambda));
		pion_lambda->Fill(t_lambda);
		kaon_cerenkov->Fill(dirc_model->get_cerenkov_angle_rand(kaon_beta,0,t_lambda));
		kaon_lambda->Fill(t_lambda);
	}
	if (output_peak_lambda == true)
	{
		printf("Pion gamma wavelength : kaon gamma wavelength:\n");
		printf("%12.04f %12.04f %12.04f %12.04f\n",energy/pimass,pion_lambda->GetMean(),energy/kmass,kaon_lambda->GetMean());
	}
	tfile->cd();

	pion_dist_geant_xy->Write();
	pion_dist_geant_xt->Write();
	pion_dist_geant_yt->Write();
	pion_dist_geant_x->Write();
	pion_dist_geant_y->Write();
	pion_dist_geant_t->Write();

	pion_dist_x->Write();
	pion_dist_y->Write();
	kaon_dist_x->Write();
	kaon_dist_y->Write();

	box_check_x->Write();
	box_check_y->Write();
	box_check_d->Write();
	box_check_xy->Write();

	pion_dist_xy->Write();
	pion_dist_xy_geantbins->Write();
	pion_dist_xy_geantbins_cm->Write();
	kaon_dist_xy->Write();
	pion_dist_xt->Write();
	kaon_dist_xt->Write();
	pion_dist_yt->Write();
	kaon_dist_yt->Write();
	pion_dist_t->Write();
	kaon_dist_t->Write();
	ll_diff_pion->Write();
	ll_diff_kaon->Write();
	ll_diff_pion_up->Write();
	ll_diff_kaon_up->Write();
	ll_diff_pion_down->Write();
	ll_diff_kaon_down->Write();
	phot_found_pion->Write();
	phot_found_kaon->Write();
	phot_found_pion_up->Write();
	phot_found_kaon_up->Write();
	phot_found_pion_down->Write();
	phot_found_kaon_down->Write();
	pion_cerenkov->Write();
	kaon_cerenkov->Write();
	pion_coverage_xy->Write();
	kaon_coverage_xy->Write();
	kinematics_yields->Write();
	liquid_dist->Write();
	pion_lambda->Write();
	kaon_lambda->Write();
	
	pion_unused_photons->Write();
	kaon_unused_photons->Write();

	ref_pion_before->Write();
	ref_pion_after->Write();
	ref_kaon_before->Write();
	ref_kaon_after->Write();
	ref_pion_sens_plane->Write();
	ref_kaon_sens_plane->Write();

	n_bounces_x->Write();
	n_bounces_z->Write();
	n_bounces_x_direct->Write();
	n_bounces_z_direct->Write();
	n_bounces_x_indirect->Write();
	n_bounces_z_indirect->Write();

	pion_midline_dx->Write();
	pion_midline_dy->Write();
	pion_midline_dt->Write();
	kaon_midline_dx->Write();
	kaon_midline_dy->Write();
	kaon_midline_dt->Write();

	pion_lut_vals->Write();
	kaon_lut_vals->Write();
	pion_lut_means->Write();
	kaon_lut_means->Write();
	pion_lut_angles->Write();
	kaon_lut_angles->Write();
	pion_lut_dt_v_dang->Write();
	kaon_lut_dt_v_dang->Write();

	upper_wedge_reflect_angle->Write();
	focus_reflect_angle->Write();
	side_reflect_angle->Write();
	large_flat_reflect_angle->Write();

	ref_theta_cphi_pion->Write();
	ref_theta_cphi_kaon->Write();

	simulation_time->Write();
	tfile->Close();

	int status = 0;
	return status;

}
