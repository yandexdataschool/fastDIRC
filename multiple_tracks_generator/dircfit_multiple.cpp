#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <utility>

#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "../include/dirc_optical_sim.h"
#include "../include/dirc_threesegbox_sim.h"
#include "../include/dirc_point.h"
#include "../include/dirc_probability_spread.h"
#include "../include/dirc_probability_separation.h"
#include "../include/dirc_spread_radius.h"
#include "../include/dirc_spread_relative.h"
#include "../include/dirc_spread_linear_soft.h"
#include "../include/dirc_spread_gaussian.h"
#include "../include/dirc_digitizer.h"
#include "../include/dirc_rect_digitizer.h"
#include "../include/dirc_progressive_separation.h"
#include "../include/dirc_gluex_lut_enum.h"
#include "../include/dirc_lut_enum.h"
#include "../include/dirc_lut.h"
#include <TFile.h>
#include <TTree.h>
#include <TH3.h>
#include <TH2.h>
#include <TH1.h>
#include <TF1.h>
#include <TRandom3.h>
#include <TMinuit.h>

int main(int nargs, char* argv[])
{  
	double energy = 5.0;
	double energy_mean = energy;
	double energy_spread = 0;
	double kmass = .4937;
	double pimass = .1396;
	//double mumass = .1057;

	double particle_x = 0;
	double particle_y = 0;
	double particle_x_mean = particle_x;
	double particle_y_mean = particle_y;
	double particle_x_spread = 0;
	double particle_y_spread = 0;
	double particle_theta = 4;
	double particle_theta_mean = particle_theta;
	double particle_theta_spread = 0;
	double particle_phi = 40;
	double const_track_off = 0;

	double particle_flight_distance = 0;

	bool use_moliere_scattering = false;
	int num_runs = 1000;
	double mean_n_phot = 40;
	double spread_n_phot = 0;

	double wedge_uncertainty = 0/57.3;
	double mirror_angle_change = 0;
	double mirror_angle_change_unc = 0;
	double mirror_angle_change_yunc = 0;
	double box_rot = 0;
	double box_rot_unc = 0;
	double bar_box_box_angle = 0/57.3;
	/* 1200 - 400 = 800.  Changed 05/09/2016.  Does not affect
	   threeseg mirror reconstruction as far as I can tell - this
	   was known. */
	double mirror_r_difference = 400;
	double wedge_non_uniformity = 0;
	double pmt_offset = 0;
	double main_mirror_nonuniformity = 0;
	double foc_mirror_size = 288;

	double pmt_min_z = -1000;
	double pmt_max_z = 1000;
	double large_mirror_min_z = -1000;
	double large_mirror_max_z = 1000;

	// Set boundaries for photons to optical plane in large box
	pmt_min_z = -559;
	pmt_max_z = -329;
	large_mirror_min_z = -559;
	large_mirror_max_z = -130;

	double upper_wedge_yang_spread = 0;
	int rseed = 1337;

	double tracking_unc = .0000*57.3; //mrad
	// double ckov_unc = .0077*57.3; //chromatic + optical aberation = 7.7mrad
	double ckov_unc = .003*57.3; //transport = 3mrad

	double resx = 6;
	double resy = 6;
	double minx = -1500;
	double maxx = 1500;
	double miny = -500;
	double maxy = 500;
	double t_unc = .27;
	double t_bin_size = 1;

	double digit_miny = -50;
	double digit_maxy = 300;

	digit_miny = miny;
	digit_maxy = maxy;

	// Sets the side boundaries of the distributions
	double sm_xl = -10000000;
	double sm_xr = -sm_xl;

	double s_func_x = 6;
	double s_func_y = s_func_x;
	double s_func_t = 1.0;
	double sfunc_sig = 1;

	int n_sim_phots = 40;

	int n_phi_phots = 150000;
	int n_z_phots = 4;

	bool use_quartz_for_liquid = false;
	bool three_seg_mirror = true;

	double liquid_absorbtion = 0*-log(.7)/1000;
	double liquid_index = 1.33;

	const unsigned MAX_FILE_NAME_LENGTH = 256;
	char rootfilename[MAX_FILE_NAME_LENGTH];
	sprintf(rootfilename, "fitdirc.root");	

	for (int i = 1; i < nargs; i++)
	{
		if (strcmp(argv[i], "-of") == 0)
		{
			i++;
			snprintf(rootfilename, MAX_FILE_NAME_LENGTH, "%s",argv[i]);	
		}
		else if (strcmp(argv[i], "-cylindrical_mirror") == 0)
		{
			three_seg_mirror = false;	
		}
		else if (strcmp(argv[i], "-particle_phi") == 0)
		{
			i++;
			particle_phi = atof(argv[i]);
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
		else if (strcmp(argv[i], "-side_mirror") == 0)
		{
			i++;
			sm_xl = atof(argv[i]);
			i++;
			sm_xr = atof(argv[i]);
		}
		else if (strcmp(argv[i], "-particle_theta") == 0)
		{
			i++;
			particle_theta = atof(argv[i]);
			particle_theta_mean = particle_theta;
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
		else if (strcmp(argv[i], "-slac_geometry") == 0)
		{
			//run with SLAC fdirc prototype geometry
			three_seg_mirror = false;
			mirror_r_difference = 0;
			//mean_n_phot = 31.1;
			mean_n_phot = 32.4;
			spread_n_phot = 6;
			liquid_index = 1.47;
			//	sm_xl = -300;
			//	sm_xr = sm_xl + 1000;
			miny = -1500;
			maxy = 1500;
			digit_miny = miny;
			digit_maxy = maxy;
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
		else if (strcmp(argv[i], "-wedge_uncertainty") == 0)
		{
			i++;
			wedge_uncertainty = atof(argv[i]);
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
		else
		{
			printf("Unrecognized argument: %s\n",argv[i]);
			exit(-1);
		}
	}


	double main_mirror_angle = 74.11+mirror_angle_change;
	double pdf_unc_red_fac = 1;
	TRandom3 spread_ang(rseed+3);
	auto dirc_model = std::make_unique<DircThreeSegBoxSim>(
			rseed,\
			-1200 + mirror_r_difference,\
			foc_mirror_size,\
			main_mirror_angle,\
			600,\
			47.87 + box_rot + mirror_angle_change);
	dirc_model->set_store_traveled(false);// uses LOTS of memory if set to true.
	dirc_model->set_liquid_index(liquid_index);
	dirc_model->set_wedge_mirror_rand(wedge_non_uniformity);
	dirc_model->set_three_seg_mirror(three_seg_mirror);
	dirc_model->set_pmt_plane_zs(pmt_min_z,pmt_max_z);
	dirc_model->set_large_mirror_zs(large_mirror_min_z,large_mirror_max_z);
	dirc_model->set_use_quartz_n_for_liquid(use_quartz_for_liquid);


	double pion_beta, kaon_beta/*, electron_beta:=1*/;
	pion_beta=kaon_beta=-1;
	double pion_angle, kaon_angle;
	pion_angle=kaon_angle = -1;

	// ROOT memory management
	TFile* tfile = new TFile(rootfilename, "RECREATE");
	TH1F* ll_diff_pion = new TH1F("ll_diff_pion",
				      "Difference of log likelihood real = pion",200000,-200,200);
	TH1F* ll_diff_kaon = new TH1F("ll_diff_kaon",
				      "Difference of log likelihood real = kaon",200000,-200,200);
	TH1F* phot_found_pion = new TH1F("phot_found_pion",
					 "number of photons found on pion angle", 1001,-.5,1000.5);
	TH1F* phot_found_kaon = new TH1F("phot_found_kaon",
					 "number of photons found on kaon angle", 1001,-.5,1000.5);

	maxy *= 5;

	DircRectDigitizer digitizer(\
			minx,\
			maxx,\
			resx,\
			digit_miny,\
			digit_maxy,\
			resy,\
			t_unc,\
			t_bin_size);

	printf("Beginning Run\n");
	double llc, llf, ll_diff;
	llc=llf=ll_diff=0;
	std::vector<dirc_point> sim_points;
	std::vector<dirc_point> confound_points;
	dirc_model->set_focmirror_nonuniformity(main_mirror_nonuniformity);
	if (num_runs > 0) {
		pion_beta = dirc_model->get_beta(energy,pimass);
		kaon_beta = dirc_model->get_beta(energy,kmass);

		std::vector<dirc_point> hit_points_pion;
		std::vector<dirc_point> hit_points_kaon;


		dirc_model->set_use_moliere(use_moliere_scattering);
		//assume momentum is the same for both for now - high energy;
		dirc_model->set_moliere_p(energy*1000);

		dirc_model->set_liquid_absorbtion(liquid_absorbtion);
		dirc_model->set_liquid_index(liquid_index);
		dirc_model->set_three_seg_mirror(three_seg_mirror);
		dirc_model->set_sidemirror(sm_xr,sm_xl);

		dirc_model->set_pmt_offset(pmt_offset);
		dirc_model->set_upper_wedge_angle_diff(wedge_uncertainty);
		dirc_model->set_bar_box_angle(bar_box_box_angle);

		//ns
		double pion_time = particle_flight_distance/(pion_beta*.3);
		double kaon_time = particle_flight_distance/(kaon_beta*.3);

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
		std::unique_ptr<DircSpreadGaussian> pdf_pion = std::make_unique<DircSpreadGaussian>(
				sfunc_sig,\
				hit_points_pion,\
				s_func_x,\
				s_func_y,\
				s_func_t);
		std::unique_ptr<DircSpreadGaussian> pdf_kaon = std::make_unique<DircSpreadGaussian>(
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
					particle_theta+const_track_off,\
					particle_phi,\
					tracking_unc,\
					ckov_unc,\
					kaon_beta);


			digitizer.digitize_points(sim_points);

			llc = pdf_pion->get_log_likelihood(sim_points);
			llf = pdf_kaon->get_log_likelihood(sim_points);

			ll_diff_kaon->Fill(1*(llc-llf));
			phot_found_kaon->Fill(sim_points.size());

		}

		printf("\nRun Completed\n");
	}
	tfile->cd();
	ll_diff_pion->Write();
	ll_diff_kaon->Write();
	phot_found_pion->Write();
	phot_found_kaon->Write();
	tfile->Close();
	return 0;
}
