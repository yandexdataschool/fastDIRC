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
#include "../include/dirc_babar_digitizer.h"
#include "../include/dirc_babar_sim.h"
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


//#define USE_BABAR_BOX

int main(int nargs, char* argv[])
{  
	bool out_csv = false;

	float energy = 5.0;
	float energy_mean = energy;
	float energy_spread = 0;
	float kmass = .4937;
	float pimass = .1396;
	//float mumass = .1057;

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
	float const_track_off = 0;

	float particle_flight_distance = 0;


	bool monochrome_plot = false;
	bool flatten_time = false;
	bool sep_updown = false;
	bool output_peak_lambda = false;
	bool lut_slac = false;
	bool perform_chromatic_correction = true;
	bool use_moliere_scattering = false;

	int num_runs = 1000;
	int sparse_sim_n = -1;
	int line_output_n = -1;
	int lut_sim_n = -1;
	int gaus_ll_n = -1;
	int sim_time_test_n = -1;

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



	float pmt_min_z = -1000;
	float pmt_max_z = 1000;
	float large_mirror_min_z = -1000;
	float large_mirror_max_z = 1000;

	//Set boundaries for photons to optical plane in large box
	pmt_min_z = -559;
	pmt_max_z = -329;
	large_mirror_min_z = -559;
	large_mirror_max_z = -130;

	float upper_wedge_yang_spread = 0;
	int rseed = 1337;

	float tracking_unc = .0000*57.3; //mrad
	// 	float ckov_unc = .0077*57.3; //chromatic + optical aberation = 7.7mrad
	float ckov_unc = .003*57.3; //transport = 3mrad

	float resx = 6;
	float resy = 6;
	float rest = 1;
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

	//int n_phi_phots = 900000;
	int n_phi_phots = 150000;
	int n_z_phots = 4;

	bool use_quartz_for_liquid = false;
	bool three_seg_mirror = true;
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
	sprintf(rootfilename,"fitdirc.root");	

	printf("Arguments Passed=%d\n",nargs);

	for (int i = 1; i < nargs; i++)
	{
		if (strcmp(argv[i], "-of") == 0)
		{
			i++;
			sprintf(rootfilename,"%s",argv[i]);	
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
		else if (strcmp(argv[i], "-side_mirror") == 0)
		{
			i++;
			sm_xl = atof(argv[i]);
			i++;
			sm_xr = atof(argv[i]);
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
		else if (strcmp(argv[i], "-use_quartz_for_liquid") == 0)
		{
			use_quartz_for_liquid = true;
		}
		else if (strcmp(argv[i], "-use_moliere_scattering") == 0)
		{
			use_moliere_scattering = true;
			printf("Enabling Moliere Scattering - only implemented in loop mode\n");
		}
		else if (strcmp(argv[i], "-out_csv") == 0)
		{
			out_csv = true;
		}
		else if (strcmp(argv[i], "-coverage_plot") == 0)
		{
			coverage_plot = true;
		}
		else if (strcmp(argv[i], "-no_chromatic_correction") == 0)
		{
			perform_chromatic_correction = false;
		}
		else if (strcmp(argv[i], "-monochrome_light") == 0)
		{
			monochrome_plot = true;
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
		}
	}


	float main_mirror_angle = 74.11+mirror_angle_change;

	float rad_to_deg = 57.2958;

	float res_enhance = 1;

	if (flatten_time == true)
	{
		s_func_t = 100000000;
	}	

	float outcsv_x,outcsv_y;
	outcsv_x = 0*35;//bars are 35mm wide
	outcsv_y = 0;//mm


	if(out_csv){
		n_phi_phots = 4000;
		n_z_phots = 8;
		particle_y += outcsv_y;}
	float pdf_unc_red_fac = 1;


	TRandom3 spread_ang(rseed+3);

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
	dirc_model->set_pmt_plane_zs(pmt_min_z,pmt_max_z);
	dirc_model->set_large_mirror_zs(large_mirror_min_z,large_mirror_max_z);
	dirc_model->set_use_quartz_n_for_liquid(use_quartz_for_liquid);


	float pion_beta, kaon_beta/*, electron_beta:=1*/;
	pion_beta=kaon_beta=-1;
	float pion_angle, kaon_angle;
	pion_angle=kaon_angle = -1;


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

	TH2F *ref_theta_cphi_pion = new TH2F("ref_theta_cphi_pion","Emitted angle of Pion Photons versus angle into interface", 7200,0,360,1800,0,90);
	TH2F *ref_theta_cphi_kaon = new TH2F("ref_theta_cphi_kaon","Emitted angle of Kaon Photons versus angle into interface", 7200,0,360,1800,0,90);

	//TH2F *pion_dist_xy = new TH2F("pion_dist_xy","xy val of intercepted points - pion",(maxx-minx)/(res_enhance*resx),minx,maxx,(maxy-miny)/(res_enhance*resy),miny,maxy);
	TH2F *pion_dist_xy = new TH2F("pion_dist_xy","xy val of intercepted points - pion",500,-1500,1500,67,-150,252);
	TH2F *kaon_dist_xy = new TH2F("kaon_dist_xy","xy val of intercepted points - kaon",(maxx-minx)/(res_enhance*resx),minx,maxx,(maxy-miny)/(res_enhance*resy),miny,maxy);

	TH2F *box_check_xy = new TH2F("box_check_xy","xy val of points from top of wedge",(maxx-minx)/(res_enhance*resx),minx,maxx,(maxy-miny)/(res_enhance*resy),miny,maxy);

	TH2F *pion_dist_xt = new TH2F("pion_dist_xt","xt val of intercepted points - pion",(maxx-minx)/(res_enhance*resx),minx,maxx,(maxt-mint)/(res_enhance*rest),mint,maxt);
	TH2F *kaon_dist_xt = new TH2F("kaon_dist_xt","xt val of intercepted points - kaon",(maxx-minx)/(res_enhance*resx),minx,maxx,(maxt-mint)/(res_enhance*rest),mint,maxt);

	TH2F *pion_dist_yt = new TH2F("pion_dist_yt","yt val of intercepted points - pion",(maxy-miny)/(res_enhance*resy),miny,maxy,(maxt-mint)/(res_enhance*rest),mint,maxt);
	TH2F *kaon_dist_yt = new TH2F("kaon_dist_yt","yt val of intercepted points - kaon",(maxy-miny)/(res_enhance*resy),miny,maxy,(maxt-mint)/(res_enhance*rest),mint,maxt);


	maxy *= 5;

	TH2F *pion_coverage_xy = new TH2F("pion_coverage_xy","xy val of generated points - pion",(maxx-minx)/(res_enhance*resx),minx,maxx,(maxy-miny)/(res_enhance*resy),miny,maxy);
	TH2F *kaon_coverage_xy = new TH2F("kaon_coverage_xy","xy val of generated points - kaon",(maxx-minx)/(res_enhance*resx),minx,maxx,(maxy-miny)/(res_enhance*resy),miny,maxy);

	TH1F *pion_cerenkov = new TH1F("pion_cerenkov","pion cerenkov angle distribution",300,45,48);
	TH1F *kaon_cerenkov = new TH1F("kaon_cerenkov","kaon cerenkov angle distribution",300,45,48);
	TH1F *pion_lambda = new TH1F("pion_lambda","pion wavelength distribution",450,250,700);
	TH1F *kaon_lambda = new TH1F("kaon_lambda","kaon wavelength distribution",450,250,700);

	TH1F *pion_unused_photons = new TH1F("pion_unused_photons","pion photonons not used on pion data, pion pdf",100,-.5,99.5);
	TH1F *kaon_unused_photons = new TH1F("kaon_unused_photons","kaon photonons not used on kaon data, kaon pdf",100,-.5,99.5);

	TH1F *liquid_dist = new TH1F("liquid_dist","distance travelled in liquid (mm)",1500,0,1500);

	TH1F *simulation_time = new TH1F("simulation_time","Simulation time per particle",4001,-.5,4000.5);

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
	std::vector<dirc_point> sim_points;
	std::vector<dirc_point> confound_points;
	dirc_model->set_focmirror_nonuniformity(main_mirror_nonuniformity);
	if (num_runs > 0)
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
		printf("Performing LUT analysis with: %d points\n", (int) lut_size);

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

			//Do an iterative mean finding for ceter of oval?
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
					0,//t
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

	pion_dist_x->Write();
	pion_dist_y->Write();
	kaon_dist_x->Write();
	kaon_dist_y->Write();

	box_check_x->Write();
	box_check_y->Write();
	box_check_d->Write();
	box_check_xy->Write();

	pion_dist_xy->Write();
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

	pion_lut_vals->Write();
	kaon_lut_vals->Write();
	pion_lut_means->Write();
	kaon_lut_means->Write();
	pion_lut_angles->Write();
	kaon_lut_angles->Write();
	pion_lut_dt_v_dang->Write();
	kaon_lut_dt_v_dang->Write();

	ref_theta_cphi_pion->Write();
	ref_theta_cphi_kaon->Write();

	simulation_time->Write();
	tfile->Close();

	int status = 0;
	return status;

}
