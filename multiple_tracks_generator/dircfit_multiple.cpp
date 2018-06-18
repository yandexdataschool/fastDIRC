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
#include <TH1.h>
#include <TRandom3.h>
#include <TMinuit.h>

// TODO(kazeevn) make this block a proper class
// the code below relies on the particle types
// being from 0 to PARTICLE_NUMBER -1
const unsigned PARTICLE_NUMBER = 5;
enum ParticleTypes {
    Muon = 0,
    Pion = 1,
    Kaon = 2,
    Proton = 3,
    Electron = 4
};

const int PARTICLE_ANGLE = -1;

// GeV/c^2
const std::array<float, PARTICLE_NUMBER> masses {
    .1057, .1396, .4937, .9382720813, 0.5109989461e-3};

const std::array<unsigned int, PARTICLE_NUMBER> particle_frequencies {
    5, 75, 15, 5, 2};

int main(int nargs, char* argv[]) {  
	float energy = 5.0;
	float energy_mean = energy;
	float energy_spread = 0.1;
	std::array<float, PARTICLE_NUMBER> betas;
	std::array<float, PARTICLE_NUMBER> times;
	std::array<std::unique_ptr<DircSpreadGaussian>, PARTICLE_NUMBER> pdfs;
	std::mt19937 random_generator;
	std::discrete_distribution<> particle_type_generator(
            particle_frequencies.begin(), particle_frequencies.end());
	
	// We have tracking!
	const float particle_x = 0;
	const float particle_y = 0;
	const float particle_x_mean = particle_x;
	const float particle_y_mean = particle_y;
	// Only for particle two
	float particle_x_spread = 1000.;
	float particle_y_spread = 1000.;
	float particle_theta = 4;
	float particle_theta_mean = particle_theta;
	float particle_theta_spread = 0;
	float particle_phi = 40;
	float const_track_off = 0;

	float particle_flight_distance = 0;

	bool use_moliere_scattering = false;
	int num_runs = 1000;
	float mean_n_phot = 40;
	float spread_n_phot = 0;

	float wedge_uncertainty = 0/57.3;
	float mirror_angle_change = 0;
	float mirror_angle_change_unc = 0;
	float mirror_angle_change_yunc = 0;
	float box_rot = 0;
	float box_rot_unc = 0;
	float bar_box_box_angle = 0/57.3;
	/* 1200 - 400 = 800.  Changed 05/09/2016.  Does not affect
	   threeseg mirror reconstruction as far as I can tell - this
	   was known. */
	float mirror_r_difference = 400;
	float wedge_non_uniformity = 0;
	float pmt_offset = 0;
	float main_mirror_nonuniformity = 0;
	float foc_mirror_size = 288;

	float pmt_min_z = -1000;
	float pmt_max_z = 1000;
	float large_mirror_min_z = -1000;
	float large_mirror_max_z = 1000;

	// Set boundaries for photons to optical plane in large box
	pmt_min_z = -559;
	pmt_max_z = -329;
	large_mirror_min_z = -559;
	large_mirror_max_z = -130;

	float upper_wedge_yang_spread = 0;
	int rseed = 1337;

	float tracking_unc = .0000*57.3; //mrad
	// float ckov_unc = .0077*57.3; //chromatic + optical aberation = 7.7mrad
	float ckov_unc = .003*57.3; //transport = 3mrad

	float resx = 6;
	float resy = 6;
	float minx = -1500;
	float maxx = 1500;
	float miny = -500;
	float maxy = 500;
	float t_unc = .27;
	float t_bin_size = 1;

	float digit_miny = -50;
	float digit_maxy = 300;

	digit_miny = miny;
	digit_maxy = maxy;

	// Sets the side boundaries of the distributions
	float sm_xl = -10000000;
	float sm_xr = -sm_xl;

	float s_func_x = 6;
	float s_func_y = s_func_x;
	float s_func_t = 1.0;
	float sfunc_sig = 1;

	int n_sim_phots = 40;

	int n_phi_phots = 150000;
	int n_z_phots = 4;

	bool use_quartz_for_liquid = false;
	bool three_seg_mirror = true;

	float liquid_absorbtion = 0*-log(.7)/1000;
	float liquid_index = 1.33;

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
			// run with SLAC fdirc prototype geometry
			three_seg_mirror = false;
			mirror_r_difference = 0;
			mean_n_phot = 32.4;
			spread_n_phot = 6;
			liquid_index = 1.47;
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

	float main_mirror_angle = 74.11+mirror_angle_change;
	float pdf_unc_red_fac = 1;
	TRandom3 spread_ang(rseed+3);
	auto dirc_model = std::make_unique<DircThreeSegBoxSim>(
			rseed,
			-1200 + mirror_r_difference,
			foc_mirror_size,
			main_mirror_angle,
			600,
			47.87 + box_rot + mirror_angle_change);
	dirc_model->set_store_traveled(false); // uses LOTS of memory if set to true.
	dirc_model->set_liquid_index(liquid_index);
	dirc_model->set_wedge_mirror_rand(wedge_non_uniformity);
	dirc_model->set_three_seg_mirror(three_seg_mirror);
	dirc_model->set_pmt_plane_zs(pmt_min_z, pmt_max_z);
	dirc_model->set_large_mirror_zs(large_mirror_min_z, large_mirror_max_z);
	dirc_model->set_use_quartz_n_for_liquid(use_quartz_for_liquid);

	// news are handled by the ROOT memory management
	TFile* tfile = new TFile(rootfilename, "RECREATE");
	TH1I* particle_one_type_th1i = new TH1I("particle_one_type",
					   "Type of the signal particle which is being measured", 
					   PARTICLE_NUMBER, 0, PARTICLE_NUMBER - 1);
	TH1I* particle_two_type_th1i = new TH1I("particle_two_type",
					   "Type of the noise particle",
					   PARTICLE_NUMBER, 0, PARTICLE_NUMBER - 1);
	TH1F* distance = new TH1F("distance",
				  "Distance between the signal and noise racks", 10, 0, 4000);
	std::array<TH1F*, PARTICLE_NUMBER> dlls_th1f;

        dlls_th1f[ParticleTypes::Kaon] = new TH1F("dll_kaon",
						 "LL(kaon) - LL(pion)", 10, -20, 20);
	dlls_th1f[ParticleTypes::Proton] = new TH1F("dll_proton",
						    "LL(proton) - LL(pion)", 10, -1000, 1000);
	dlls_th1f[ParticleTypes::Muon] = new TH1F("dll_muon",
				  "LL(muon) - LL(pion)", 10, -20, 20);
	dlls_th1f[ParticleTypes::Electron] = new TH1F("dll_electron",
				    "LL(electron) - LL(pion)", 10, -20, 20);
	dlls_th1f[ParticleTypes::Pion] = nullptr;
	maxy *= 5;
	DircRectDigitizer digitizer(
			minx,
			maxx,
			resx,
			digit_miny,
			digit_maxy,
			resy,
			t_unc,
			t_bin_size);

	printf("Beginning Run\n");
	dirc_model->set_focmirror_nonuniformity(main_mirror_nonuniformity);
	dirc_model->set_use_moliere(use_moliere_scattering);
	// assume momentum is the same for both for now - high energy;
	dirc_model->set_moliere_p(energy*1000);

	dirc_model->set_liquid_absorbtion(liquid_absorbtion);
	dirc_model->set_liquid_index(liquid_index);
	dirc_model->set_three_seg_mirror(three_seg_mirror);
	dirc_model->set_sidemirror(sm_xr,sm_xl);

	dirc_model->set_pmt_offset(pmt_offset);
	dirc_model->set_upper_wedge_angle_diff(wedge_uncertainty);
	dirc_model->set_bar_box_angle(bar_box_box_angle);

	// conmpute and intialize the pdfs
	for (size_t particle = 0; particle < PARTICLE_NUMBER; ++particle) {
	    betas[particle] = dirc_model->get_beta(energy, masses[particle]);
	    // ns
	    times[particle] = particle_flight_distance/(betas[particle]*.3);
	    std::vector<dirc_point> hit_points;

	    dirc_model->sim_reg_n_photons(hit_points,
					  n_phi_phots,
					  n_z_phots,
					  -1,
					  1,
					  particle_x,
					  particle_y,
					  times[particle],
					  particle_theta,
					  particle_phi,
					  0,
					  ckov_unc/pdf_unc_red_fac,
					  betas[particle]);
	    pdfs[particle] = std::make_unique<DircSpreadGaussian>(
	        sfunc_sig, hit_points, s_func_x, s_func_y, s_func_t);
	}
	for (int i = 0; i < num_runs; i++) {
	    printf("\r                                                    ");
	    printf("\rrunning iter %8d/%d  ", i+1, num_runs);
	    fflush(stdout);
	    dirc_model->set_focus_mirror_angle(
	        spread_ang.Gaus(main_mirror_angle, mirror_angle_change_unc),
		spread_ang.Gaus(0, mirror_angle_change_yunc));
	    dirc_model->set_upper_wedge_angle_diff(spread_ang.Gaus(0,wedge_uncertainty),
						   spread_ang.Gaus(0,upper_wedge_yang_spread));
	    dirc_model->set_bar_box_angle(spread_ang.Gaus(0,box_rot_unc));

	    if (particle_theta_mean < .01) {
		particle_phi = spread_ang.Uniform(0, 360);
	    }
	    particle_theta = spread_ang.Gaus(particle_theta_mean, particle_theta_spread);
	    energy = spread_ang.Gaus(energy_mean,energy_spread);
	    n_sim_phots = spread_ang.Gaus(mean_n_phot,spread_n_phot);
	    // We want more or less the same number of 
	    // signal particles of each type
	    const unsigned int particle_one_type = spread_ang.Integer(PARTICLE_NUMBER);
	    particle_one_type_th1i->Fill(particle_one_type);
	    // assume its a middle bar
	    // The first particle is always (0, 0) - we have tracking!
	    std::vector<dirc_point> sim_points;
	    dirc_model->sim_rand_n_photons(sim_points,
					   n_sim_phots,
					   PARTICLE_ANGLE,
					   1,
					   0.,
					   0.,
					   times[particle_one_type],
					   particle_theta+const_track_off,
					   particle_phi,
					   tracking_unc,
					   ckov_unc,
					   betas[particle_one_type]);
	    std::vector<dirc_point> hits_second;
	    particle_theta = spread_ang.Gaus(particle_theta_mean, particle_theta_spread);
	    const float particle_two_x = spread_ang.Gaus(particle_x_mean, particle_x_spread);
	    const float particle_two_y = spread_ang.Gaus(particle_y_mean, particle_y_spread);
	    // TODO(kazeevn) square?
	    distance->Fill(sqrt(particle_two_x*particle_two_x + particle_two_y*particle_two_y));
	    energy = spread_ang.Gaus(energy_mean, energy_spread);
	    n_sim_phots = spread_ang.Gaus(mean_n_phot, spread_n_phot);
	    // For the noise particle, we want an LHCb-like distribution
	    // also, since TRandom3 doesn't provide weights, we use the
	    // standard library
	    const unsigned int particle_two_type = particle_type_generator(random_generator);
	    particle_two_type_th1i->Fill(particle_two_type);
	    dirc_model->sim_rand_n_photons(hits_second,
					   n_sim_phots,
					   PARTICLE_ANGLE,
					   1,
					   particle_two_x,
					   particle_two_y,
					   times[ParticleTypes::Kaon],
					   particle_theta + const_track_off,
					   particle_phi,
					   tracking_unc,
					   ckov_unc,
					   betas[ParticleTypes::Kaon]);
	    // TODO(kazeevn) switch from copying memory
	    sim_points.insert(sim_points.end(), hits_second.begin(), hits_second.end());
	    digitizer.digitize_points(sim_points);
	    const float ll_pion = pdfs[ParticleTypes::Pion]->get_log_likelihood(sim_points);
	    for (size_t particle = 0; particle < PARTICLE_NUMBER; ++particle) {
		if (particle == ParticleTypes::Pion) {
		    continue;
		}
		dlls_th1f[particle]->Fill(pdfs[particle]->get_log_likelihood(sim_points) - ll_pion);
	    }
	}
	printf("\nRun Completed\n");
	tfile->cd();
	distance->Write();
	particle_one_type_th1i->Write();
	particle_two_type_th1i->Write();
	for (size_t particle = 0; particle < PARTICLE_NUMBER; ++particle) {
	    if (particle == ParticleTypes::Pion) {
		continue;
	    }
	    dlls_th1f[particle]->Write();
	}
	tfile->Close();
	return 0;
}
