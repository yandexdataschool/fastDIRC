#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <iterator>
#include <utility>
#include <stdlib.h>
#include <math.h>

#include <TFile.h>
#include <TMath.h>
#include <TTree.h>
#include <TRandom3.h>

#include "../include/dirc_threesegbox_sim.h"
#include "../include/dirc_spread_gaussian.h"
#include "../include/dirc_rect_digitizer.h"

// TODO(kazeevn) make this block a proper class
// the code in main relies on the particle types
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

// LHCb-esque
const std::array<unsigned int, PARTICLE_NUMBER> particle_frequencies {
    5, 75, 15, 5, 2};

int main(int nargs, char* argv[]) {  
	float energy_mean = 5.0;
	float energy_spread = 0.;
	//const float eta_min = -0.10001;
	//const float eta_max = -0.1;
	const float eta_mean = -0.07;
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
	float particle_x_spread = 50000.;
	float particle_y_spread = 50000.;
	float particle_phi = 40;
	float const_track_off = 0;
	
	// meters
	float particle_flight_distance = 0;

	unsigned int num_runs = 1000;
	float mean_n_phot = 40;
	float spread_n_phot = 0;

	float mirror_angle_change = 0;
	float box_rot = 0;
	float bar_box_box_angle = 0/57.3;
	/* 1200 - 400 = 800.  Changed 05/09/2016.  Does not affect
	   threeseg mirror reconstruction as far as I can tell - this
	   was known. */
	float mirror_r_difference = 400;
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

	int rseed = 1337;

	float tracking_unc = .0000*57.3; //mrad
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

	int n_phi_phots = 100000;
	int n_z_phots = 4;
	//const unsigned int kde_generation_iterations = 40;

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
		else if (strcmp(argv[i], "-use_quartz_for_liquid") == 0)
		{
			use_quartz_for_liquid = true;
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
		else if (strcmp(argv[i], "-energy_spread") == 0)
		{
			i++;
			energy_spread = atof(argv[i]);
		}
		else if (strcmp(argv[i], "-E") == 0)
		{
			i++;
			energy_mean = atof(argv[i]);
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
		else if (strcmp(argv[i], "-box_rot") == 0)
		{
			i++;
			box_rot = atof(argv[i]);
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
		else if (strcmp(argv[i], "-rseed") == 0)
		{
			i++;
			rseed = atof(argv[i]);
		}
		else
		{
			printf("Unrecognized argument: %s\n",argv[i]);
			exit(-1);
		}
	}
	float main_mirror_angle = 74.11 + mirror_angle_change;
	float pdf_unc_red_fac = 1;
	std::unique_ptr<TRandom3> spread_ang = std::make_unique<TRandom3>(rseed + 3);
	auto dirc_model = std::make_unique<DircThreeSegBoxSim>(
			rseed,
			-1200 + mirror_r_difference,
			foc_mirror_size,
			main_mirror_angle,
			600,
			47.87 + box_rot + mirror_angle_change);
	dirc_model->set_store_traveled(false); // uses LOTS of memory if set to true.
	dirc_model->set_liquid_index(liquid_index);
	dirc_model->set_wedge_mirror_rand(0.);
	dirc_model->set_three_seg_mirror(three_seg_mirror);
	dirc_model->set_pmt_plane_zs(pmt_min_z, pmt_max_z);
	dirc_model->set_large_mirror_zs(large_mirror_min_z, large_mirror_max_z);
	dirc_model->set_use_quartz_n_for_liquid(use_quartz_for_liquid);

	// news are handled by the ROOT memory management
	// trying to unique_ptr them breaks things
	TFile* tfile = new TFile(rootfilename, "RECREATE");
	TTree* tree = new TTree("fastDIRC", "PID simulation via fastDIRC");
	// TODO(kazeevn) make UCHar
	UChar_t particle_one_type, particle_two_type;
	tree->Branch("particle_one_type", &particle_one_type, 
		     "Type of the signal particle which is being measured/b");
	tree->Branch("particle_two_type", &particle_two_type, 
		     "Type of the noise particle/b");
	Float_t tracks_distance;
	tree->Branch("tracks_distance", &tracks_distance,
		     "Distance between the signal and noise tracks/F");
	Float_t particle_one_energy, particle_one_eta;
	tree->Branch("particle_one_energy", &particle_one_energy,
		     "Energy of the signal particle, GeV/F");
	tree->Branch("particle_one_eta", &particle_one_eta,
		     "Pseudorapidity of particle one/F");
	std::array<Float_t, PARTICLE_NUMBER> dlls;
	tree->Branch("dll_kaon", &(dlls[ParticleTypes::Kaon]), "LL(kaon) - LL(pion)/F");
	tree->Branch("dll_muon", &(dlls[ParticleTypes::Muon]), "LL(muon) - LL(pion)/F");
	tree->Branch("dll_electron", &(dlls[ParticleTypes::Electron]), "LL(electron) - LL(pion)/F");
	tree->Branch("dll_proton", &(dlls[ParticleTypes::Proton]), "LL(proton) - LL(pion)/F");
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

	dirc_model->set_focmirror_nonuniformity(main_mirror_nonuniformity);

	dirc_model->set_liquid_absorbtion(liquid_absorbtion);
	dirc_model->set_liquid_index(liquid_index);
	dirc_model->set_three_seg_mirror(three_seg_mirror);
	dirc_model->set_sidemirror(sm_xr,sm_xl);

	dirc_model->set_pmt_offset(pmt_offset);
	dirc_model->set_upper_wedge_angle_diff(0.);
	dirc_model->set_bar_box_angle(bar_box_box_angle);
	dirc_model->set_bar_box_offsets(0., 0., 0.);
	dirc_model->set_focus_mirror_angle(main_mirror_angle, 0);
	dirc_model->set_upper_wedge_angle_diff(0., 0.);
 	dirc_model->set_bar_box_angle(0.);

	// conmpute and intialize the pdfs
	for (size_t particle = 0; particle < PARTICLE_NUMBER; ++particle) {
	    std::vector<dirc_point> hit_points;
	    std::back_insert_iterator<std::vector<dirc_point>> fill_hit_points = \
		std::back_inserter(hit_points);

	    //	    for (size_t kde_iteration = 0; kde_iteration < kde_generation_iterations;
	    //	 ++kde_iteration) 
	    {
		const float energy = energy_mean;
		//const float particle_eta = 0.5*(eta_max + eta_min);
		const float particle_eta = eta_mean;
		// degrees
		const float particle_theta = 90 - TMath::RadToDeg()*2*atan(exp(-particle_eta));
		const float beta = dirc_model->get_beta(energy, masses[particle]);
		// ns
		const float time = particle_flight_distance/(beta*.3);

		dirc_model->fill_reg_phi(fill_hit_points,
					 n_phi_phots,
					 n_z_phots,
					 -1,
					 1,
					 particle_x,
					 particle_y,
					 time,
					 particle_theta,
					 particle_phi,
					 0,
					 ckov_unc/pdf_unc_red_fac,
					 beta);
	    }
	    pdfs[particle] = std::make_unique<DircSpreadGaussian>(
	        sfunc_sig, hit_points, s_func_x, s_func_y, s_func_t);
	}
	
	printf("Beginning Run\n");
	for (unsigned int i = 0; i < num_runs; ++i) {
	    const int particle_one_n_sim_phots = spread_ang->Gaus(mean_n_phot, spread_n_phot);
	    // We want more or less the same number of
	    // signal particles of each type
	    particle_one_type = spread_ang->Integer(PARTICLE_NUMBER);
	    particle_one_energy = energy_mean; //spread_ang->Gaus(energy_mean, energy_spread);
	    particle_one_eta = eta_mean;//spread_ang->Uniform(eta_min, eta_max);
	    // degrees
	    const float particle_one_theta = 90 - TMath::RadToDeg()*2*atan(exp(-particle_one_eta));
	    const float particle_one_beta = dirc_model->get_beta(
	         particle_one_energy, masses[particle_one_type]);
	    // ns
	    const float particle_one_time = particle_flight_distance/(particle_one_beta*.3);
	    // assume its a middle bar
	    // The first particle is always (0, 0) - we have tracking!
	    std::vector<dirc_point> sim_points;
	    std::back_insert_iterator<std::vector<dirc_point>> fill_sim_points = \
		std::back_inserter(sim_points);
	    dirc_model->fill_rand_phi(fill_sim_points,
				      particle_one_n_sim_phots,
				      PARTICLE_ANGLE,
				      1,
				      0.,
				      0.,
				      particle_one_time,
				      particle_one_theta + const_track_off,
				      particle_phi,
				      tracking_unc,
				      ckov_unc,
				      particle_one_theta);
	    const float particle_two_n_sim_phots = spread_ang->Gaus(mean_n_phot, spread_n_phot);
	    const float particle_two_x = spread_ang->Gaus(particle_x_mean, particle_x_spread);
	    const float particle_two_y = spread_ang->Gaus(particle_y_mean, particle_y_spread);
	    // TODO(kazeevn) square?
	    tracks_distance = sqrt(particle_two_x*particle_two_x + particle_two_y*particle_two_y);
	    // const float particle_two_energy = spread_ang->Gaus(energy_mean, energy_spread);
	    // For the noise particle, we want an LHCb-like distribution
	    // Since TRandom3 doesn't provide weights, we use the
	    // standard library
	    particle_two_type = particle_type_generator(random_generator);
	    // const float particle_two_eta = spread_ang->Uniform(eta_min, eta_max);
	    // const float particle_two_theta = 90 - TMath::RadToDeg()*2*atan(exp(-particle_two_eta));
	    // const float particle_two_beta = dirc_model->get_beta(
	    //     particle_two_energy, masses[particle_two_type]);
	    // // ns
	    // const float particle_two_time = particle_flight_distance/(particle_two_beta*.3);
	    // dirc_model->fill_rand_phi(fill_sim_points,
	    // 			      particle_two_n_sim_phots,
	    // 			      PARTICLE_ANGLE,
	    // 			      1,
	    // 			      particle_two_x,
	    // 			      particle_two_y,
	    // 			      particle_two_time,
	    // 			      particle_two_theta + const_track_off,
	    // 			      particle_phi,
	    // 			      tracking_unc,
	    // 			      ckov_unc,
	    // 			      particle_two_beta);
	    digitizer.digitize_points(sim_points);
	    const float ll_pion = pdfs[ParticleTypes::Pion]->get_log_likelihood(sim_points);
	    for (size_t particle = 0; particle < PARTICLE_NUMBER; ++particle) {
		if (particle == ParticleTypes::Pion) {
		    continue;
		}
		dlls[particle] = pdfs[particle]->get_log_likelihood(sim_points) - ll_pion;
	    }
	    tree->Fill();
	}
	printf("\nRun Completed\n");
	tfile->cd();
	tree->Write();
	tfile->Close();
	return 0;
}
