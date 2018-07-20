#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <iterator>
#include <utility>
#include <stdlib.h>
#include <math.h>

#include <TFile.h>
#include <TH2F.h>
#include <TMath.h>
#include <TTree.h>
#include <TRandom3.h>

#include "../include/dirc_threesegbox_sim.h"
#include "../include/dirc_spread_gaussian.h"
#include "../include/dirc_rect_digitizer.h"

// TODO(kazeevn) make this block a proper class
// the code in main relies on the particle types
// being from 0 to PARTICLE_NUMBER -1
const unsigned int PARTICLE_NUMBER = 5;
enum ParticleTypes {
    Electron = 0,
    Muon = 1,
    Pion = 2,
    Kaon = 3,
    Proton = 4,
};

const int PARTICLE_ANGLE = -1;

// GeV/c^2
const std::array<float, PARTICLE_NUMBER> masses {
    0.5109989461e-3, .1057, .1396, .4937, .9382720813};

// LHCb-esque
const std::array<unsigned int, PARTICLE_NUMBER> particle_frequencies {
    2, 5, 75, 15, 5};


// mm, BaBar
const float interaction_point_height = 810;
const float interaction_point_x = 0.5*35.;
const float interaction_point_y = 0.5*4900;

const inline float sq(const float x) {
    return x*x;
}


void compute_geometry(const float eta, const float x,
		      float* theta_degrees,
		      float* phi_degrees,
		      float* flight_distance, // mm
		      float* y // mm
		      ) {
    const float exp_eta = exp(-eta);
    const float exp_eta_2 = exp_eta * exp_eta;
    // std::cerr << "90 - Beam theta " << 90 - TMath::RadToDeg()*2*atan(exp_eta) << std::endl;
    const float cot_beam_theta = 0.5/(exp_eta_2 + 1) * (1/exp_eta - exp_eta*exp_eta_2);
    const float dX = interaction_point_x  - x;
    const float L = cot_beam_theta * sqrt(sq(interaction_point_height) + sq(dX));
    const float a_2 = sq(dX) + sq(L);
    // always positive, as designed
    *flight_distance = sqrt(a_2 + sq(interaction_point_height));
    const float a = copysign(sqrt(a_2), L);
    *theta_degrees = TMath::RadToDeg()*atan(a/interaction_point_height);
    *phi_degrees = 90 + TMath::RadToDeg()*asin(dX/a);
    *y = interaction_point_y + L;
}

int main(int nargs, char* argv[]) {  
	float energy_mean = 6.0;
	float energy_spread = 1.5;
	// For BaBar max_theta = 71 deg.
	// thus maximum eta is 1.8
	// and we'll be a bit conservative
	const float eta_min = -1.5;
	const float eta_max = 1.5;
	std::array<std::unique_ptr<DircSpreadGaussian>, PARTICLE_NUMBER> pdfs;
	std::mt19937 random_generator;
	std::discrete_distribution<> particle_type_generator(
            particle_frequencies.begin(), particle_frequencies.end());
	
	// We have BaBar!
	const float particle_x_min = 0;
	const float particle_x_max = 35;
	
	unsigned int num_runs = 1000;
	const unsigned int num_runs_with_params = 50;
	float mean_n_phot = 400000;
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

	const float s_func_x = 6;
	const float s_func_y = s_func_x;
	const float s_func_t = 1.0;
	const float sfunc_sig = 1;

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
		else if (strcmp(argv[i], "-tracking_unc") == 0)
		{
			i++;
			tracking_unc = atof(argv[i]);
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
	const float pdf_unc_red_fac = 1.;
	std::unique_ptr<TRandom3> spread_ang = std::make_unique<TRandom3>(rseed + 3);
	random_generator.seed(rseed + 13442);
	auto dirc_model = std::make_unique<DircThreeSegBoxSim>(
			rseed,
			-1200 + mirror_r_difference,
			foc_mirror_size,
			main_mirror_angle,
			600,
			47.87 + box_rot + mirror_angle_change,
			4900,
			35,
                        17,
                        178.6);
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
	UChar_t particle_one_type, particle_two_type;
	tree->Branch("particle_one_type", &particle_one_type, 
		     "Type of the signal particle which is being measured/b");
	tree->Branch("particle_two_type", &particle_two_type, 
		     "Type of the noise particle/b");
	Float_t particle_one_energy, particle_two_energy, particle_one_eta, particle_two_eta;
	tree->Branch("particle_one_energy", &particle_one_energy,
		     "Energy of the signal particle, GeV/F");
	tree->Branch("particle_two_energy", &particle_two_energy,
		     "Energy of the noise particle, GeV/F");
	tree->Branch("particle_one_eta", &particle_one_eta,
		     "Pseudorapidity of particle one/F");
	tree->Branch("particle_two_eta", &particle_two_eta,
		     "Pseudorapidity of particle two/F");
	Float_t particle_one_x, particle_two_x;
	tree->Branch("particle_one_x", &particle_one_x,
		     "Particle one x, degrees/F");
	tree->Branch("particle_two_x", &particle_two_x,
		     "Particle two x, degrees/F");
	std::array<Float_t, PARTICLE_NUMBER> dlls;
	tree->Branch("dll_electron", &(dlls[ParticleTypes::Electron]), "LL(electron) - LL(pion)/F");
	tree->Branch("dll_kaon", &(dlls[ParticleTypes::Kaon]), "LL(kaon) - LL(pion)/F");
	tree->Branch("dll_muon", &(dlls[ParticleTypes::Muon]), "LL(muon) - LL(pion)/F");
	tree->Branch("dll_proton", &(dlls[ParticleTypes::Proton]), "LL(proton) - LL(pion)/F");
	std::array<UInt_t, PARTICLE_NUMBER> kde_support_n_photons;
	tree->Branch("support_electron", &(kde_support_n_photons[ParticleTypes::Electron]),
		     "Number of photons in electron KDE support");
	tree->Branch("support_kaon", &(kde_support_n_photons[ParticleTypes::Kaon]),
		     "Number of photons in kaon KDE support");
	tree->Branch("support_muon", &(kde_support_n_photons[ParticleTypes::Muon]),
		     "Number of photons in muon KDE support");
	tree->Branch("support_proton", &(kde_support_n_photons[ParticleTypes::Proton]),
		     "Number of photons in proton KDE support");
	tree->Branch("support_pion", &(kde_support_n_photons[ParticleTypes::Pion]),
		     "Number of photons in pion KDE support");
	Float_t dirc_bt;
	tree->Branch("dll_bt", &dirc_bt, "LL(Below threshold) - LL(pion)/F");
	UInt_t particle_one_n_photons, total_n_photons;
	tree->Branch("particle_one_n_photons", &particle_one_n_photons,
		     "Number of photons reaching detector for the signal particle, before digitization");
	tree->Branch("total_n_photons", &total_n_photons,
		     "Number of photons reaching detector in total, after digitization");
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

	const DircSpreadGaussian pdf_bt(
	    sfunc_sig, std::vector<dirc_point>(), s_func_x, s_func_y, s_func_t);

	printf("Beginning Run\n");
	for (unsigned int i = 0; i < num_runs; ++i) {
	    if (i % 100 == 0) {
		std::cout << "Iteration: " << i << std::endl;
	    }
	    // We want more or less the same number of
	    // signal particles of each type
	    particle_one_energy = spread_ang->Gaus(energy_mean, energy_spread);
	    particle_one_eta = spread_ang->Uniform(eta_min, eta_max);
	    particle_one_x = spread_ang->Uniform(particle_x_min, particle_x_max);
	    float particle_one_theta, particle_one_phi;
	    float particle_one_flight_distance, particle_one_y;
	    compute_geometry(particle_one_eta, particle_one_x,
			     &particle_one_theta,
			     &particle_one_phi,
			     &particle_one_flight_distance,
			     &particle_one_y);
	    // std::cerr << particle_one_eta << " " << particle_one_x << " " <<  \
	    // 	particle_one_theta << " " << particle_one_phi << " " <<	\
	    // 	particle_one_flight_distance << " " << particle_one_y << std::endl;
	    for (size_t particle = 0; particle < PARTICLE_NUMBER; ++particle) {
		std::vector<dirc_point> hit_points;
		std::back_insert_iterator<std::vector<dirc_point>> fill_hit_points = \
		    std::back_inserter(hit_points);
		const float beta = dirc_model->get_beta(particle_one_energy, masses[particle]);
		// ns, flight distance in mm unlike the original fastDIRC
		const float time = 1e-3 * particle_one_flight_distance/(beta*.3);
		dirc_model->fill_reg_phi(fill_hit_points,
					 n_phi_phots,
					 n_z_phots,
					 PARTICLE_ANGLE,
					 1,
					 particle_one_x,
					 particle_one_y,
					 time,
					 particle_one_theta,
					 particle_one_phi,
					 0,
					 ckov_unc/pdf_unc_red_fac,
					 beta);

		pdfs[particle] = std::make_unique<DircSpreadGaussian>(
		    sfunc_sig, hit_points, s_func_x, s_func_y, s_func_t);
		kde_support_n_photons[particle] = hit_points.size();
	    }
	    for (unsigned int j = 0; j < num_runs_with_params; ++j) {
		std::vector<dirc_point> sim_points;
		std::back_insert_iterator<std::vector<dirc_point>> fill_sim_points = \
		    std::back_inserter(sim_points);

		particle_one_type = spread_ang->Integer(PARTICLE_NUMBER);
		const float particle_one_beta = dirc_model->get_beta(
		    particle_one_energy, masses[particle_one_type]);
		// ns, flight distance in mm unlike the original fastDIRC
		const float particle_one_time = 1e-3*particle_one_flight_distance / \
		    (particle_one_beta*.3);
		const int particle_one_n_sim_phots = spread_ang->Gaus(mean_n_phot, spread_n_phot);
		dirc_model->fill_rand_phi(fill_sim_points,
					  particle_one_n_sim_phots,
					  PARTICLE_ANGLE,
					  1,
					  particle_one_x,
					  particle_one_y,
					  particle_one_time,
					  particle_one_theta,
					  particle_one_phi,
					  tracking_unc,
					  ckov_unc,
					  particle_one_beta);
		particle_one_n_photons = sim_points.size();
		const float particle_two_n_sim_phots = spread_ang->Gaus(mean_n_phot, spread_n_phot);
		particle_two_energy = spread_ang->Gaus(energy_mean, energy_spread);
		// For the noise particle, we want an LHCb-like distribution
		// Since TRandom3 doesn't provide weights, we use the
		// standard library
		particle_two_type = particle_type_generator(random_generator);
		particle_two_eta = spread_ang->Uniform(eta_min, eta_max);
		particle_two_x = spread_ang->Uniform(particle_x_min, particle_x_max);
		float particle_two_theta, particle_two_phi;
		float particle_two_flight_distance, particle_two_y;
		compute_geometry(particle_two_eta, particle_two_x,
				 &particle_two_theta,
				 &particle_two_phi,
				 &particle_two_flight_distance,
				 &particle_two_y);
		const float particle_two_beta = dirc_model->get_beta(
		    particle_two_energy, masses[particle_two_type]);
		const float particle_two_time = 1e-3*particle_two_flight_distance / \
		    (particle_two_beta*.3);
		dirc_model->fill_rand_phi(fill_sim_points,
					  particle_two_n_sim_phots,
					  PARTICLE_ANGLE,
					  1,
					  particle_two_x,
					  particle_two_y,
					  particle_two_time,
					  particle_two_theta,
					  particle_two_phi,
					  tracking_unc,
					  ckov_unc,
					  particle_two_beta);
		digitizer.digitize_points(sim_points);
		total_n_photons = sim_points.size();
		// TODO(kazeevn) a better model
		// TODO(kazeevn) blend the models
		if (sim_points.size() == 0) {
		    for (size_t particle = 0; particle < PARTICLE_NUMBER; ++particle) {
			const float p = (float)pdfs[particle]->get_support_size() / \
	        	    (float)(n_phi_phots * n_z_phots);
			dlls[particle] = particle_one_n_sim_phots * log(1 - p);
		    }
		    for (size_t particle = 0; particle < PARTICLE_NUMBER; ++particle) {
			if (particle == ParticleTypes::Pion) {
			    continue;
			}
			dlls[particle] -= dlls[ParticleTypes::Pion];
		    }
		    dirc_bt = log(1.) - dlls[ParticleTypes::Pion];
		} else {
		    const float ll_pion = pdfs[ParticleTypes::Pion]->get_log_likelihood(sim_points);
		    for (size_t particle = 0; particle < PARTICLE_NUMBER; ++particle) {
			if (particle == ParticleTypes::Pion) {
			    continue;
			}
			dlls[particle] = pdfs[particle]->get_log_likelihood(sim_points) - ll_pion;
		    }
		    dirc_bt = pdf_bt.get_log_likelihood(sim_points) - ll_pion;
		}
		tree->Fill();
	    }
	}
	printf("\nRun Completed\n");
	tfile->cd();
	tree->Write();
	tfile->Close();
	return 0;
}
