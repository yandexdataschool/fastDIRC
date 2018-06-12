#include <vector>
#include <utility>
#include <TRandom3.h>

#include "dirc_base_sim.h"
#include "dirc_point.h"

#ifndef DIRC_BABAR_SIM
#define DIRC_BABAR_SIM 
class DircBaBarSim : public DircBaseSim
{
protected:
	float sens_r;
	float sens_subtend_angle;

	float sensCylMinZ;
	float sensCylY;
	float sensCylZ;

	float boxCloseZ;
	
	float sidemirror_xr;
	float sidemirror_xl;
	float sidemirror_reflectivity;
	
	float liquidIndex;
	float quartzLiquidY;
	
	float liquidAbsorbtion;
	std::vector<float> dist_traveled;
	bool store_traveled;
	bool kaleidoscope_plot;
	
	bool store_refraction;
	std::vector<float> refraction_before;
	std::vector<float> refraction_after;

	bool storeOpticalAngles;
	std::vector<float> focus_photon_angles;
	std::vector<float> side_photon_angles;
	std::vector<float> large_flat_photon_angles;
	
	float min_QE,max_QE,sep_QE;
	int num_QE;
	std::vector<float> vals_QE;


	float min_transmittance,max_transmittance,sep_transmittance;
	int num_transmittance;

	bool absorbtion_mc(float dx, float dy);
	void build_readout_box();
	void sidemirror_reflect_points(std::vector<dirc_point> &points);
	void spread_wedge_mirror();

	void warp_readout_box(\
		dirc_point &out_val,\
		int particle_bar,\
		float &mm_index,\
		float &x,\
		float &y,\
		float &z,\
		float &dx,\
		float &dy,\
		float &dz);
	//Technically one of the "warp" functions, but can and should be optimized somehow
	void warp_sens_cyl(\
		dirc_point &fill_val,\
		float &mm_index,\
		float &x,\
		float &y,\
		float &z,\
		float dx,\
		float dy,\
		float dz);
	float warp_box(\
                float &x,\
                float &y,\
                float &z,\
                float &dx,\
                float &dy,\
                float &dz);



public:
	const float get_cerenkov_angle_rand(float beta, float additional_spread, float &wavelength);
	float get_sens_r();
	float get_sens_subtend_angle();

	//Note many of these functions are implemented for convient interfacing sith dircfit
	//should not be required in the end
	void set_focmirror_nonuniformity(float nonuni_deg);
	void set_foc_mirror_r(float ifoc_r);
	void set_sidemirror(float ixr, float ixl);
	void set_sidemirror_reflectivity(float isr);
	void sidemirror_reflect_point(dirc_point &ipt);
	void set_three_seg_mirror(bool itsm);
	void set_pmt_offset(float r);
	void set_liquid_absorbtion(float iabs);
	std::vector<float> get_dist_traveled();
	void set_store_traveled(bool sst = true);
	void set_focus_mirror_angle(float ang,float yang = 0, float zang = 0);
	void set_pmt_angle(float ang);
	void set_pmt_plane_zs(float imin, float imax);
	void set_large_mirror_zs(float imin, float imax);
	
	void set_store_optical_angles(bool ibool);
	std::vector<float> get_focus_photon_angles();
	std::vector<float> get_side_photon_angles();
	std::vector<float> get_large_flat_photon_angles();
	DircBaBarSim(\
		int rand_seed=4357,\
                float isens_r=540.66, \
                float isens_subtend_angle=52.4, \
                float ibar_length=4900,\
                float ibar_width=35,\
                float ibar_depth=17.25,
                float iupper_wedge_top=178.6); 

};
#endif
