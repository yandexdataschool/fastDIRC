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
	double sens_r;
	double sens_subtend_angle;

	double sensCylMinZ;
	double sensCylY;
	double sensCylZ;

	double boxCloseZ;
	
	double sidemirror_xr;
	double sidemirror_xl;
	double sidemirror_reflectivity;
	
	double liquidIndex;
	double quartzLiquidY;
	
	double liquidAbsorbtion;
	std::vector<double> dist_traveled;
	bool store_traveled;
	bool kaleidoscope_plot;
	
	bool store_refraction;
	std::vector<double> refraction_before;
	std::vector<double> refraction_after;

	bool storeOpticalAngles;
	std::vector<double> focus_photon_angles;
	std::vector<double> side_photon_angles;
	std::vector<double> large_flat_photon_angles;
	
	double min_QE,max_QE,sep_QE;
	int num_QE;
	std::vector<double> vals_QE;


	double min_transmittance,max_transmittance,sep_transmittance;
	int num_transmittance;

	bool absorbtion_mc(double dx, double dy);
	void build_readout_box();
	void sidemirror_reflect_points(std::vector<dirc_point> &points);
	void spread_wedge_mirror();

	void warp_readout_box(\
		dirc_point &out_val,\
		int particle_bar,\
		double &mm_index,\
		double &x,\
		double &y,\
		double &z,\
		double &dx,\
		double &dy,\
		double &dz);
	//Technically one of the "warp" functions, but can and should be optimized somehow
	void warp_sens_cyl(\
		dirc_point &fill_val,\
		double &mm_index,\
		double &x,\
		double &y,\
		double &z,\
		double dx,\
		double dy,\
		double dz);
	double warp_box(\
                double &x,\
                double &y,\
                double &z,\
                double &dx,\
                double &dy,\
                double &dz);



public:
	const double get_cerenkov_angle_rand(double beta, double additional_spread, double &wavelength);
	double get_sens_r();
	double get_sens_subtend_angle();

	//Note many of these functions are implemented for convient interfacing sith dircfit
	//should not be required in the end
	void set_focmirror_nonuniformity(double nonuni_deg);
	void set_foc_mirror_r(double ifoc_r);
	void set_sidemirror(double ixr, double ixl);
	void set_sidemirror_reflectivity(double isr);
	void sidemirror_reflect_point(dirc_point &ipt);
	void set_three_seg_mirror(bool itsm);
	void set_pmt_offset(double r);
	void set_liquid_absorbtion(double iabs);
	std::vector<double> get_dist_traveled();
	void set_store_traveled(bool sst = true);
	void set_focus_mirror_angle(double ang,double yang = 0, double zang = 0);
	void set_pmt_angle(double ang);
	void set_pmt_plane_zs(double imin, double imax);
	void set_large_mirror_zs(double imin, double imax);
	
	void set_store_optical_angles(bool ibool);
	std::vector<double> get_focus_photon_angles();
	std::vector<double> get_side_photon_angles();
	std::vector<double> get_large_flat_photon_angles();
	DircBaBarSim(\
		int rand_seed=4357,\
                double isens_r=540.66, \
                double isens_subtend_angle=52.4, \
                double ibar_length=4900,\
                double ibar_width=35,\
                double ibar_depth=17.25,
                double iupper_wedge_top=178.6); 

};
#endif
