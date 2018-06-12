#include <vector>
#include <utility>
#include <TRandom3.h>

#include "dirc_base_sim.h"
#include "dirc_point.h"

#ifndef DIRC_THREESEGBOX_SIM
#define DIRC_THREESEGBOX_SIM 
class DircThreeSegBoxSim : public DircBaseSim
{
protected:
	float foc_r;
	float foc_mirror_size;
	float foc_rot;
	float foc_yrot;
	float foc_zrot;
	float sens_size;
	float sens_rot;

	float largePlanarMirrorNx;
	float largePlanarMirrorNy;
	float largePlanarMirrorNz;
	float largePlanarMirrorD;
	float largePlanarMirrorMinZ;
	float largePlanarMirrorMaxZ;
	float pmtPlaneMinZ;
	float pmtPlaneMaxZ;
	
	float focMirrorBottom;
	float focMirrorTop;
	float focMirrorZDim;
	//Multiseg?  probably not.  If it goes up again, use arrays
	float threeSeg1Nx,threeSeg1Ny,threeSeg1Nz,threeSeg1D;
	float threeSeg2Nx,threeSeg2Ny,threeSeg2Nz,threeSeg2D;
	float threeSeg3Nx,threeSeg3Ny,threeSeg3Nz,threeSeg3D;
	
	//Not yet used - implement to branch faster on the threeseg
	float threeSeg1_2dny,threeSeg1_2dnz,threeSeg1_2dd;
	float threeSeg2_2dny,threeSeg2_2dnz,threeSeg2_2dd;
	float threeSeg3_2dny,threeSeg3_2dnz,threeSeg3_2dd;
	
	float threeSeg1Y,threeSeg1Z;
	float threeSeg2Y,threeSeg2Z;
	float threeSeg3Y,threeSeg3Z;
	
	bool nonUniformFocMirror;
	float foc_mirror_nonuni;
	
	
	float sensPlaneNx;
	float sensPlaneNy;
	float sensPlaneNz;
	float sensPlaneD;
	float sensPlaneY;
	float sensPlaneZ;
	
	float unReflSensPlaneNx;
	float unReflSensPlaneNy;
	float unReflSensPlaneNz;
	float unReflSensPlaneD;
	float unReflSensPlaneY;
	float unReflSensPlaneZ;

        float focPlaneNx;
        float focPlaneNy;
        float focPlaneNz;
        float focPlaneD;
        float focPlaneMinZ;

	float focYoff;
	float focZoff;

	float sensPlaneYdistConversion;
	float sensPlaneZdistConversion;

	float boxCloseZ;
	float reflOff;
	float baseReflOff;
	
	float focMirrorY;
	float focMirrorZ;
	
	bool three_seg_mirror;
	
	float sidemirror_xr;
	float sidemirror_xl;
	float sidemirror_reflectivity;
	
	float quartzLiquidY;
	
	float box_angle_off_cval;
	float box_angle_off_sval;

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
	void fill_sens_plane_vecs();
	void fill_threeseg_plane_vecs();
	void fill_foc_mirror_vecs();
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
	float cylindrical_reflect(\
		float &x,\
		float &y,\
		float &z,\
		float &dx,\
		float &dy,\
		float &dz);
	float three_seg_reflect(\
                float &x,\
                float &y,\
                float &z,\
                float &dx,\
                float &dy,\
                float &dz); 


	//Technically one of the "warp" functions, but can and should be optimized somehow
	float warp_sens_plane(\
		dirc_point &fill_val,\
		float &x,\
		float &y,\
		float &z,\
		float &dx,\
		float &dy,\
		float &dz);
	float warp_box(\
                float &x,\
                float &y,\
                float &z,\
                float &dx,\
                float &dy,\
                float &dz);



public:
	const float get_cerenkov_angle_rand(float beta, float additional_spread, float &wavelength);
	
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
	void set_mirror_plane_offsets(float off_y, float off_z);


	
	void set_store_optical_angles(bool ibool);
	std::vector<float> get_focus_photon_angles();
	std::vector<float> get_side_photon_angles();
	std::vector<float> get_large_flat_photon_angles();
	DircThreeSegBoxSim(\
		int rand_seed = 4357,\
		float ifoc_r = -1200, \
		float ifoc_mirror_size = 288, \
		float ifoc_rot = 74.11, \
		float isens_size = 600, \
		float isens_rot = 47.87,\
                float ibar_length=4900,\
		float ibar_width=35,\
                float ibar_depth=17.25,\
                float iupper_wedge_top = 178.6);
};
#endif
