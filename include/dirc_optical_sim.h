#include <vector>
#include <utility>
#include <TRandom3.h>

#include "dirc_point.h"
//#include "dirc_optical_components.h"

#ifndef DIRC_OPTICAL_SIM
#define DIRC_OPTICAL_SIM 
class DircOpticalSim
{
private:
	float foc_r;
	float foc_mirror_size;
	float foc_rot;
	float foc_yrot;
	float foc_zrot;
	float sens_size;
	float sens_rot;
	
	float barLength;
	float barWidth;
	float barDepth;
	float wedgeWidthOff;
	float wedgeDepthOff;
	float wedgeFarAngle;
	float wedgeCloseAngle;
	float wedgeWidth;
	float wedgeDepthHigh;
	float wedgeHeight;
	float upperWedgeDepthHigh;
	float upperWedgeTop;
	float upperWedgeHeight;
	float upperWedgeBottom;

	float windowThickness;
	
	float wedgeClosePlaneNx;
	float wedgeClosePlaneNy;
	float wedgeClosePlaneNz;
	float wedgeClosePlaneD;

	float largePlanarMirrorNx;
	float largePlanarMirrorNy;
	float largePlanarMirrorNz;
	float largePlanarMirrorD;
	float largePlanarMirrorMinZ;
	float largePlanarMirrorMaxZ;
	float pmtPlaneMinZ;
	float pmtPlaneMaxZ;

	float upperWedgeClosePlaneNx;
	float upperWedgeClosePlaneNy;
	float upperWedgeClosePlaneNz;
	float upperWedgeClosePlaneD;
	float lowerWedgeExtensionZ;

	float upperWedgeGap;
	
	
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
	
	float sensPlaneYdistConversion;
	float sensPlaneZdistConversion;
	
	bool upperWedgeNonUniform;
	float upperWedgeNonUniformSpread;
	
	float wedgeFarPlaneNx;
	float wedgeFarPlaneNy;
	float wedgeFarPlaneNz;
	float wedgeFarPlaneD;
	
	float upperWedgeFarPlaneNx;
	float upperWedgeFarPlaneNy;
	float upperWedgeFarPlaneNz;
	float upperWedgeFarPlaneD;
	
	float upperWedgeFarZ;
	
	float boxCloseZ;
	float reflOff;
	
	float focMirrorY;
	float focMirrorZ;
	
	bool three_seg_mirror;
	
	float sidemirror_xr;
	float sidemirror_xl;
	float sidemirror_reflectivity;
	
	float quartzIndex;
	float liquidIndex;
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
	
	float min_QE,max_QE,sep_QE;
	int num_QE;
	std::vector<float> vals_QE;
	float min_transmittance,max_transmittance,sep_transmittance;
	int num_transmittance;
	std::vector<float> quartz_transmittance;
	
	TRandom3 *rand_gen;
	
	void rotate_2d(float &x, float &y, float cval, float sval);
	
	
	void build_system();
	void fill_sens_plane_vecs();
	void fill_threeseg_plane_vecs();
	void fill_foc_mirror_vecs();
	void sidemirror_reflect_points(std::vector<dirc_point> &points);
	void spread_wedge_mirror();
	bool quartz_transmission_mc(float R, float lambda);
	bool absorbtion_mc(float dx, float dy);
	

	void fill_rand_phi(\
		std::vector<dirc_point> &ovals,\
		int n_photons, \
		float ckov_theta /*= 47*/, \
        	float particle_bar /*=0*/, \
		float particle_x /*= 0*/, \
		float particle_y /*= 0*/, \
		float particle_t /*= 0*/, \
		float particle_theta /*= 0*/, \
		float particle_phi /*= 0*/,\
		float phi_theta_unc /*= .0015*57.3*/,\
		float ckov_theta_unc /* = .0055*57.3*/,\
		float beta /* = -1*/);
	void fill_reg_phi(\
		std::vector<dirc_point> &fill_points,\
		int n_photons_phi, \
		int n_photons_z,\
		float ckov_theta /*= 47*/, \
	        float particle_bar /*=0*/, \
		float particle_x /*= 0*/, \
		float particle_y /*= 0*/, \
		float particle_t /*= 0*/, \
		float particle_theta /*= 0*/, \
		float particle_phi /*= 0*/,\
		float phi_theta_unc, /*= 0*/
		float ckov_theta_unc /* = 0*/,\
		float beta /* = -1*/);
	
	float get_quartz_n(float lambda);
	bool optical_interface_z(\
		float n1,\
		float n2,\
		float &dx,\
		float &dy,\
		float &dz);
	float warp_ray(\
		float &x,\
		float &y,\
		float &z,\
		float &dx,\
		float &dy,\
		float &dz,\
		float cos_critical_angle);
	float warp_wedge(\
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
	bool x_wedge_coerce_check(\
		float &x,\
		float &y,\
		float &z,\
		float &dx,\
		float &dy,\
		float &dz,\
		float dt);
	float three_seg_reflect(\
		float &x,\
		float &y,\
		float &z,\
		float &dx,\
		float &dy,\
		float &dz);
	//Utility function - should definitely be inlined
	void plane_reflect(\
		float Nx,\
		float Ny,\
		float Nz,\
		float D,\
		float &x,\
		float &y,\
		float &z,\
		float &dx,\
		float &dy,\
		float &dz,\
		float &dt,\
		float offang = 0);
	//Utility function - should combine this with above for some speed
	float get_z_intercept(\
		float Nx,\
		float Ny,\
		float Nz,\
		float D,\
		float x,\
		float y,\
		float z,\
		float dx,\
		float dy,\
		float dz);
	//yet another inlinable utility function (make a second with not return?)
	float get_intercept_plane(\
		float Nx,\
		float Ny,\
		float Nz,\
		float D,\
		float &x,\
		float &y,\
		float &z,\
		float dx,\
		float dy,\
		float dz);
	//yet another inlinable utility function
	float cylindrical_reflect(\
		float &x,\
		float &y,\
		float &z,\
		float &dx,\
		float &dy,\
		float &dz);
	//The compiler should be inlining this without our help
	float sgn(float val);
	//Technically one of the "warped functions, but can and should be optimized somehow
	float warp_sens_plane(\
		dirc_point &fill_val,\
		float &x,\
		float &y,\
		float &z,\
		float &dx,\
		float &dy,\
		float &dz);
public:
	void set_focmirror_nonuniformity(float nonuni_deg);
	void set_foc_mirror_r(float ifoc_r);
	void set_sidemirror(float ixr, float ixl);
	void set_sidemirror_reflectivity(float isr);
	void sidemirror_reflect_point(dirc_point &ipt);
	void set_three_seg_mirror(bool itsm);
	void set_kaleidoscope_plot(bool ikp);
	void set_pmt_offset(float r);
	void set_liquid_absorbtion(float iabs);
	std::vector<float> get_dist_traveled();
	void set_store_traveled(bool sst = true);
	void set_liquid_index(float li);
	void set_bar_box_angle(float ang);
	void set_focus_mirror_angle(float ang,float yang = 0, float zang = 0);
	void set_pmt_angle(float ang);
	void set_wedge_mirror_rand(float ispread);
	void set_pmt_plane_zs(float imin, float imax);
	void set_large_mirror_zs(float imin, float imax);
	float get_cerenkov_angle_rand(float beta, float additional_spread, float &wavelength);
	float get_beta(float E, float m);
	void set_upper_wedge_angle_diff(float rads, float radsy_y = 0);	
	float get_bar_offset(int bar);
	int get_bar_from_x(float x);
	
	DircOpticalSim(\
		int rand_seed = 4357,\
		float ifoc_r = -1200, \
		float ifoc_mirror_size = 300.38, \
		float ifoc_rot = 74.11, \
		float isens_size = 600, \
		float isens_rot = 47.87);
	std::vector<std::pair<float,float> > get_refraction_rand_phi(\
		std::vector<float> &before_interface,\
		std::vector<float> &after_interface,\
		std::vector<float> &pmt_incidence,\
		int n_photons, \
		float ckov_theta = 47, \
		float particle_x = 0, \
		float particle_y = 0, \
		float particle_theta = 0, \
		float particle_phi = 0,\
		float phi_theta_unc = .0015*57.3,\
		float ckov_theta_unc = .0055*57.3,\
		float beta = -1);
	void sim_rand_n_photons(\
		std::vector<dirc_point> &out_points,\
		int n_photons,\
		float ckov_theta = 47, \
        	float particle_bar= 0, \
		float particle_x = 0, \
		float particle_y = 0, \
		float particle_t = 0, \
		float particle_theta = 0, \
		float particle_phi = 0,\
		float phi_theta_unc = .08594,\
		float ckov_theta_unc = .3151,\
		float beta = -1);	
	void sim_reg_n_photons(\
		std::vector<dirc_point> &out_points,\
		int n_photons_phi,\
		int n_photons_z,\
		float ckov_theta = 47, \
        	float particle_bar= 0, \
		float particle_x = 0, \
		float particle_y = 0, \
		float particle_t = 0, \
		float particle_theta = 0, \
		float particle_phi = 0,\
		float phi_theta_unc = 0,\
		float ckov_theta_unc = 0,\
		float beta = -1);	
	void test_from_wedge_top(\
                std::vector<dirc_point> &ovals,\
                int n_photons, \
                float particle_bar /*= 1*/, \
                float particle_x /*= 0*/, \
                float phot_theta /*= 0*/, \
                float phot_phi /*= 0*/,\
                float theta_unc, /*= 0*/
                float phi_unc /* = 0*/,\
		float overall_theta); 
};
#endif
