#include <vector>
#include <utility>
#include <TRandom3.h>

#include "dirc_point.h"

#ifndef DIRC_BASE_SIM
#define DIRC_BASE_SIM 
struct dirc_base_sim_tracking_step
{
	//position at start of step
	float x;
	float y;
	float z;

	//radians
	float sin_theta;
	float cos_theta;
	float sin_phi;
	float cos_phi;
};
class DircBaseSim
{
protected:
	float barLength;
	float barWidth;
	float barDepth;
	float windowThickness;

	float wedgeTop;
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

	float wedgeClosePlaneNx;
	float wedgeClosePlaneNy;
	float wedgeClosePlaneNz;
	float wedgeClosePlaneD;

	float upperWedgeClosePlaneNx;
	float upperWedgeClosePlaneNy;
	float upperWedgeClosePlaneNz;
	float upperWedgeClosePlaneD;
	float lowerWedgeExtensionZ;

	float upperWedgeGap;
	
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
	
	float box_angle_off_cval;
	float box_angle_off_sval;
	float bar_box_xoff;
	float bar_box_yoff;
	float bar_box_zoff;
	float quartzIndex;
	float quartzLiquidY;

	int wedge_bounces;
	int lastWallX;
	int wedgeBeforeInterface;

	float moliereP;
	bool useMoliere;

	float liquidIndex;
	float liquidAbsorbtion;
	bool use_liquid_n;
	bool use_quartz_n_for_liquid;
      
	std::vector<float> dist_traveled;
	std::vector<float> refraction_before;
	std::vector<float> refraction_after;

	bool store_bounces;
	std::vector<int> x_bounces;
	std::vector<int> z_bounces;
	std::vector<int> x_direct_bounces;
	std::vector<int> z_direct_bounces;
	std::vector<int> x_indirect_bounces;
	std::vector<int> z_indirect_bounces;

	
	float min_transmittance,max_transmittance,sep_transmittance;
	int num_transmittance;
	std::vector<float> quartz_transmittance;
	TRandom3 *rand_gen;

	bool midLineMode;
	int midLineWedgeWallFlip;

	bool upperWedgeAngleStore;
	std::vector<float> upper_wedge_incident;	
	
	void build_system();
	void spread_wedge_mirror();

	bool quartz_transmission_mc(float R, float lambda);
	bool absorbtion_mc(float dx, float dy);
	
	void rotate_2d(float &x, float &y, float cval, float sval);

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


	float generate_cos_moliere_angle(\
                float rad_length);

	void bar_box_interface(\
		float &x,\
		float &y,\
		float &z,\
		float &dx,\
		float &dy,\
		float &dz);

	void fill_moliere_tracking_steps(\
		std::vector<dirc_base_sim_tracking_step> &rsteps,\
		float &travel_distance,\
		float step_length,\
		float start_theta,\
		float start_phi,\
		float start_x,\
		float start_y,\
		float start_z);
		

	
	float get_quartz_n(float lambda);
	float get_liquid_n(float lambda);
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
	bool x_wedge_coerce_check(\
		float &x,\
		float &y,\
		float &z,\
		float &dx,\
		float &dy,\
		float &dz,\
		float dt);
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
	//yet another inlinable utility function (make a second with no return?)
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


	//The compiler should be inlining this without our help
	float sgn(float val);
	

	//inherit and change this function for your implementation
	virtual void warp_readout_box(\
		dirc_point &out_val,\
		int particle_bar,\
		float &mm_index,\
		float &x,\
		float &y,\
		float &z,\
		float &dx,\
		float &dy,\
		float &dz) {}
public:
	// Also inherit and change this - allows for custom quantum efficiency.
	// See source/dirc_threeseg_box_sim.cpp for sample implementation
	virtual const float get_cerenkov_angle_rand(
	        float beta, float additional_spread, float &wavelength) = 0;
	

	void set_store_bounces(bool isb);
	void fill_bounces_vecs(\
		std::vector<int> &fxbounces,\
		std::vector<int> &fzbounces,\
		std::vector<int> &fxdirbounces,\
		std::vector<int> &fzdirbounces,\
		std::vector<int> &fxindirbounces,\
		std::vector<int> &fzindirbounces);
	std::vector<float> get_dist_traveled();
	void set_upper_wedge_angle_store(bool istore);
	std::vector<float> get_upper_wedge_incident();
	void set_liquid_index(float li);
	void set_bar_box_angle(float ang);
	void set_bar_box_offsets(float x, float y, float z);
	void set_wedge_mirror_rand(float ispread);
	float get_beta(float E, float m);
	void set_upper_wedge_angle_diff(float rads, float radsy_y = 0);	
	float get_bar_offset(int bar);
	int get_bar_from_x(float x);
	
	void set_use_quartz_n_for_liquid(bool iu);

	void set_moliere_p(float ip);
	void set_use_moliere(bool ium);
	//Random seed chosen arbitrarily
	//default parameters correspond to babar dirc bars
	//default upper wedge top is for gluex implementation.  Set to 0 to remove upper wedge
	DircBaseSim(\
		int rand_seed = 4357,\
		float ibarLength = 4900,\
		float ibarWidth = 35,\
		float ibarDepth = 17,\
		float iupperWedgeTop = 178.6);
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
	bool track_single_photon(\
	        dirc_point &out_val,\
	        float emit_theta,\
	        float emit_phi,\
       		float particle_theta,\
	        float particle_phi,\
	        float particle_x,\
	        float particle_y,\
	        float particle_z,\
	        float particle_t,\
		int particle_bar);
	bool track_single_photon_beta(\
	        dirc_point &out_val,\
	        float particle_beta,\
	        float emit_phi,\
       		float particle_theta,\
	        float particle_phi,\
	        float particle_x,\
	        float particle_y,\
	        float particle_z,\
	        float particle_t,\
		int particle_bar);
	bool track_line_photon(\
	        dirc_point &out_val,\
	        float particle_beta,\
	        float emit_phi,\
       		float particle_theta,\
	        float particle_phi,\
	        float particle_x,\
	        float particle_y,\
	        float particle_z,\
	        float particle_t,\
		int particle_bar,\
		float z_at_top = 1);
	bool track_all_line_photons(\
                std::vector<dirc_point> &left_vals,\
                std::vector<dirc_point> &right_vals,\
                int points_per_side,\
                float emit_theta,\
                float particle_theta,\
                float particle_phi,\
                float particle_x,\
                float particle_y,\
                float particle_z,\
                float particle_t,\
                int particle_bar,\
                float z_at_top =1);
	void sim_lut_points(\
                std::vector<dirc_point> &ovals,\
                std::vector<float> &phis,\
                std::vector<float> &thetas,\
                int n_photons, \
                float particle_bar /*= 0*/);

	void test_from_wedge_top(\
                std::vector<dirc_point> &ovals,\
                int n_photons, \
                float particle_bar = 1, \
                float particle_x = 0, \
                float phot_theta = 0, \
                float phot_phi = 0,\
                float theta_unc = 0,\
                float phi_unc = 0,\
		float overall_theta = 0); 
};
#endif
