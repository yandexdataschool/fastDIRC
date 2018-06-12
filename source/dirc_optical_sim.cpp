#include <vector>

#include "../include/dirc_optical_sim.h"
#include "../include/dirc_point.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
//
#include <stdlib.h>
#include <math.h>

#include <TRandom3.h>



DircOpticalSim::DircOpticalSim(
		int rand_seed /*=4357*/,\
		float ifoc_r/*=540.66*/, \
		float ifoc_mirror_size/*=300.38*/, \
		float ifoc_rot/*=-74.11*/, \
		float isens_size/*=600*/, \
		float isens_rot/*=90*/) {
	foc_r = ifoc_r;
	foc_mirror_size = ifoc_mirror_size;
	foc_rot = ifoc_rot;
	foc_yrot = 0;
	sens_size = isens_size;
	sens_rot = isens_rot;

	float prism_height = 0;
	bool remove_prism = true;
	if (remove_prism == true)
	{
		prism_height = 20;
		upperWedgeGap = prism_height;
		printf("Removing the small wedge (prism) with a size of %12.02f mm\n",prism_height);
	}

	kaleidoscope_plot = false;
	barLength=4900;
//	barLength=100;
	barWidth=35;
//	barWidth=17.25;
	barDepth=17.25;
//	barDepth=5.25;
	wedgeWidthOff = 1.75;
	wedgeDepthOff = 10;
	wedgeFarAngle = .006*57.3;
	wedgeCloseAngle = 30;
	wedgeWidth=barWidth - wedgeWidthOff;
	wedgeDepthHigh = 79;
	wedgeHeight = 91;
//	upperWedgeDepthHigh = 130;
	upperWedgeTop = 178.6;
	upperWedgeHeight = 78;
	upperWedgeBottom = upperWedgeTop-upperWedgeHeight;

//	printf("window thickness calc: %12.04f\n",upperWedgeBottom-wedgeHeight);

	//not, this is about 4mm off from the original BaBar schematic - the angles and 
        upperWedgeDepthHigh = wedgeDepthHigh + (upperWedgeTop-wedgeHeight)*sin(wedgeCloseAngle/57.3);


	lowerWedgeExtensionZ = -wedgeDepthHigh\
			       - tan(wedgeCloseAngle/57.296)*(upperWedgeBottom - wedgeHeight);

	printf("Lower Wedge Top: %12.04f Upper Wedge Bottom: %12.04f\n",wedgeHeight, upperWedgeBottom);

	//Variables used for plane intersection
	wedgeClosePlaneNx = 0; //Shouldn't be needed
	wedgeClosePlaneNy = sin(wedgeCloseAngle/57.296);
	wedgeClosePlaneNz = cos(wedgeCloseAngle/57.296);

	upperWedgeClosePlaneNx = 0; //Shouldn't be needed
	upperWedgeClosePlaneNy = sin(wedgeCloseAngle/57.296);
	upperWedgeClosePlaneNz = cos(wedgeCloseAngle/57.296);

	upperWedgeNonUniform = false;
	upperWedgeNonUniformSpread = 0;

	//only used for checking collision
	largePlanarMirrorNx = 0; //Shouldn't be needed
	largePlanarMirrorNy = 1;
	largePlanarMirrorNz = 0;
	largePlanarMirrorD = upperWedgeTop + barLength/2;
	largePlanarMirrorMinZ = -559;
	largePlanarMirrorMaxZ = -130;

	pmtPlaneMinZ = -559;
	pmtPlaneMaxZ = -329;

/*
	pmtPlaneMinZ = -1000;
	pmtPlaneMaxZ = 1000;
	largePlanarMirrorMinZ = -1000;
	largePlanarMirrorMaxZ = 1000;
*/	
	wedgeClosePlaneD = barLength/2*wedgeClosePlaneNy - wedgeClosePlaneNz * (barDepth+wedgeDepthOff);

	upperWedgeClosePlaneD = (barLength/2 + upperWedgeBottom)*upperWedgeClosePlaneNy + upperWedgeClosePlaneNz*lowerWedgeExtensionZ;

	upperWedgeFarPlaneNx = 0; //Shouldn't be needed
	upperWedgeFarPlaneNy = 0;
	upperWedgeFarPlaneNz = -1;

	upperWedgeFarPlaneD = (barLength/2 + upperWedgeBottom)*upperWedgeFarPlaneNy;

	focMirrorBottom = 139 + upperWedgeTop + barLength/2;




	sidemirror_xl = -1000000;
	sidemirror_xr = 1000000;
	sidemirror_reflectivity = .9;


	//Take average
	quartzIndex = 1.47;
	liquidIndex = 1.47;
	quartzLiquidY = upperWedgeBottom;


	//Negative to make reflections easier
	wedgeFarPlaneNx = 0; //Shouldn't be needed
	wedgeFarPlaneNy = -sin(wedgeFarAngle/57.296);
	wedgeFarPlaneNz = -cos(wedgeFarAngle/57.296);

	wedgeFarPlaneD = barLength/2*wedgeFarPlaneNy;

	upperWedgeFarZ = 0;//Change based on geometry

	boxCloseZ = -614;

	reflOff = 9;

	three_seg_mirror = false;

	rand_gen = new TRandom3(rand_seed);

	box_angle_off_cval = 1;
	box_angle_off_sval = 0;


	num_QE = 31;
	min_QE = 300;
	max_QE = 600;
	sep_QE = (max_QE - min_QE)/(num_QE - 1);

	float t_QE[31] = {\
		0.016415, 0.074064, 0.141658, 0.184219, 0.20634,  0.217191, 0.223244,
	       0.222296, 0.215232, 0.206302, 0.195574, 0.183007, 0.169403, 0.155447,
	       0.140014, 0.127696, 0.115716, 0.104086, 0.092256, 0.084083, 0.075418,
	       0.067311, 0.060243, 0.053588, 0.047765, 0.04344,  0.037999, 0.034177,
	       0.030869, 0.027848, 0.024141
	};


	// Transmittance of quartz per 1m



	for (int i = 0; i < num_QE; i++) {
		vals_QE.push_back(t_QE[i]);
	}

	num_transmittance = 36;
	min_transmittance = 300;
	max_transmittance = 660;
	sep_transmittance = (max_transmittance - min_transmittance)/(num_transmittance - 1);

	float tmp_quartz_transmittance[36] = {\
		0.999572036,0.999544661,0.999515062,0.999483019,0.999448285,\
			0.999410586,0.999369611,0.999325013,0.999276402,0.999223336,\
			0.999165317,0.999101778,0.999032079,0.998955488,0.998871172,\
			0.998778177,0.99867541 ,0.998561611,0.998435332,0.998294892,\
			0.998138345,0.997963425,0.997767484,0.997547418,0.99729958 ,\
			0.99701966 ,0.99670255 ,0.996342167,0.995931242,0.995461041,\
			0.994921022,0.994298396,0.993577567,0.992739402,0.991760297,\
			0.990610945\
	};
	for (int i = 0; i < num_transmittance; i++) {
		quartz_transmittance.push_back(tmp_quartz_transmittance[i]);
	}
	set_focmirror_nonuniformity(0);
	fill_threeseg_plane_vecs();
	fill_foc_mirror_vecs();
	fill_sens_plane_vecs();
	build_system();
}
void DircOpticalSim::set_pmt_plane_zs(float imin, float imax)
{
	pmtPlaneMinZ = imin;
	pmtPlaneMaxZ = imax;
}
void DircOpticalSim::set_large_mirror_zs(float imin, float imax)
{
	largePlanarMirrorMinZ = imin;
	largePlanarMirrorMaxZ = imax;
}
void DircOpticalSim::fill_sens_plane_vecs() {
	float adjusted_sens_size = 312;

	sensPlaneNx = 0;
	sensPlaneNy = sin(sens_rot/57.3);
	sensPlaneNz = cos(sens_rot/57.3);

	sensPlaneY = -adjusted_sens_size*sin(sens_rot/57.3)/2-reflOff+barLength/2;
	sensPlaneZ = boxCloseZ + sens_size*cos(sens_rot/57.3)/2;

	sensPlaneD = sensPlaneNy*sensPlaneY + sensPlaneNz*sensPlaneZ;


	sensPlaneYdistConversion = 1/cos(sens_rot/57.3);
	sensPlaneZdistConversion = 1/sin(sens_rot/57.3);
}
void DircOpticalSim::set_sidemirror_reflectivity(float isr) {
	sidemirror_reflectivity = isr;
}
void DircOpticalSim::set_foc_mirror_r(float ifoc_r) {
	foc_r = ifoc_r;
	build_system();
}
void DircOpticalSim::fill_foc_mirror_vecs() {
	//is off by pi/2 to reduce rounding errors and such
	float foc_center_ang = foc_rot/57.3 + acos(-foc_mirror_size/(2*foc_r));
	focMirrorY = focMirrorBottom - foc_r*cos(foc_center_ang);
	focMirrorZ = foc_r*sin(foc_center_ang);
}
void DircOpticalSim::fill_threeseg_plane_vecs() {
	focMirrorTop = focMirrorBottom + foc_mirror_size*cos(foc_rot/57.3);
	focMirrorZDim = foc_mirror_size*sin(foc_rot/57.3);
	//If we ever go to more than 3 segments, use a loop

	float theta_m = foc_rot/57.3;//radians of rotation to go through
	float theta_c = fabs(2*asin(focMirrorZDim/(2*foc_r)));//angle subtended by the mirror
	float seg_h = fabs(2*foc_r*sin(theta_c/6));//length of each segment;

	//I had to do some geometry and algebra to get these numbers, but in hindsight, it's obvious.  Always that way.
	float theta_1 = theta_m - theta_c/3;
	float theta_2 = theta_m;
	float theta_3 = theta_m + theta_c/3;

	threeSeg1Y = focMirrorBottom;
	threeSeg1Z = 0;

	threeSeg2Y = threeSeg1Y + seg_h*cos(theta_1);
	threeSeg2Z = threeSeg1Z - seg_h*sin(theta_1);

	threeSeg3Y = threeSeg2Y + seg_h*cos(theta_2);
	threeSeg3Z = threeSeg2Z - seg_h*sin(theta_2);


	threeSeg1Nx = 0;
	threeSeg1Ny = sin(theta_1);
	threeSeg1Nz = cos(theta_1);
	rotate_2d(threeSeg1Nx,threeSeg1Nz,cos(foc_yrot/57.3),sin(foc_yrot/57.3));//I think this is a slightly wrong rotation if the mirrors are carved out of a solid block, but it should be good enough at small angles
	rotate_2d(threeSeg1Nx,threeSeg1Ny,cos(foc_zrot/57.3),sin(foc_zrot/57.3));//I think this is a slightly wrong rotation if the mirrors are carved out of a solid block, but it should be good enough at small angles
	threeSeg1D = threeSeg1Ny*threeSeg1Y + threeSeg1Nz*threeSeg1Z;//Use point x=0 as reference

	threeSeg2Nx = 0;
	threeSeg2Ny = sin(theta_2);
	threeSeg2Nz = cos(theta_2);
	rotate_2d(threeSeg2Nx,threeSeg2Nz,cos(foc_yrot/57.3),sin(foc_yrot/57.3));
	rotate_2d(threeSeg2Nx,threeSeg2Ny,cos(foc_zrot/57.3),sin(foc_zrot/57.3));//I think this is a slightly wrong rotation if the mirrors are carved out of a solid block, but it should be good enough at small angles
	threeSeg2D = threeSeg2Ny*threeSeg2Y + threeSeg2Nz*threeSeg2Z;//Use point x=0 as reference

	threeSeg3Nx = 0;
	threeSeg3Ny = sin(theta_3);
	threeSeg3Nz = cos(theta_3);
	rotate_2d(threeSeg3Nx,threeSeg3Nz,cos(foc_yrot/57.3),sin(foc_yrot/57.3));
	rotate_2d(threeSeg3Nx,threeSeg3Ny,cos(foc_zrot/57.3),sin(foc_zrot/57.3));//I think this is a slightly wrong rotation if the mirrors are carved out of a solid block, but it should be good enough at small angles
	threeSeg3D = threeSeg3Ny*threeSeg3Y + threeSeg3Nz*threeSeg3Z;//Use point x=0 as reference

}
void DircOpticalSim::set_focmirror_nonuniformity(float nonuni_deg) {
	foc_mirror_nonuni = nonuni_deg;
	nonUniformFocMirror = (fabs(nonuni_deg) > .001);
}
void DircOpticalSim::set_kaleidoscope_plot(bool ikp) {
	kaleidoscope_plot = ikp;
}
void DircOpticalSim::sidemirror_reflect_points(std::vector<dirc_point> &points) {
	float tmpx = 0;
	for (unsigned int i = 0; i < points.size(); i++) {
		tmpx = points[i].x;
		while (tmpx < sidemirror_xl || tmpx > sidemirror_xr) {
			if (tmpx < sidemirror_xl) {
				tmpx = 2*sidemirror_xl - tmpx;
			}
			if (tmpx > sidemirror_xr) {
				tmpx = 2*sidemirror_xr - tmpx;
			}
		}
		points[i].x = tmpx;
	}
}
void DircOpticalSim::sidemirror_reflect_point(dirc_point &point) {
	float tmpx = 0;
	tmpx = point.x;
	while (tmpx < sidemirror_xl || tmpx > sidemirror_xr) {
		//printf("%12.04f ",tmpx);
		if (rand_gen->Uniform(0,1) > sidemirror_reflectivity)
		{
			point.t = -1337;
			break;
		}
		//printf("%12.04f ",tmpx);
		if (tmpx < sidemirror_xl) {
			tmpx = 2*sidemirror_xl - tmpx;
		}
		if (tmpx > sidemirror_xr) {
			tmpx = 2*sidemirror_xr - tmpx;
		}
	}
	point.x = tmpx;
	//printf("%12.04f \n",tmpx);
}
void DircOpticalSim::set_sidemirror(float ixr, float ixl) {
	sidemirror_xr = ixr;
	sidemirror_xl = ixl;
}
void DircOpticalSim::set_three_seg_mirror(bool itsm) {
	if (three_seg_mirror != itsm) {
		three_seg_mirror = itsm;
		build_system();
	}
}
void DircOpticalSim::set_pmt_offset(float r) {
	//positive r increases path length
	boxCloseZ = -614 - r*cos(sens_rot/57.3);
	reflOff = 9 + r*sin(sens_rot/57.3);
	build_system();
}
void DircOpticalSim::set_liquid_absorbtion(float iabs) {
	liquidAbsorbtion = iabs;
}
std::vector<float> DircOpticalSim::get_dist_traveled() {
	return dist_traveled;
}
void DircOpticalSim::set_store_traveled(bool sst/* = true*/) {
	store_traveled = sst;
}
void DircOpticalSim::set_liquid_index(float li) {
	liquidIndex = li;
	printf("Liquid Index: %12.04f\n",liquidIndex);
}
void DircOpticalSim::rotate_2d(float &x, float &y, float cval, float sval) {
	//Standard rotatation allows precomputation of matrix elements
	//Only store one variable for speed
	//Sin should be actual sin (not negative sin)
	float tx = x;
	x = cval*tx - sval*y;
	y = sval*tx + cval*y;
}
void DircOpticalSim::set_bar_box_angle(float ang) {
	//expect degrees
	box_angle_off_cval = cos(ang/57.3);
	box_angle_off_sval = sin(ang/57.3);
}
void DircOpticalSim::set_upper_wedge_angle_diff(float rads, float rads_y) {
	upperWedgeClosePlaneNx = 0; //Shouldn't be needed
	upperWedgeClosePlaneNy = sin(wedgeCloseAngle/57.296 + rads);
	upperWedgeClosePlaneNz = cos(wedgeCloseAngle/57.296 + rads);

	rotate_2d(upperWedgeClosePlaneNx,upperWedgeClosePlaneNz,cos(rads_y),sin(rads_y));

	upperWedgeClosePlaneD = (barLength/2 + upperWedgeBottom)*upperWedgeClosePlaneNy + upperWedgeClosePlaneNz*lowerWedgeExtensionZ;

	upperWedgeFarPlaneNx = 0; //Shouldn't be needed
	upperWedgeFarPlaneNy = -sin(rads);
	upperWedgeFarPlaneNz = -cos(rads);

	upperWedgeFarPlaneD = upperWedgeBottom*upperWedgeFarPlaneNy;
}
void DircOpticalSim::set_focus_mirror_angle(float ang,float yang, float zang) {
	foc_rot = ang;
	foc_yrot = yang;
	foc_zrot = zang;
	build_system();
}
void DircOpticalSim::set_pmt_angle(float ang) {
	sens_rot = ang;
	build_system();
}
void DircOpticalSim::set_wedge_mirror_rand(float ispread) {
	if (ispread > 0) {
		upperWedgeNonUniform = true;
		upperWedgeNonUniformSpread = ispread;
	} else {
		upperWedgeNonUniform = false;
		upperWedgeNonUniformSpread = 0;
	}
	spread_wedge_mirror();
}
void DircOpticalSim::spread_wedge_mirror() {

	if (upperWedgeNonUniform == true) {
		float spread_ang = wedgeCloseAngle;
		spread_ang += rand_gen->Gaus(0,upperWedgeNonUniformSpread);
		upperWedgeClosePlaneNx = 0; //Shouldn't be needed
		upperWedgeClosePlaneNy = sin(spread_ang/57.296);
		upperWedgeClosePlaneNz = cos(spread_ang/57.296);
		upperWedgeClosePlaneD = (barLength/2 + upperWedgeBottom)*upperWedgeClosePlaneNy + upperWedgeClosePlaneNz*lowerWedgeExtensionZ;
	}
}
void DircOpticalSim::build_system() {
	fill_foc_mirror_vecs();
	spread_wedge_mirror();
	fill_threeseg_plane_vecs();
	fill_sens_plane_vecs();
}
float DircOpticalSim::get_quartz_n(float lambda) {
	float B1,B2,B3,C1,C2,C3;
	B1 = 0.6961663;             // B1
	B2 = 0.4079426;             // B2
	B3 = 0.8974794;             // B3
	C1 = 0.0046791;             // C1
	C2 = 0.0135121;             // C2
	C3 = 97.9340025;          // C3
	float lam2;
	float n_lam;
	lam2 = lambda*lambda/1000000;

	n_lam = sqrt(1 + B1*lam2/(lam2-C1) + B2*lam2/(lam2-C2) + B3*lam2/(lam2-C3));

	return n_lam;
}
float DircOpticalSim::get_cerenkov_angle_rand(float beta, float additional_spread, float &wavelength) {
	//May be slow enough to consider approximating in distribution generation
	float out_ang = 0;
	float tmp_lam = 0;
	float tmp_QE_val;
	float above_ind;
	int ind_QE;
	float n_lam;

	while (true) {
		tmp_lam = rand_gen->Uniform(min_QE,max_QE);
		wavelength = tmp_lam;

		//Ignoring that the QE points are right on multiples of 10.  Assuming they are for speed.
		//This may not be neccessary, but I doubt it matters.
		ind_QE = (tmp_lam - min_QE)/sep_QE;

		above_ind = tmp_lam - (min_QE + sep_QE*ind_QE);

		//Simple linear interp between values.  5th order poly fit looked like it worked too
		tmp_QE_val = vals_QE[ind_QE]*(sep_QE-above_ind)/sep_QE + vals_QE[ind_QE+1]*above_ind/sep_QE;

		//Max QE val is ~.23, this saves lot of loops
		if (rand_gen->Uniform(0,.25) > tmp_QE_val) continue;

		//Test emission distribution, second b/c it's a less stringent cut
		if (rand_gen->Uniform(0,1/(min_QE*min_QE)) > 1/(tmp_lam*tmp_lam)) continue;


		n_lam = get_quartz_n(tmp_lam);

		out_ang = 57.3*acos(1/(beta*n_lam));
		break;
	}

	out_ang += rand_gen->Gaus(0,additional_spread);

	return out_ang;
}

bool DircOpticalSim::quartz_transmission_mc(float R, float lambda) {
	int ind_transmittance;
	float above_ind;
	float tmp_transmittance;
	float tmp_transmittance_constant;
	float trans_prob;

	ind_transmittance = (lambda - min_transmittance)/sep_transmittance;

	above_ind = lambda - (min_transmittance + sep_transmittance*ind_transmittance);

	//Simple linear interp between values.  5th order poly fit looked like it worked too
	tmp_transmittance = quartz_transmittance[ind_transmittance]*(sep_transmittance-above_ind)/sep_transmittance\
			    + quartz_transmittance[ind_transmittance+1]*above_ind/sep_transmittance;

	tmp_transmittance_constant = log(tmp_transmittance);



	trans_prob = exp(R/1000.*tmp_transmittance_constant);


	return (rand_gen->Uniform(0,1) < trans_prob);

}
float DircOpticalSim::get_beta(float E, float m) {
	float gam = E/m;
	if (gam > 5) {
		//Approximate and possibly save time
		return 1-1/(2*gam*gam);
	} else {
		return sqrt(1-1/(gam*gam));
	}

}
float DircOpticalSim::get_bar_offset(int bar)
{
	return fabs(bar)/bar*((150-0.5*(barWidth))+bar*(barWidth+.015));
}
int DircOpticalSim::get_bar_from_x(float x)
{

	if (x < -150)
	{
		return (x+150)/(barWidth+.015) - 1;
	}
	else if (x > 150)
	{
		return  (x-150)/(barWidth+.015) + 1;
	}
	else
	{
		return 0;
	}

}
void DircOpticalSim::sim_rand_n_photons(\
		std::vector<dirc_point> &out_points,\
		int n_photons, \
		float ckov_theta /*= 47*/, \
		float particle_bar /*=0*/, \
		float particle_x /*= 0*/, \
		float particle_y /*= 0*/, \
		float particle_t /*=0*/,\
		float particle_theta /*= 0*/, \
		float particle_phi /*= 0*/,\
		float phi_theta_unc /*= .0015*57.3*/,\
		float ckov_theta_unc /* = .0055*57.3*/,\
		float beta /* = -1*/) {

	//     std::vector<dirc_point> out_points;
	out_points.clear();
	fill_rand_phi(\
			out_points,\
			n_photons,\
			ckov_theta,\
			particle_bar,\
			particle_x,\
			particle_y,\
			particle_t /*=0*/,\
			particle_theta,\
			particle_phi,\
			phi_theta_unc,\
			ckov_theta_unc,\
			beta);
	//Reflection inside loops
	//    sidemirror_reflect_points(out_points);


	//     return out_points;
}
void DircOpticalSim::sim_reg_n_photons(\
		std::vector<dirc_point> &out_points,\
		int n_photons_phi, \
		int n_photons_z,\
		float ckov_theta /*= 47*/, \
		float particle_bar /*=0*/, \
		float particle_x /*= 0*/, \
		float particle_y /*= 0*/, \
		float particle_t /*=0*/,\
		float particle_theta /*= 0*/, \
		float particle_phi /*= 0*/,\
		float phi_theta_unc /*= 0*/,\
		float ckov_theta_unc /* = 0*/,\
		float beta /* = -1*/) {
	//     std::vector<dirc_point> out_points;
	out_points.clear();
	fill_reg_phi(\
			out_points,\
			n_photons_phi,\
			n_photons_z,\
			ckov_theta,\
			particle_bar,\
			particle_x,\
			particle_y,\
			particle_t,\
			particle_theta,\
			particle_phi,\
			phi_theta_unc,\
			ckov_theta_unc,\
			beta);
	//Reflection inside loops
	//    sidemirror_reflect_points(out_points);
	//     return out_points;
}
void DircOpticalSim::test_from_wedge_top(\
		std::vector<dirc_point> &ovals,\
		int n_photons, \
		float particle_bar /*= 1*/, \
		float particle_x /*= 0*/, \
		float phot_theta /*= 0*/, \
		float phot_phi /*= 0*/,\
		float theta_unc, /*= 0*/
		float phi_unc /* = 0*/,\
		float overall_theta) {

	ovals.clear();
	//Note that Theta and Phi are defined along the bar axis, not as they are elsewhere
	//negative bar number is - x axis?
	int numPhots = n_photons;

	float sourcez = -barDepth/2;
	sourcez = -wedgeDepthHigh/2;
	sourcez = 0;
	float sourcey = barLength/2 + upperWedgeTop + .1;
	float sourcex = particle_x;

	float temit, randPhi;

	float x,y,z,dx,dy,dz;

	float mm_index = 0;

	for (int i = 0; i < numPhots; i++) {
		mm_index = 0;
		randPhi = phot_phi + rand_gen->Gaus(0,phi_unc);
		temit = phot_theta + rand_gen->Gaus(0,theta_unc);

		x = sourcex;
		y = sourcey;
		z = sourcez;
		//Note, the geometry change again
		dx = sin(temit/57.3)*sin(randPhi/57.3);
		dy = cos(temit/57.3);
		dz = sin(temit/57.3)*cos(randPhi/57.3);

		rotate_2d(dy,dz,cos(overall_theta/57.3),sin(overall_theta/57.3));

		optical_interface_z(quartzIndex,liquidIndex,dx,dz,dy);

		dirc_point out_val;
		warp_readout_box(out_val,particle_bar,mm_index,x,y,z,dx,dy,dz);
		
		if (out_val.t < 0)
		{
			continue;
		}	

		ovals.push_back(out_val);
	}
}
void DircOpticalSim::fill_rand_phi(\
		std::vector<dirc_point> &ovals,\
		int n_photons, \
		float ckov_theta /*= 47*/, \
		float particle_bar /*= 0*/, \
		float particle_x /*= 0*/, \
		float particle_y /*= 0*/, \
		float particle_t /*= 0*/, \
		float particle_theta /*= 0*/, \
		float particle_phi /*= 0*/,\
		float phi_theta_unc, /*= .0015*57.3*/
		float ckov_theta_unc /* = .0055*57.3*/,\
		float beta/* = -1*/) {
	// 	float sDepth = .95*barDepth;

	//cdd
	//negative bar number is - x axis?
	float emitAngle = ckov_theta;
	float particleTheta = particle_theta + rand_gen->Gaus(0,phi_theta_unc);
	float particlePhi = particle_phi + rand_gen->Gaus(0,phi_theta_unc);
	int numPhots = n_photons/cos(particle_theta/57.3);

	float sourceOff,randPhi;

	float temit, rand_add;
	float wavelength = 400;

	float x,y,z,dx,dy,dz;

	float cos_ptheta = cos(particleTheta/57.3);
	float sin_ptheta = sin(particleTheta/57.3);
	float cos_pphi = cos(particlePhi/57.3);
	float sin_pphi = sin(particlePhi/57.3);

	float mm_index = 0;

	int tmp_updown = 0;

	for (int i = 0; i < numPhots; i++) {
		randPhi = rand_gen->Uniform(0,2*3.14159265);
		sourceOff = -rand_gen->Uniform(0,barDepth);
		if (kaleidoscope_plot == true)
		{	
			sourceOff = -barDepth/2;
		}

		if (beta < 0) {
			rand_add = rand_gen->Gaus(0,ckov_theta_unc);
			temit = emitAngle + rand_add;
		} else {
			temit = get_cerenkov_angle_rand(beta,ckov_theta_unc,wavelength);
		}
		mm_index = (sourceOff - barDepth)*1.47;

		x = 0;
		y = 0;
		z = sourceOff;
/*             
		x = 0;
		y = 0;
		z = barDepth/2;
 
		temit = 47;
                randPhi = 2.7635;
*/
		dx = sin(temit/57.3)*cos(randPhi);
		dy = sin(temit/57.3)*sin(randPhi);
		dz = cos(temit/57.3);

		rotate_2d(z,y,cos_ptheta,sin_ptheta);
		rotate_2d(x,y,cos_pphi,sin_pphi);

		rotate_2d(dz,dy,cos_ptheta,sin_ptheta);
		rotate_2d(dy,dx,cos_pphi,sin_pphi);

		z -= barDepth;
		x += particle_x;
		y += particle_y;

		//photon is now defined as up or down
		if (dy > 0)
		{
			tmp_updown = 1;
		}
		else
		{//Assume that zero doesn't count
			tmp_updown = -1;
		}

		mm_index += warp_ray(\
				x,\
				y,\
				z,\
				dx,\
				dy,\
				dz,\
				sqrt(1-1/(1.47*1.47)));

		if (z > 0)
		{
			continue;
		}

		spread_wedge_mirror();

		mm_index += warp_wedge(\
				x,\
				y,\
				z,\
				dx,\
				dy,\
				dz);

		//account (quickly) for the bar box having a different angle than the readout
		rotate_2d(dy,dz,box_angle_off_cval,box_angle_off_sval);

		if (z > 0) {
			continue;
		}
		dirc_point out_val;
              
		warp_readout_box(out_val,particle_bar,mm_index,x,y,z,dx,dy,dz);
	        
		if (out_val.t < 0)
		{
			continue;
		}

		out_val.t += particle_t;
		out_val.updown = tmp_updown;
		ovals.push_back(out_val);
	}
}
void DircOpticalSim::warp_readout_box(
	dirc_point &out_val,\
	int particle_bar,\
	float &mm_index,\
	float &x,\
	float &y,\
	float &z,\
	float &dx,\
	float &dy,\
	float &dz)
{
	float c_mm_ns = 300;
	mm_index += warp_box(\
			x,\
			y,\
			z,\
			dx,\
			dy,\
			dz);

	if (z > 0) {
		out_val.t = -1337;
		return;
	}

	//check absorbtion
	if (!(absorbtion_mc(dx,dy))) {
		out_val.t = -1337;
		return;
	}

	mm_index += warp_sens_plane(\
			out_val,\
			x,\
			y,\
			z,\
			dx,\
			dy,\
			dz);

	if (z > 0) {
		out_val.t = -1337;
		return;
	}

	out_val.t = mm_index/(c_mm_ns);
	out_val.x += get_bar_offset(particle_bar);
	
	//Must reflect after offset....
	sidemirror_reflect_point(out_val);
}
std::vector<std::pair<float,float> > DircOpticalSim::get_refraction_rand_phi(\
		std::vector<float> &before_interface,\
		std::vector<float> &after_interface,\
		std::vector<float> &pmt_incidence,\
		int n_photons, \
		float ckov_theta /*= 47*/, \
		float particle_x /*= 0*/, \
		float particle_y /*= 0*/, \
		float particle_theta /*= 0*/, \
		float particle_phi /*= 0*/,\
		float phi_theta_unc, /*= .0015*57.3*/
		float ckov_theta_unc /* = .0055*57.3*/,\
		float beta/* = -1*/) {
	//returns theta versus cerenkov phi_theta_unc

	std::vector<std::pair<float,float> > rval;

	before_interface.clear();
	after_interface.clear();
	refraction_before.clear();
	refraction_after.clear();
	pmt_incidence.clear();
	store_refraction = true;

	float emitAngle = ckov_theta;
	float particleTheta = particle_theta + rand_gen->Gaus(0,phi_theta_unc);
	float particlePhi = particle_phi + rand_gen->Gaus(0,phi_theta_unc);

	int numPhots = n_photons/cos(particle_theta/57.3);

	// 	float sourcez = -sDepth;
	float sourcey = particle_y-barDepth*tan(particleTheta/57.3);
	float sourcex = particle_x;
	float tsy = sourcey;
	float tsx = sourcex;
	sourcey = tsy*cos(particlePhi/57.3)-tsx*sin(particlePhi/57.3);
	sourcex = tsy*sin(particlePhi/57.3)+tsx*cos(particlePhi/57.3);


	float sourceOff,randPhi;

	float temit, rand_add;
	float wavelength = 400;

	float x,y,z,dx,dy,dz;

	float cos_ptheta = cos(particle_theta/57.3);
	float sin_ptheta = sin(particle_theta/57.3);
	float cos_pphi = cos(particle_phi/57.3);
	float sin_pphi = sin(particle_phi/57.3);

	float mm_index = 0;
	float c_mm_ns = 300;

	for (int i = 0; i < numPhots; i++) {
		randPhi = rand_gen->Uniform(0,2*3.14159265);
		sourceOff = -rand_gen->Uniform(0,barDepth);

		if (beta < 0) {
			rand_add = rand_gen->Gaus(0,ckov_theta_unc);
			temit = emitAngle + rand_add;
		} else {
			temit = get_cerenkov_angle_rand(beta,ckov_theta_unc,wavelength);
		}
		mm_index = (sourceOff - barDepth)*1.47;

		x = 0;
		y = 0;
		z = sourceOff;

		dx = sin(temit/57.3)*cos(randPhi);
		dy = sin(temit/57.3)*sin(randPhi);
		dz = cos(temit/57.3);

		rotate_2d(z,y,cos_ptheta,sin_ptheta);
		rotate_2d(x,y,cos_pphi,sin_pphi);

		rotate_2d(dz,dy,cos_ptheta,sin_ptheta);
		rotate_2d(dy,dx,cos_pphi,sin_pphi);

		z -= barDepth;
		x += particle_x;
		y += particle_y;


		mm_index += warp_ray(\
				x,\
				y,\
				z,\
				dx,\
				dy,\
				dz,\
				asin(1/1.47));

		if (z > 0) {
			continue;
		}


		spread_wedge_mirror();

		mm_index += warp_wedge(\
				x,\
				y,\
				z,\
				dx,\
				dy,\
				dz);

		//account (quickly) for the bar box having a different angle than the readout
		rotate_2d(dy,dz,box_angle_off_cval,box_angle_off_sval);

		if (z > 0) {
			continue;
		}

		mm_index += warp_box(\
				x,\
				y,\
				z,\
				dx,\
				dy,\
				dz);

		if (z > 0) {
			continue;
		}

		dirc_point out_val;
		
		//check absorbtion
		//We have the distance now - should use the real version
		if (!(absorbtion_mc(dx,dy))) {
			continue;
		}

		mm_index += warp_sens_plane(\
				out_val,\
				x,\
				y,\
				z,\
				dx,\
				dy,\
				dz);

		if (z > 0) {
			continue;
		}
		sidemirror_reflect_point(out_val);	
		if (out_val.t < 0)
		{
			continue;
		}
	
//		printf("%12.04f\n",sensPlaneNx*dx + sensPlaneNy*dy + sensPlaneNz*dz);
	
		pmt_incidence.push_back(acos(fabs(sensPlaneNx*dx + sensPlaneNy*dy + sensPlaneNz*dz)));

		//should be threading time information into this soon
		out_val.t = mm_index/(c_mm_ns);
		std::pair<float, float> add_val;
		
		add_val.first = randPhi;
		add_val.second = refraction_before[refraction_before.size() - 1];
		rval.push_back(add_val);
	}

	for (unsigned int i = 0; i < refraction_before.size(); i++) {
		before_interface.push_back(refraction_before[i]);
		after_interface.push_back(refraction_after[i]);
	}

	store_refraction = false;
	refraction_before.clear();
	refraction_after.clear();

	return rval;
}

void DircOpticalSim::fill_reg_phi(\
		std::vector<dirc_point> &ovals,\
		int n_photons_phi, \
		int n_photons_z,\
		float ckov_theta /*= 47*/, \
		float particle_bar /*= 0*/,\
		float particle_x /*= 0*/, \
		float particle_y /*= 0*/, \
		float particle_t /*= 0*/, \
		float particle_theta /*= 0*/, \
		float particle_phi /*= 0*/,\
		float phi_theta_unc /*= 0*/,\
		float ckov_theta_unc /* = 0*/,\
		float beta /* = -1*/)
{
	float sDepth = .95*barDepth;
	float emitAngle = ckov_theta;
//	float particleTheta = particle_theta;
//	float particlePhi = particle_phi;

	int tmp_updown = 0;

	float sourceOff,regPhi;

	float temit;
	float rand_add;
	float wavelength = 400;

	float x,y,z,dx,dy,dz;

	float cos_ptheta = cos(particle_theta/57.3);
	float sin_ptheta = sin(particle_theta/57.3);
	float cos_pphi = cos(particle_phi/57.3);
	float sin_pphi = sin(particle_phi/57.3);

	float mm_index = 0;

	float sin_emit;
	float cos_emit;
	float sin_regphi;
	float cos_regphi;

	int adj_n_photons_phi = n_photons_phi/cos(particle_theta/57.3);

	for (int i = 0; i < n_photons_z; i++) {
		sourceOff = (i+.5)*sDepth/(n_photons_z);

		if (kaleidoscope_plot == true)
		{
			sourceOff = -sDepth/2;
		}

		for (int j = 0; j < adj_n_photons_phi; j++) {
			regPhi = j*2*3.14159265357/(adj_n_photons_phi);

			if (beta < 0) {
				rand_add = rand_gen->Gaus(0,ckov_theta_unc);
				temit = emitAngle + rand_add;
			} else {
				temit = get_cerenkov_angle_rand(beta,ckov_theta_unc,wavelength);
			}

			mm_index = (sourceOff - barDepth)*1.47;

			x = 0;
			y = 0;
			z = sourceOff;

			//save some time ~30ms per 400k
			//could compute sin even faster with a taylor series
			sin_emit = sin(temit/57.2957795131);
			cos_emit = sqrt(1-sin_emit*sin_emit);
			cos_regphi = cos(regPhi);
			sin_regphi = sgn(3.14159265359 - regPhi)*sqrt(1-cos_regphi*cos_regphi);

			dx = sin_emit*cos_regphi;
			dy = sin_emit*sin_regphi;
			dz = cos_emit;

			rotate_2d(z,y,cos_ptheta,sin_ptheta);
			rotate_2d(x,y,cos_pphi,sin_pphi);

			rotate_2d(dz,dy,cos_ptheta,sin_ptheta);
			rotate_2d(dy,dx,cos_pphi,sin_pphi);

			z -= barDepth;
			x += particle_x;
			y += particle_y;

			//photon is now defined as up or down
			if (dy > 0)
			{
				tmp_updown = 1;
			}
			else
			{//Assume that zero doesn't count
				tmp_updown = -1;
			}

			mm_index += warp_ray(\
					x,\
					y,\
					z,\
					dx,\
					dy,\
					dz,\
					sqrt(1-1/(1.47*1.47)));

			if (z > 0)
			{
				continue;
			}

			spread_wedge_mirror();

			mm_index += warp_wedge(\
					x,\
					y,\
					z,\
					dx,\
					dy,\
					dz);

			//account (quickly) for the bar box having a different angle than the readout
			rotate_2d(dy,dz,box_angle_off_cval,box_angle_off_sval);

			if (z > 0)
			{
				continue;
			}

			dirc_point out_val;
			warp_readout_box(out_val,particle_bar,mm_index,x,y,z,dx,dy,dz);
			
			if (out_val.t < 0)
			{
				continue;
			}

			out_val.t += particle_t;
			out_val.updown = tmp_updown;
			ovals.push_back(out_val);
		}
	}

}
float DircOpticalSim::warp_ray(\
		float &x,\
		float &y,\
		float &z,\
		float &dx,\
		float &dy,\
		float &dz,\
		float cos_critical_angle) {
	//Implemented to avoid total internal reflection computations
	//Expects ray to be "prerotated" by particle theta and phi - just propagates the vectors
	//Uses bar geometry.
	//Modifys the ray.
	//Be careful about x,y,and z - can't straight rotate that
	//returns distance traveled (mm) times index

	float rval = 0; //stores total distance;

	float delx, dely,delz;
	// 	float xzR;
	//Start the rays in the quartz for simulation reasons
	float grace_room = 0;

	if (dy > 0) {
		//Going up
		dely = barLength*(0.5-grace_room) - y;
	} else {
		//going down, flip dy
		dely = barLength*(1.5-grace_room) + y;
		dy = -dy;
	}

	rval += dely*dely;

	y = barLength*(.5-grace_room);

	if (fabs(dz) > cos_critical_angle ||\
			fabs(dx) > cos_critical_angle ||\
			(dy < 0 && fabs(dy) > cos_critical_angle) ||
			dy*dy < 1e-4) {
		//If it's not totally internally reflected, assume it escapes
		//also assume failure if it isn't going up fast enough (absorbed)
		//Positive z origin signals fail
		z = 1337;
		return -1;
	}

	//Del on x and z refers to distance from first reflection
	//I sincerly hope dz is positive - might be worth checking

	int nbouncesx;
	// 	int nbouncesy;
	int nbouncesz;

	float remainderx;
	// 	float remaindery;
	float lrx;
	// 	float lrz = 1;
	float remainderz = 0;
	/*Not used right now, so don't branch if you don't have to
	  if(dy > 0)
	  {
	  nbouncesy = 0;
	  }
	  else
	  {
	  nbouncesy = 1;
	  }*/

	//deterimines if particle is going left or right

	lrx = sgn(dx);

	delx = fabs(dely*dx/dy) - lrx*(barWidth/2-x);
	delz = fabs(dely*dz/dy) + z;

	rval += delx*delx + delz*delz;

	// 	if (dz < 0) printf("dz: %12.04f : Negative and confusing\n",dz);


	nbouncesx = delx/barWidth + 1;
	remainderx = delx - barWidth*(nbouncesx-1);

	if (nbouncesx % 2 == 1) {
		dx = -dx;
		x = lrx*(barWidth/2-remainderx);
	} else {
		x = lrx*(-barWidth/2+remainderx);
	}

	if (x > barWidth/2 - wedgeWidthOff) {
		//Hit the 5mm where the wedge is not.  Assume lost
		z = 1337;
		return -1;
	}
	nbouncesz = delz/barDepth + 1;
	remainderz = delz - barDepth*(nbouncesz-1);

	if (nbouncesz % 2 == 1) {
		dz = -dz;
		z = -remainderz;
	} else {
		z = remainderz-barDepth;
	}

	return sqrt(rval)*quartzIndex;
}
float DircOpticalSim::warp_wedge(\
		float &x,\
		float &y,\
		float &z,\
		float &dx,\
		float &dy,\
		float &dz) {
	//No Critical angle checking - shouldn't need it
	//I am abusing the x angle being the same in the bar and wedge here
	//may have to rotate after to the mirror
	//Starts y at 0

	float mm_index = 0;
	bool passed_interface = false;
	// 	float x0 = x;
	float dt = 0; //dummy variable representing parametric line
	float n_dot_v = 0; //dummy variable to save dot product
	float n_dot_v0 = 0;

	//deal with yz first

	//Check for reflection from far wedge plane - max 1 bounce
	if (dz > 0) {
		n_dot_v = -(dy*wedgeFarPlaneNy + dz*wedgeFarPlaneNz);
		n_dot_v0 = -(y*wedgeFarPlaneNy + z*wedgeFarPlaneNz);

		dt = -(wedgeFarPlaneD+n_dot_v0)/n_dot_v;
		//pretending y=0 at bar top for speed
		if (dt*dy < wedgeHeight) {
			//reflect it off bottom wedge
			//Does not pass through optical interface

			//Will always be true - just use for propagation (short circuit var?)

			mm_index += dt*quartzIndex;
			x_wedge_coerce_check(x,y,z,dx,dy,dz,dt);


			dy += 2*n_dot_v*wedgeFarPlaneNy;
			dz += 2*n_dot_v*wedgeFarPlaneNz;
		} else {
			n_dot_v = -(dy*upperWedgeFarPlaneNy + dz*upperWedgeFarPlaneNz);
			n_dot_v0 = -(y*upperWedgeFarPlaneNy + z*upperWedgeFarPlaneNz);

			dt = -(upperWedgeFarPlaneD+n_dot_v0)/n_dot_v;
			//never reached??? probably not  can possibly remove if statement
			//No bottom wedge, try top
			if (dt*dy < upperWedgeTop) { //Should always be true... I hope (remove later?)
				//Does pass through optical interface

				//Following statement performs the propagation if it does not fail
				mm_index += dt*quartzIndex;
				if (!(x_wedge_coerce_check(x,y,z,dx,dy,dz,dt))) {
					//Goes out the window, fail and return
					z = 1337;
					return -1;
				}

				//Reflect off top wedge
				dy += 2*n_dot_v*upperWedgeFarPlaneNy;
				dz += 2*n_dot_v*upperWedgeFarPlaneNz;
			}
		}
	}
	//Now dz < 0 or we have a new starting vector.  Either way, we intersect with the "close" wedge now
	n_dot_v = -(dy*wedgeClosePlaneNy + dz*wedgeClosePlaneNz);
	n_dot_v0 = -(y*wedgeClosePlaneNy + z*wedgeClosePlaneNz);

	dt = -(wedgeClosePlaneD+n_dot_v0)/n_dot_v;

	//Assume close enough to determine which wedge it hit
	float ty = dt*dy + y - barLength/2;

	if ((ty < wedgeHeight)) {
		//reflect it off bottom wedge

		//Again, always true
		mm_index += dt*quartzIndex;
		if (!(x_wedge_coerce_check(x,y,z,dx,dy,dz,dt))) {
			//Goes out the window, fail and return
			//Shouldn't
			z = 1337;
			return -1;
		}

		dy += 2*n_dot_v*wedgeClosePlaneNy;
		dz += 2*n_dot_v*wedgeClosePlaneNz;

	} else if (ty > upperWedgeBottom && ty < upperWedgeTop) {
		//passed interface before reflection
		//get dt to propagate to interface
		dt = (quartzLiquidY + barLength/2 - y)/dy;
		mm_index += dt*quartzIndex;
		if (!(x_wedge_coerce_check(x,y,z,dx,dy,dz,dt))) {
			//Goes out the window, fail and return
			z = 1337;
			return -1;
		}
		//correct ordering below - interface is xz plane, so dx and dz go first
		optical_interface_z(quartzIndex,liquidIndex,dx,dz,dy);
		passed_interface = true;

		//NOW intersect with upper wedge

		//Signs are weird here, it works, but I think they are twisted up with the initialization too
		n_dot_v = -(dy*upperWedgeClosePlaneNy + dz*upperWedgeClosePlaneNz + dx*upperWedgeClosePlaneNx);
		n_dot_v0 = -(y*upperWedgeClosePlaneNy + z*upperWedgeClosePlaneNz + x*upperWedgeClosePlaneNx);


		dt = -(upperWedgeClosePlaneD+n_dot_v0)/n_dot_v;

		if (dt*dy + y < upperWedgeTop + barLength/2) {
			mm_index += dt*liquidIndex;
			if (!(x_wedge_coerce_check(x,y,z,dx,dy,dz,dt))) {
				//Goes out the window, fail and return
				z = 1337;
				return -1;
			}

			if (y > wedgeHeight + barLength/2 && y < upperWedgeBottom + upperWedgeGap + barLength/2)
			{
				//We let it go out the gap on the slanted side, but not on the left and right.
				z = 1337;
				return -1;
			}

			dx += 2*n_dot_v*upperWedgeClosePlaneNx;
			dy += 2*n_dot_v*upperWedgeClosePlaneNy;
			dz += 2*n_dot_v*upperWedgeClosePlaneNz;

		} else {
			//refracted such that it no longer bounces
			dt = (upperWedgeTop + barLength/2 - y)/dy;
			mm_index += dt*liquidIndex;
			if (!(x_wedge_coerce_check(x,y,z,dx,dy,dz,dt))) {
				//Goes out the window, fail and return
				z = 1337;
				return -1;
			}
		}

	} else if (ty < upperWedgeTop) {
		//out the window
		z = 1337;
		return -1;
	}

	//Have now performed the maximum number of reflections (besides x direction)
	//Finish by taking it to the top
	if (!passed_interface == true) {
		//Will go through the interface now before it hits the top
		dt = (quartzLiquidY + barLength/2 - y)/dy;
		mm_index += dt*quartzIndex;
		if (!(x_wedge_coerce_check(x,y,z,dx,dy,dz,dt))) {
			//Goes out the window, fail and return
			z = 1337;
			return -1;
		}
		//correct ordering below - interface is xz plane, so dx and dz go first
		optical_interface_z(quartzIndex,liquidIndex,dx,dz,dy);
		passed_interface = true;
	}

	//Now we just finish
//	printf("liquidIndex: %12.04f  quartzIndex: %12.04f\n",liquidIndex,quartzIndex);
//      printf("%12.04f %12.04f %12.04f %12.04f %12.04f %12.04f\n",x,y,z,dx,dy,dz);
	dt = (upperWedgeTop + barLength/2 - y)/dy;
	mm_index += dt*liquidIndex;
	if (!(x_wedge_coerce_check(x,y,z,dx,dy,dz,dt))) {
		//I'm not sure if it can go out the window at this point, but just in case
		//Goes out the window, fail and return
		z = 1337;
		return -1;
	}


	//and we're done.  Everything should be correct now with the ray at the top of the wedge
	//return the distance traveled times the index it traveled in

	return mm_index;

}
//Possibly inline these or something for speed, but right now, leave them for sanity
bool DircOpticalSim::optical_interface_z(\
		float n1,\
		float n2,\
		float &dx,\
		float &dy,\
		float &dz) {
	//n1 is starting index of refraction, n2 is ending
	//Assume that it's going through the z plane
	//(this will change the ordering when called in our current geometry)

	//dz and dy are flipped, so this is really acos
 	float before_ang = acos(dz);
//	printf("%12.04f\n",before_ang);

	float n12rat = n1/n2;
	float n12rat2 = n12rat*n12rat;

	//takes care of trig
	dz = sqrt(1-n12rat2*(1-dz*dz));

	// 	if (dz != dz) return false;//total internal reflection

	//simpler expression than I expected
	dx *= n12rat;
	dy *= n12rat;

// 	printf("in optical interface store_refraction: %s\n",store_refraction ? "true" : "false");
 	if (store_refraction == true)
 	{
 		float after_ang = acos(dz);
 	  	refraction_before.push_back(before_ang);
 	  	refraction_after.push_back(after_ang);
 	}

	return true;
}
bool DircOpticalSim::x_wedge_coerce_check(\
		float &x,\
		float &y,\
		float &z,\
		float &dx,\
		float &dy,\
		float &dz,\
		float dt) {
	//already have dt, be sure to add it in when calling this

	//assumes x starts between the sides of the wedge
	//assumes dy/dx constant over the whole distance
	//also performs the index plane crossing

	//changes x and dx, correctly propagating them
	//propagates y and z as well
	//returns false if it goes out the window on the x step
	//Will not check for y and z directions (assumes dt corresponds to an intersection/reflection point)

	float cdt = 0;
	float tdt = 0;
	float ty = y - barLength/2;


	while (true) {
		//Get time of this bounce based on x direction
		if (dx < 0) {
			tdt = x + barWidth/2;
			tdt /= -dx;
		} else {
			tdt = -x + barWidth/2 - wedgeWidthOff;
			tdt /= dx;
		}
		if (tdt + cdt < dt) {
			//bounced
			//add to total time taken
			cdt += tdt;
			ty += dy*tdt;
			if (ty > wedgeHeight && ty < upperWedgeBottom) {
				//out the window on an edge
				//this is ok, finish propagation without bouncing
				//return false;
				x += tdt*dx;
				tdt = dt - cdt;
				break;
			}

			//not out the window, propagate x for the next iteration
			x += tdt*dx;
			if (ty < upperWedgeBottom) {
				//Must be less than Wedge height at this point
				//bounce off lower wedge
				dx = -dx;
			} else {
				//if in upper wedge, don't bounce and propagate like normal
				tdt = dt - cdt;
				break;
			}
		} else {
			//does not bounce - record remaining dt and break;
			tdt = dt - cdt;
			break;
		}
	}
	//Finish last leg of trip:
	x += dx*tdt;
	y = barLength/2 + ty + dy*tdt;//Note that y starts at zero (bar middle), but ty starts at bar top
	z += dt*dz;//Implemented for convenience - z direction irrelevant to window calculation

	return true;//not out the window for the whole trip
}
bool DircOpticalSim::absorbtion_mc(float dx, float dy) {
	//True if the particle makes it
	//expects vector pointing after bounce
	//Magic number corresponding to minimal distancel traveled
	float min_dist = 650+56;

	//approximating here - assume yz distance is independent of where on mirror it hit
	//measuring side plots with a ruler seems to mean this is good to ~5-10%
	float approx_dist = min_dist*sqrt(dx*dx+dy*dy)/dy;

	if (store_traveled == true) {
		dist_traveled.push_back(approx_dist);
	}

	float prob_transmitted = exp(-liquidAbsorbtion*approx_dist);
	if (rand_gen->Uniform(0,1) < prob_transmitted) {
		return true;
	} else {
		return false;
	}
}
float DircOpticalSim::warp_box(\
		float &x,\
		float &y,\
		float &z,\
		float &dx,\
		float &dy,\
		float &dz) {
	//propagates light from the top of the wedge until just after bouncing off of the focusing mirrors
	//Consider branch profiling or some such other optimizations
	//Also, Fast Math, but don't tell the people in the penthouse

	//does not currently include x bouncing off the side - need to add that later
	//possibly in the last "warp to plane" bit
	float rval = 0; //distance traveled

	//first reflect off the back of the bar
	float dt;

	if (dz > 0) {
		//Condition may not be needed - try without for speed later
		dt = -z/dz;

		if (y + dy*dt < focMirrorBottom) {
			//reflects off of back
			x += dx*dt;
			y += dy*dt;
			z += dz*dt;

			rval += dt*liquidIndex;

			//Abuse the fact that this vector is normalized and hope it stays that way
			dz = -dz;
			//rval += dt;
		}
	}


	float trval  = 0;

	if (three_seg_mirror == true) {
		trval = three_seg_reflect(\
				x,\
				y,\
				z,\
				dx,\
				dy,\
				dz);
	} else {
		trval = cylindrical_reflect(\
				x,\
				y,\
				z,\
				dx,\
				dy,\
				dz);
	}

	if (trval < 0) {
		z = 1337;
		return -1;
	}

	//short propagate in front of the mirrors:
	//TODO remove later for speed
	x += .1*dx;
	y += .1*dy;
	z += .1*dz;


	rval += trval;

	return rval;
}
void DircOpticalSim::plane_reflect(\
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
		float offang /*=0*/) {
	//Propagate to the intercept and reflects off a plane
	//Could use this in the Wedge, but would loose some speed

	float n_dot_v = (Nx*dx + Ny*dy + Nz*dz);
	float n_dot_v0 = (Nx*x + Ny*y + Nz*z);

	dt = (D - n_dot_v0)/n_dot_v;

	x += dt*dx;
	y += dt*dy;
	z += dt*dz;

	if (fabs(offang) > 1e-8) {
		//remove for speed obviously
		rotate_2d(Nz,Ny,cos(offang),sin(offang));
	}

	dx -= 2*n_dot_v*Nx;
	dy -= 2*n_dot_v*Ny;
	dz -= 2*n_dot_v*Nz;
}
float DircOpticalSim::get_z_intercept(\
		float Nx,\
		float Ny,\
		float Nz,\
		float D,\
		float x,\
		float y,\
		float z,\
		float dx,\
		float dy,\
		float dz) {
	//finds z intercept (without modifying x y and z)
	//Could use this in the Wedge, but would loose some speed

	float n_dot_v = (Nx*dx + Ny*dy + Nz*dz);
	float n_dot_v0 = (Nx*x + Ny*y + Nz*z);
	float dt;//optimized by math compression?
	dt = (D - n_dot_v0)/n_dot_v;

	return z + dt*dz;
}
float DircOpticalSim::get_intercept_plane(\
		float Nx,\
		float Ny,\
		float Nz,\
		float D,\
		float &x,\
		float &y,\
		float &z,\
		float dx,\
		float dy,\
		float dz) {
	//Just intersects a plane (modifying x,y,z)
	//Could use this in the Wedge, but would loose some speed
	//Non returning verion for speed when needed?

	float n_dot_v = (Nx*dx + Ny*dy + Nz*dz);
	float n_dot_v0 = (Nx*x + Ny*y + Nz*z);
	float dt;//optimized by math compression?
	dt = (D - n_dot_v0)/n_dot_v;

	x += dx*dt;
	y += dy*dt;
	z += dz*dt;

	return dt;
}
float DircOpticalSim::three_seg_reflect(\
		float &x,\
		float &y,\
		float &z,\
		float &dx,\
		float &dy,\
		float &dz) {
	//Intersects and reflects the three segment plane
	float rval = 0;
	//definitely losing time here - combuting the n_dot_v and n_dot_v0 and dt twice for the chosen plane
	//TODO fix this once the code is correct

	//I hope there's a fast way to do these reflections

	float tz = 0;

	//check first seg (again, loop this if we go more than 3)
	tz = get_z_intercept(\
			threeSeg1Nx,\
			threeSeg1Ny,\
			threeSeg1Nz,\
			threeSeg1D,\
			x,\
			y,\
			z,\
			dx,\
			dy,\
			dz);
	float offang = 0;
	if (tz > threeSeg2Z && tz < 0) {
		//reflect off mirror closest to box
		if (nonUniformFocMirror == true) {
			//obviously not in the real run
			offang = rand_gen->Gaus(0,foc_mirror_nonuni/57.3);
			// 			rotate_2d(threeSeg1Nz,threeSeg1Ny,cos(offang),sin(offang));
		}
		plane_reflect(\
				threeSeg1Nx,\
				threeSeg1Ny,\
				threeSeg1Nz,\
				threeSeg1D,\
				x,\
				y,\
				z,\
				dx,\
				dy,\
				dz,\
				rval,\
				offang);
		return rval*liquidIndex;
	}

	tz = get_z_intercept(\
			threeSeg2Nx,\
			threeSeg2Ny,\
			threeSeg2Nz,\
			threeSeg2D,\
			x,\
			y,\
			z,\
			dx,\
			dy,\
			dz);
	// 	printf("tz2: %12.04f\n",tz);
	if (tz > threeSeg3Z && tz < 0) {
		//reflect off middle mirror
		if (nonUniformFocMirror == true) {
			//obviously not in the real run
			offang = rand_gen->Gaus(0,foc_mirror_nonuni/57.3);
			// 			rotate_2d(threeSeg2Nz,threeSeg2Ny,cos(offang),sin(offang));
		}
		plane_reflect(\
				threeSeg2Nx,\
				threeSeg2Ny,\
				threeSeg2Nz,\
				threeSeg2D,\
				x,\
				y,\
				z,\
				dx,\
				dy,\
				dz,\
				rval,\
				offang);
		return rval*liquidIndex;
	}

	tz = get_z_intercept(\
			threeSeg3Nx,\
			threeSeg3Ny,\
			threeSeg3Nz,\
			threeSeg3D,\
			x,\
			y,\
			z,\
			dx,\
			dy,\
			dz);
	// 	printf("tz3: %12.04f\n",tz);
	if (tz > -focMirrorZDim && tz < 0) {
		//reflect off mirror closest to box
		if (nonUniformFocMirror == true) {
			//obviously not in the real run
			offang = rand_gen->Gaus(0,foc_mirror_nonuni/57.3);
			// 			printf("%12.04e\n",offang);
			// 			rotate_2d(threeSeg3Nz,threeSeg3Ny,cos(offang),sin(offang));
		}
		plane_reflect(\
				threeSeg3Nx,\
				threeSeg3Ny,\
				threeSeg3Nz,\
				threeSeg3D,\
				x,\
				y,\
				z,\
				dx,\
				dy,\
				dz,\
				rval,\
				offang);
		return rval*liquidIndex;
	}
	//reflects off of nothing :(
	z = 1337;
	return -1;
}
float DircOpticalSim::sgn(float val) {
	//stolen from http://stackoverflow.com/questions/1903954/is-there-a-standard-sign-function-signum-sgn-in-c-c
	return (0 < val) - (val < 0);
}
float DircOpticalSim::cylindrical_reflect(\
		float &x,\
		float &y,\
		float &z,\
		float &dx,\
		float &dy,\
		float &dz) {
	//intersects and reflects the focusing cylinder
	//Pretending the cylinder has no Nx for intersection purposes
	//Should be a valid approximation, and saves a ton of speed
	//could be implemented on the threeseg mirror (at least the branching part

	float rval = 0;

	float dydz_norm = sqrt(dz*dz+dy*dy);
	float dy_norm = dy/dydz_norm;
	float dz_norm = dz/dydz_norm;

	float localNx = 0;
	float localNy = 0;
	float localNz = 0;

	//there's gotta be a faster way to do all of this
	//different than plane intercept D.  Took this algorithm from internet (wolfram Mathworld)
	//Parametric intercept sucks
	float D = (z - focMirrorZ)*dy_norm - (y - focMirrorY)*dz_norm;
	float detD = sqrt(foc_r*foc_r - D*D);

	//dy > 0 always (or should be.  Remove sgn and fabs function calls after confirming)
	float zrel = D*dy_norm + sgn(dy_norm)*dz_norm*detD;
	float yrel = -D*dz_norm + fabs(dy_norm)*detD;

	float newy = yrel + focMirrorY;
	float newz = zrel + focMirrorZ;

	float ydiff = newy - y;
	float zdiff = newz - z;

	rval += sqrt((ydiff*ydiff+zdiff*zdiff))/dydz_norm;

	x += dx*rval;
	y = newy;
	z = newz;

	//combine later to save speed
	localNy = yrel;
	localNz = zrel;

	if (nonUniformFocMirror == true) {
		//obviously not in the real run
		float offang = rand_gen->Gaus(0,foc_mirror_nonuni/57.3);
		rotate_2d(localNz,localNy,cos(offang),sin(offang));
	}


	float norm_loc = sqrt(localNy*localNy + localNz*localNz);//there's gotta be a better way to normalize this

	localNy /= norm_loc;
	localNz /= norm_loc;

	//remove for speed in ideal case
	// 	rotate_2d(localNx,localNz,cos(foc_yrot/57.3),sin(foc_yrot/57.3));

	float n_dot_v = -(dx*localNx + dy*localNy + dz*localNz);

	// 	printf("dx: %8.04f dy: %8.04f dz: %8.04f\n",dx,dy,dz);

	dx += 2*n_dot_v*localNx;
	dy += 2*n_dot_v*localNy;
	dz += 2*n_dot_v*localNz;

	// 	printf("dx: %8.04f dy: %8.04f dz: %8.04f\n",dx,dy,dz);

	return rval*liquidIndex;
}

float DircOpticalSim::warp_sens_plane(\
		dirc_point &fill_val,\
		float &x,\
		float &y,\
		float &z,\
		float &dx,\
		float &dy,\
		float &dz) {
	//don't strictly need to modify the z, could be sped up
	//First check to see if it goes through the large planar mirror - it probably doesn't

	float tmpz = get_z_intercept(\
			largePlanarMirrorNx,\
			largePlanarMirrorNy,\
			largePlanarMirrorNz,\
			largePlanarMirrorD,\
			x,\
			y,\
			z,\
			dx,\
			dy,\
			dz);

	if (tmpz < largePlanarMirrorMinZ || tmpz > largePlanarMirrorMaxZ)
	{
//		printf("planarZ: %12.04f\n", tmpz);
		z=1337;
		return 100000;
	}
	float rval =\
		     get_intercept_plane(\
				     sensPlaneNx,\
				     sensPlaneNy,\
				     sensPlaneNz,\
				     sensPlaneD,\
				     x,\
				     y,\
				     z,\
				     dx,\
				     dy,\
				     dz);
	if (z < pmtPlaneMinZ || z > pmtPlaneMaxZ)
	{
//		printf("planarZ: %12.04f\n", tmpz);
		z=1337;
		return 100000;
	}

	fill_val.x = x;
	//fill_val.y = (y-sensPlaneY)*sensPlaneYdistConversion;
	fill_val.y = (-z-559)*sensPlaneYdistConversion + 220;

	return rval*liquidIndex;

}
