#include <vector>

#include "../include/dirc_base_sim.h"
#include "../include/dirc_point.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <stdlib.h>
#include <math.h>
#include <TRandom3.h>

DircBaseSim::DircBaseSim(
		int rand_seed /*=4357*/,\
		float ibarLength/*=4900*/, \
		float ibarWidth/*=35*/, \
		float ibarDepth/*=17.25*/,
		float iupperWedgeTop/*=178.6*/) {

	barLength=ibarLength;
	barWidth=ibarWidth;
	barDepth=ibarDepth;
	upperWedgeTop = iupperWedgeTop;
	wedgeWidthOff = 1.75;
	wedgeHeight = 91;
	wedgeDepthHigh = 79;
	upperWedgeAngleStore = false;
	upperWedgeGap = 20;
	
	wedgeWidth=barWidth - wedgeWidthOff;
	wedgeDepthHigh = wedgeDepthOff+barDepth+wedgeHeight*sin(wedgeCloseAngle/57.296);
	upperWedgeBottom = wedgeHeight + windowThickness + upperWedgeGap;
	upperWedgeHeight = upperWedgeTop - upperWedgeBottom + upperWedgeGap;
	upperWedgeDepthHigh = wedgeDepthHigh + (upperWedgeTop-wedgeHeight)*sin(wedgeCloseAngle/57.296);
	lowerWedgeExtensionZ = -wedgeDepthHigh - \
	    tan(wedgeCloseAngle/57.296)*(upperWedgeBottom - wedgeHeight);

	//Variables used for plane intersection
	wedgeClosePlaneNx = 0; //Shouldn't be needed
	wedgeClosePlaneNy = sin(wedgeCloseAngle/57.296);
	wedgeClosePlaneNz = cos(wedgeCloseAngle/57.296);

	upperWedgeClosePlaneNx = 0; //Shouldn't be needed
	upperWedgeClosePlaneNy = sin(wedgeCloseAngle/57.296);
	upperWedgeClosePlaneNz = cos(wedgeCloseAngle/57.296);

	upperWedgeNonUniform = false;
	upperWedgeNonUniformSpread = 0;

	wedgeClosePlaneD = barLength/2*wedgeClosePlaneNy - wedgeClosePlaneNz * (barDepth+wedgeDepthOff);

	//printf("Wedge Nxyz: %12.04f %12.04f %12.04f\n",wedgeClosePlaneNx,wedgeClosePlaneNy,wedgeClosePlaneNz);

	//upperWedgeClosePlaneD = (barLength/2 + upperWedgeBottom)*upperWedgeClosePlaneNy + upperWedgeClosePlaneNz*lowerWedgeExtensionZ;
	//printf("uwb: %12.04f lwez: %12.04f\n",upperWedgeBottom,lowerWedgeExtensionZ);
	upperWedgeClosePlaneD = wedgeClosePlaneD; //should be in the same plane/;

	upperWedgeFarPlaneNx = 0; //Shouldn't be needed
	upperWedgeFarPlaneNy = 0;
	upperWedgeFarPlaneNz = -1;

	upperWedgeFarPlaneD = (barLength/2 + upperWedgeBottom)*upperWedgeFarPlaneNy;

	//Take average
	quartzIndex = 1.47;
	liquidIndex = 1.47;
	quartzLiquidY = upperWedgeBottom - upperWedgeGap;

	use_liquid_n = false;
	use_quartz_n_for_liquid = false;


	//Negative to make reflections easier
	wedgeFarPlaneNx = 0; //Shouldn't be needed
	wedgeFarPlaneNy = -sin(wedgeFarAngle/57.296);
	wedgeFarPlaneNz = -cos(wedgeFarAngle/57.296);


	wedgeFarPlaneD = barLength/2*wedgeFarPlaneNy;

	upperWedgeFarZ = 0;//Change based on geometry

	rand_gen = std::make_unique<TRandom3>(rand_seed);

	box_angle_off_cval = 1;
	box_angle_off_sval = 0;

	midLineMode = false;
	midLineWedgeWallFlip = 1;

	store_bounces = false;
	x_bounces.clear();
	z_bounces.clear();
	x_direct_bounces.clear();
	z_direct_bounces.clear();
	x_indirect_bounces.clear();
	z_indirect_bounces.clear();

	build_system();
}
void DircBaseSim::set_store_bounces(bool isb)
{
	store_bounces = isb;
}
void DircBaseSim::rotate_2d(float &x, float &y, float cval, float sval) {
	//Standard rotatation allows precomputation of matrix elements
	//Only store one variable for speed
	//Sin should be actual sin (not negative sin)
	float tx = x;
	x = cval*tx - sval*y;
	y = sval*tx + cval*y;
}

void DircBaseSim::set_bar_box_angle(float ang) {
	//expect degrees
	//sets angle between readout box and bars - rotates angle coming out of bars
	box_angle_off_cval = cos(ang/57.3);
	box_angle_off_sval = sin(ang/57.3);
}
void DircBaseSim::set_bar_box_offsets(float x, float y, float z) {
	//expect degrees
	//sets angle between readout box and bars - rotates angle coming out of bars
	bar_box_xoff = x;
	bar_box_yoff = y;
	bar_box_zoff = z;
}
std::vector<float> DircBaseSim::get_dist_traveled() {
	return dist_traveled;
}
void DircBaseSim::set_liquid_index(float li) {
	//Sets the index of the upper quartz wedge
	liquidIndex = li;
	//printf("Liquid index set: %12.04f\n",liquidIndex);

}
void DircBaseSim::set_upper_wedge_angle_diff(float rads, float rads_y) {
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
void DircBaseSim::set_wedge_mirror_rand(float ispread) {
	if (ispread > 0) {
		upperWedgeNonUniform = true;
		upperWedgeNonUniformSpread = ispread;
	} else {
		upperWedgeNonUniform = false;
		upperWedgeNonUniformSpread = 0;
	}
	spread_wedge_mirror();
}
void DircBaseSim::spread_wedge_mirror() {

	if (upperWedgeNonUniform == true) {
		float spread_ang = wedgeCloseAngle;
		spread_ang += rand_gen->Gaus(0,upperWedgeNonUniformSpread);
		upperWedgeClosePlaneNx = 0; //Shouldn't be needed
		upperWedgeClosePlaneNy = sin(spread_ang/57.296);
		upperWedgeClosePlaneNz = cos(spread_ang/57.296);
		upperWedgeClosePlaneD = (barLength/2 + upperWedgeBottom)*upperWedgeClosePlaneNy + upperWedgeClosePlaneNz*lowerWedgeExtensionZ;
	}
}
void DircBaseSim::build_system() {
	spread_wedge_mirror();
}
float DircBaseSim::get_quartz_n(float lambda) {
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
void DircBaseSim::set_use_quartz_n_for_liquid(bool iu)
{
	use_quartz_n_for_liquid = iu;
}
float DircBaseSim::get_liquid_n(float lambda) 
{
	if (use_liquid_n == true)
	{
		return liquidIndex;
	}
	if (use_quartz_n_for_liquid == true)
	{
		return get_quartz_n(lambda);
	}
	else
	{
		//2 line approximation
		float ln1 = 1.3786 - 0.1*lambda/1000;
		float ln2 = 1.3529 - 0.0353*lambda/1000;
		float rval = std::max(ln1,ln2);
		return rval;
	}
}
bool DircBaseSim::quartz_transmission_mc(float R, float lambda) {
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
float DircBaseSim::get_beta(float E, float m) {
	float gam = E/m;
	if (gam > 5) {
		//Approximate and possibly save time
		return 1-1/(2*gam*gam);
	} else {
		return sqrt(1-1/(gam*gam));
	}

}
float DircBaseSim::get_bar_offset(int bar)
{
	return fabs(bar)/bar*((150-0.5*(barWidth))+bar*(barWidth+.015));
}
int DircBaseSim::get_bar_from_x(float x)
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

void DircBaseSim::set_upper_wedge_angle_store(bool istore)
{
	upperWedgeAngleStore = istore;
}

std::vector<float> DircBaseSim::get_upper_wedge_incident()
{
	return upper_wedge_incident;
}

void DircBaseSim::test_from_wedge_top(\
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
void DircBaseSim::bar_box_interface(
		float &x,
		float &y,
		float &z,
		float &dx,
		float &dy,
		float &dz)
{
	x += bar_box_xoff;
	y += bar_box_yoff;
	z += bar_box_zoff;
	rotate_2d(dy,dz,box_angle_off_cval,box_angle_off_sval);
}
void DircBaseSim::sim_lut_points(\
		std::vector<dirc_point> &ovals,\
		std::vector<float> &phis,\
		std::vector<float> &thetas,\
		int n_photons, \
		float particle_bar /*= 0*/){

	ovals.clear();
	phis.clear();
	thetas.clear();


	float x,y,z,dx,dy,dz;

	float randPhi;
	float randTheta;	

	float mm_index = 0;

//	float c_mm_ns = 300;

	float saveGeneralQuartzIndex = quartzIndex;
	float saveGeneralLiquidIndex = liquidIndex;
	//Using approximate peak wavelength
	quartzIndex = get_quartz_n(380);//Painful way of doing this - saving and correcting is inelegant
	liquidIndex = get_liquid_n(380);//Painful way of doing this - saving and correcting is inelegant

	//printf("LUT q and L: %12.04f %12.04f\n",quartzIndex,liquidIndex);

	for (int i = 0; i < n_photons; i++) {
		randPhi = rand_gen->Uniform(0,2*3.14159265);
		randTheta = acos(2*rand_gen->Uniform(.5,1) - 1);
		//Useful for directly testing box, not fo lut
		//		randTheta = 45/57.3;
		//		randPhi = -3.1415926535*(180+45)/180.;
		//This timing won't be correct
		mm_index = 0;

		x = 0;
		y = barLength/2;
		z = -barDepth/2;


		dx = sin(randTheta)*sin(randPhi);
		dy = cos(randTheta);
		dz = sin(randTheta)*cos(randPhi); 

		//printf("lutstart: %12.04f %12.04f %12.04f\n",dx,dy,dz);
		/*
		   dx = 0;
		   dy = 1;
		   dz = 0;
		 */

		mm_index += warp_wedge(\
				x,\
				y,\
				z,\
				dx,\
				dy,\
				dz);
		//account (quickly) for the bar box having a different angle than the readout
		//rotate_2d(dy,dz,box_angle_off_cval,box_angle_off_sval);
		bar_box_interface(x,y,z,dx,dy,dz);
		if (z > 0) {
			continue;
		}
		dirc_point out_val;



		warp_readout_box(out_val,particle_bar,mm_index,x,y,z,dx,dy,dz);

		if (out_val.t < 0)
		{
			continue;
		}

		//Timing is hard...
		//out_val.t = mm_index/c_mm_ns;
		out_val.updown = 0;
		ovals.push_back(out_val);
		phis.push_back(randPhi);
		thetas.push_back(randTheta);
	}
	quartzIndex = saveGeneralQuartzIndex;
	liquidIndex = saveGeneralLiquidIndex;
}
void DircBaseSim::fill_rand_phi(
		std::back_insert_iterator<std::vector<dirc_point>> &ovals,
		int n_photons,
		float ckov_theta /*= 47*/,
		float particle_bar /*= 0*/,
		float particle_x /*= 0*/,
		float particle_y /*= 0*/,
		float particle_t /*= 0*/,
		float particle_theta /*= 0*/,
		float particle_phi /*= 0*/,
		float phi_theta_unc, /*= .0015*57.3*/
		float ckov_theta_unc /* = .0055*57.3*/,
		float beta/* = -1*/) {

	float emitAngle = ckov_theta;
	float particleTheta = particle_theta + rand_gen->Gaus(0,phi_theta_unc);
	float particlePhi = particle_phi + rand_gen->Gaus(0,phi_theta_unc);
	int numPhots = n_photons/cos(particle_theta/57.3);

	float randPhi;

	float temit, rand_add;
	float wavelength = 400;

	float x,y,z,dx,dy,dz;

	float cos_ptheta = cos(particleTheta/57.3);
	float sin_ptheta = sin(particleTheta/57.3);
	float cos_pphi = cos(particlePhi/57.3);
	float sin_pphi = sin(particlePhi/57.3);

	float mm_index = 0;

	int tmp_updown = 0;

	float saveGeneralQuartzIndex = quartzIndex;
	float saveGeneralLiquidIndex = liquidIndex;

	std::vector<dirc_base_sim_tracking_step> track_steps;
	float dist_traveled = -1;
	//hardcode step length to 1mm for now
	//note: These are general tracking vectors - you can implement your own MC with them.  Each refers to position and direction at the start of the step

	dirc_base_sim_tracking_step step1;
	step1.x = particle_x;	
	step1.y = particle_y;	
	step1.z = -barDepth;

	step1.sin_theta = sin_ptheta;
	step1.cos_theta = cos_ptheta;
	step1.sin_phi = sin_pphi;
	step1.cos_phi = cos_pphi;

	dist_traveled = barDepth/cos_ptheta;
	//take_1_step
	const float step_length = barDepth/cos_ptheta;

	dirc_base_sim_tracking_step step2;
	step2.x = particle_x + step_length*sin_ptheta*cos_pphi;	
	step2.y = particle_y + step_length*sin_ptheta*sin_pphi;	
	step2.z = -barDepth + step_length*cos_ptheta;

	step2.sin_theta = sin_ptheta;
	step2.cos_theta = cos_ptheta;
	step2.sin_phi = sin_pphi;
	step2.cos_phi = cos_pphi;

	track_steps.push_back(step1);
	track_steps.push_back(step2);	

	float track_loc = -1;

	//define numbers to linearly interp
	int low_ind = 0;
	float above_ind = 0;

	for (int i = 0; i < numPhots; i++) {
		randPhi = rand_gen->Uniform(0,2*3.14159265);
		track_loc = rand_gen->Uniform(0,dist_traveled);
		if (beta < 0) {
			rand_add = rand_gen->Gaus(0,ckov_theta_unc);
			temit = emitAngle + rand_add;
		} else {
			temit = get_cerenkov_angle_rand(beta,ckov_theta_unc,wavelength);
			quartzIndex = get_quartz_n(wavelength);//Painful way of doing this - saving and correcting is inelegant
			liquidIndex = get_liquid_n(wavelength);//Painful way of doing this - saving and correcting is inelegant
		}
//		mm_index = (sourceOff - barDepth)*quartzIndex/cos(particleTheta/57.3);

		
		//this will actually floor - track_loc > 0
		low_ind = (int) (track_loc/step_length);
		above_ind = (track_loc - low_ind*step_length);

		cos_ptheta = track_steps[low_ind].cos_theta;
		sin_ptheta = track_steps[low_ind].sin_theta;
		cos_pphi = track_steps[low_ind].cos_phi;
		sin_pphi = track_steps[low_ind].sin_phi;

		x = track_steps[low_ind].x + above_ind*sin_ptheta*cos_pphi;
		y = track_steps[low_ind].y + above_ind*sin_ptheta*sin_pphi;
		z = track_steps[low_ind].z + above_ind*cos_ptheta;

		mm_index = (barDepth+z)*quartzIndex/cos(particleTheta/57.3);
//		printf("particle xyz: %12.04f %12.04f %12.04f\n",x,y,z);

		dx = sin(temit/57.3)*cos(randPhi);
		dy = sin(temit/57.3)*sin(randPhi);
		dz = cos(temit/57.3);

		rotate_2d(dz,dy,cos_ptheta,sin_ptheta);
		rotate_2d(dy,dx,cos_pphi,sin_pphi);


		//photon is now defined as up or down
		if (dy > 0)
		{
			tmp_updown = 1;
		}
		else
		{//Assume that zero doesn't count
			tmp_updown = -1;
		}

		mm_index += warp_ray(
				x,
				y,
				z,
				dx,
				dy,
				dz,
				sqrt(1-1/(1.47*1.47)));

		if (z > 0)
		{
			continue;
		}
		spread_wedge_mirror();

		mm_index += warp_wedge(
				x,
				y,
				z,
				dx,
				dy,
				dz);
		//account (quickly) for the bar box having a different angle than the readout
		//rotate_2d(dy,dz,box_angle_off_cval,box_angle_off_sval);
		bar_box_interface(x,y,z,dx,dy,dz);

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
		// ovals is a push_back_iterator
		ovals = out_val;
	}
	quartzIndex = saveGeneralQuartzIndex;
	liquidIndex = saveGeneralLiquidIndex;
}
std::vector<std::pair<float,float> > DircBaseSim::get_refraction_rand_phi(\
		std::vector<float> &before_interface,\
		std::vector<float> &after_interface,\
		std::vector<float> &pmt_incidence,\
		int n_photons, \
		float ckov_theta /*= 47*/, \
		float particle_x /*= 0*/, \
		float particle_y /*= 0*/, \
		float particle_theta /*= 0*/, \
		float particle_phi /*= 0*/,\
		float phi_theta_unc /*= .0015*57.3*/,\
		float ckov_theta_unc /* = .0055*57.3*/,\
		float beta/* = -1*/) {
	//returns theta versus cerenkov phi_theta_unc

	std::vector<std::pair<float,float> > rval;

	before_interface.clear();
	after_interface.clear();
	refraction_before.clear();
	refraction_after.clear();
	pmt_incidence.clear();

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
		mm_index = (sourceOff - barDepth)*quartzIndex;

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
		//rotate_2d(dy,dz,box_angle_off_cval,box_angle_off_sval);
		bar_box_interface(x,y,z,dx,dy,dz);

		if (z > 0) {
			continue;
		}
		dirc_point out_val;

		warp_readout_box(\
				out_val,\
				0,\
				mm_index,\
				x,\
				y,\
				z,\
				dx,\
				dy,\
				dz);
		if (out_val.t < 0)
		{
			continue;
		}

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

	refraction_before.clear();
	refraction_after.clear();

	return rval;
}

void DircBaseSim::fill_reg_phi(
		std::back_insert_iterator<std::vector<dirc_point>> &fill_points,
		int n_photons_phi,
		int n_photons_z,
		float ckov_theta /*= 47*/,
		float particle_bar /*= 0*/,
		float particle_x /*= 0*/,
		float particle_y /*= 0*/,
		float particle_t /*= 0*/,
		float particle_theta /*= 0*/,
		float particle_phi /*= 0*/,
		float phi_theta_unc /*= 0*/,
		float ckov_theta_unc /* = 0*/,
		float beta /* = -1*/)
{
	float sDepth = .95*barDepth;
	float emitAngle = ckov_theta;

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
	float saveGeneralQuartzIndex = quartzIndex;
	float saveGeneralLiquidIndex = liquidIndex;

	for (int i = 0; i < n_photons_z; i++) {
		sourceOff = (i+.5)*sDepth/(n_photons_z);
		for (int j = 0; j < adj_n_photons_phi; j++) {
			regPhi = j*2*3.14159265357/(adj_n_photons_phi);
			if (beta < 0) {
				rand_add = rand_gen->Gaus(0,ckov_theta_unc);
				temit = emitAngle + rand_add;
			} else {
				temit = get_cerenkov_angle_rand(beta,ckov_theta_unc,wavelength);
				//Painful way of doing this - saving and correcting is inelegant
				quartzIndex = get_quartz_n(wavelength);
				//Painful way of doing this - saving and correcting is inelegant
				liquidIndex = get_liquid_n(wavelength);
			}

			mm_index = (sourceOff - barDepth)*quartzIndex/cos_ptheta;

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
			mm_index += warp_ray(
					x,
					y,
					z,
					dx,
					dy,
					dz,
					sqrt(1-1/(1.47*1.47)));
			if (z > 0)
			{
				continue;
			}
			mm_index += warp_wedge(
					x,
					y,
					z,
					dx,
					dy,
					dz);
			//account (quickly) for the bar box having a different angle than the readout
			bar_box_interface(x,y,z,dx,dy,dz);

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
			// Not the niciest syntax
			// but this insterts into the output vector
			fill_points = out_val;
		}
	}
	quartzIndex = saveGeneralQuartzIndex;
	liquidIndex = saveGeneralLiquidIndex;
}

bool DircBaseSim::track_single_photon(\
		dirc_point &out_val,\
		float emit_theta,\
		float emit_phi,\
		float particle_theta,\
		float particle_phi,\
		float particle_x,\
		float particle_y,\
		float particle_z,\
		float particle_t,\
		int particle_bar)
{
	float x,y,z,dx,dy,dz;
	float mm_index;

	float temit = emit_theta;
	float regPhi = emit_phi;

	float cos_emit, sin_emit;
	float cos_regphi, sin_regphi;
	int tmp_updown = 0;

	float cos_ptheta = cos(particle_theta/57.3);
	float sin_ptheta = sin(particle_theta/57.3);
	float cos_pphi = cos(particle_phi/57.3);
	float sin_pphi = sin(particle_phi/57.3);

	float sourceOff = particle_z + barDepth;

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
		return false;
	}

	spread_wedge_mirror();


	if (dx < 0)
	{
		lastWallX = 1;
	}
	else
	{
		lastWallX = -1;
	}
	//int beforeWedgeLastWall = lastWallX;

	mm_index += warp_wedge(\
			x,\
			y,\
			z,\
			dx,\
			dy,\
			dz);

	if (dx < 0)	
	{
		lastWallX = 1;
	}
	else
	{
		lastWallX = -1;
	}
	out_val.last_wall_x = lastWallX;
	//	printf("lastwall agrees: %d\n",beforeWedgeLastWall==lastWallX);
	//account (quickly) for the bar box having a different angle than the readout
	//rotate_2d(dy,dz,box_angle_off_cval,box_angle_off_sval);
	bar_box_interface(x,y,z,dx,dy,dz);

	if (z > 0)
	{
		return false;
	}

	out_val.wedge_before_interface = -1;
	warp_readout_box(out_val,particle_bar,mm_index,x,y,z,dx,dy,dz);

	out_val.wedge_before_interface = wedgeBeforeInterface;

	if (out_val.t < 0)
	{
		return false;
	}

	out_val.t += particle_t;
	out_val.updown = tmp_updown;

	return true;
}
bool DircBaseSim::track_line_photon(\
		dirc_point &out_val,\
		float emit_theta,\
		float emit_phi,\
		float particle_theta,\
		float particle_phi,\
		float particle_x,\
		float particle_y,\
		float particle_z,\
		float particle_t,\
		int particle_bar,\
		float z_at_top /*=1*/)
{
	//I'm being bad and breaking encapsulation here, but it's for the greater good
	//Think of the children
	float saveBarWidth = barWidth;
	float saveWedgeWidthOff = wedgeWidthOff;

	midLineMode = true;
	midLineMode = false;
	float width_fraction = 1;
	barWidth *= width_fraction;
	wedgeWidthOff *= width_fraction;
	wedgeWidthOff = 0;
	build_system();

	float x,y,z,dx,dy,dz;
	float mm_index;

	float temit = emit_theta;
	float regPhi = emit_phi;

	float cos_emit, sin_emit;
	float cos_regphi, sin_regphi;
	int tmp_updown = 0;

	float cos_ptheta = cos(particle_theta/57.3);
	float sin_ptheta = sin(particle_theta/57.3);
	float cos_pphi = cos(particle_phi/57.3);
	float sin_pphi = sin(particle_phi/57.3);

	float sourceOff = particle_z + barDepth;

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
		barWidth = saveBarWidth;
		wedgeWidthOff = saveWedgeWidthOff;
		build_system();
		midLineMode = false;
		return false;
	}

	spread_wedge_mirror();


	if (dx < 0)
	{
		lastWallX = 1;
	}
	else
	{
		lastWallX = -1;
	}
//	int beforeWedgeLastWall = lastWallX;

/*
	float target_z = -17.25/2;

	if (z_at_top < 0)
	{
		target_z = z_at_top;
	}

	//z = target_z;
	if (dz > 0)
	{
		//y += dy/dz*(-z);
		//	dz = -dz;
	}
*/

	mm_index += warp_wedge(\
			x,\
			y,\
			z,\
			dx,\
			dy,\
			dz);

	if (dx < 0)	
	{
		lastWallX = 1;
	}
	else
	{
		lastWallX = -1;
	}
	out_val.last_wall_x = lastWallX;
	//	printf("lastwall agrees: %d\n",beforeWedgeLastWall==lastWallX);
	//account (quickly) for the bar box having a different angle than the readout
	//rotate_2d(dy,dz,box_angle_off_cval,box_angle_off_sval);
	bar_box_interface(x,y,z,dx,dy,dz);

	if (z > 0)
	{
		barWidth = saveBarWidth;
		wedgeWidthOff = saveWedgeWidthOff;
		build_system();
		midLineMode = false;
		return false;
	}
	if (lastWallX == 1)
	{
		//	x += saveBarWidth/2 - 3;
	}
	else if (lastWallX == -1)
	{
		//x += saveBarWidth/2 - saveWedgeWidthOff;
		//	x += saveBarWidth/2 - 3;
	}	
	out_val.wedge_before_interface = -1;
	warp_readout_box(out_val,particle_bar,mm_index,x,y,z,dx,dy,dz);

	out_val.wedge_before_interface = wedgeBeforeInterface;

	if (out_val.t < 0)
	{
		barWidth = saveBarWidth;
		wedgeWidthOff = saveWedgeWidthOff;
		build_system();
		midLineMode = false;
		return false;
	}

	out_val.t += particle_t;
	out_val.updown = tmp_updown;

	barWidth = saveBarWidth;
	wedgeWidthOff = saveWedgeWidthOff;
	build_system();
	midLineMode = false;
	return true;
}
bool DircBaseSim::track_single_photon_beta(\
		dirc_point &out_val,\
		float particle_beta,\
		float emit_phi,\
		float particle_theta,\
		float particle_phi,\
		float particle_x,\
		float particle_y,\
		float particle_z,\
		float particle_t,\
		int particle_bar)
{
	float x,y,z,dx,dy,dz;
	float mm_index;

	float wavelength = 400;

	float temit = get_cerenkov_angle_rand(particle_beta,0,wavelength);
	float regPhi = emit_phi;

	float cos_emit, sin_emit;
	float cos_regphi, sin_regphi;
	int tmp_updown = 0;

	float cos_ptheta = cos(particle_theta/57.3);
	float sin_ptheta = sin(particle_theta/57.3);
	float cos_pphi = cos(particle_phi/57.3);
	float sin_pphi = sin(particle_phi/57.3);

	float sourceOff = particle_z + barDepth;

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
		return false;
	}

	spread_wedge_mirror();


	if (dx < 0)
	{
		lastWallX = 1;
	}
	else
	{
		lastWallX = -1;
	}
//	int beforeWedgeLastWall = lastWallX;

	mm_index += warp_wedge(\
			x,\
			y,\
			z,\
			dx,\
			dy,\
			dz);

	if (dx < 0)	
	{
		lastWallX = 1;
	}
	else
	{
		lastWallX = -1;
	}
	out_val.last_wall_x = lastWallX;
	//	printf("lastwall agrees: %d\n",beforeWedgeLastWall==lastWallX);
	//account (quickly) for the bar box having a different angle than the readout
	//rotate_2d(dy,dz,box_angle_off_cval,box_angle_off_sval);
	bar_box_interface(x,y,z,dx,dy,dz);

	if (z > 0)
	{
		return false;
	}

	out_val.wedge_before_interface = -1;
	warp_readout_box(out_val,particle_bar,mm_index,x,y,z,dx,dy,dz);

	out_val.wedge_before_interface = wedgeBeforeInterface;

	if (out_val.t < 0)
	{
		return false;
	}

	out_val.t += particle_t;
	out_val.updown = tmp_updown;

	return true;
}
bool DircBaseSim::track_all_line_photons(\
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
		float z_at_top /*=1*/)
{
	//I'm being bad and breaking encapsulation here, but it's for the greater good
	//Think of the children
	float saveBarWidth = barWidth;
	float saveWedgeWidthOff = wedgeWidthOff;

	left_vals.clear();
	right_vals.clear();

	midLineMode = true;
	float width_fraction = .001;
	barWidth *= width_fraction;
	wedgeWidthOff *= width_fraction;
	wedgeWidthOff = 0;
	build_system();

	std::vector<float> mustache_phi;

	float x,y,z,dx,dy,dz;
	//righthand side versions - right hand not needed, just reuse
	//float x_r,y_r,z_r,dx_r,dy_r,dz_r;
	float mm_index;
//	float mm_index_r;

	float temit = emit_theta;
	float regPhi;

	float cos_emit, sin_emit;
	float cos_regphi, sin_regphi;
	int tmp_updown = 0;

	float cos_ptheta = cos(particle_theta/57.3);
	float sin_ptheta = sin(particle_theta/57.3);
	float cos_pphi = cos(particle_phi/57.3);
	float sin_pphi = sin(particle_phi/57.3);

	float sourceOff = particle_z + barDepth;

	float phi_inc = (2*3.14159265359)/points_per_side;
	regPhi = 0;
	float saveGeneralQuartzIndex = quartzIndex;
	float saveGeneralLiquidIndex = liquidIndex;
	//Using approximate peak wavelength
	quartzIndex = get_quartz_n(380);//Painful way of doing this - saving and correcting is inelegant
	liquidIndex = get_liquid_n(380);//Painful way of doing this - saving and correcting is inelegant
	for (int i = 0; i < points_per_side; i++)
	{
		x = 0;
		y = 0;
		z = sourceOff;
		mm_index = (sourceOff - barDepth)*1.47;
		regPhi += phi_inc;

		//save some time ~31ms per 400k
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


//		int beforeWedgeLastWall = lastWallX;


		float target_z = -17.25/2;

		if (z_at_top < 0)
		{
			target_z = z_at_top;
		}

		z = target_z;
		if (dz > 0)
		{
			dz = -dz;
		}

		dirc_point out_val;

		//branch here

		if (upperWedgeBottom/(upperWedgeBottom*tan(wedgeCloseAngle/57.3)+1.5*barDepth+wedgeDepthOff) < -dy/dz)
		{
			mustache_phi.push_back(regPhi);
			//	printf("%12.04f\n",regPhi);
		}

		mm_index += warp_wedge(\
				x,\
				y,\
				z,\
				dx,\
				dy,\
				dz);

//		mm_index_r = mm_index;

		int dx_lr = sgn(dx);

		out_val.last_wall_x = -1;
		//rotate_2d(dy,dz,box_angle_off_cval,box_angle_off_sval);
		bar_box_interface(x,y,z,dx,dy,dz);

		x += saveBarWidth/2-3;

		//And this set of ifs is why we try and do one photon at a time normally
		if (z < 0)
		{
			out_val.wedge_before_interface = -1;
			warp_readout_box(out_val,particle_bar,mm_index,x,y,z,dx,dy,dz);
			out_val.wedge_before_interface = wedgeBeforeInterface;
			if (out_val.t > 0)
			{
				out_val.t += particle_t;
				out_val.updown = tmp_updown;
				out_val.init_phi = regPhi;
				if (dx_lr < 0)
				{
					left_vals.push_back(out_val);
				}
				else
				{
					right_vals.push_back(out_val);
				}
			}
		}
	}
	for (unsigned int i = 0; i < mustache_phi.size(); i++)
	{
		x = 0;
		y = 0;
		z = sourceOff;
		mm_index = (sourceOff - barDepth)*1.47;
		regPhi = mustache_phi[i];


		//save some time ~31ms per 400k
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


//		int beforeWedgeLastWall = lastWallX;

		//now in the mustache => at wall
		//And going up
		float target_z = 0;

		if (z_at_top < 0)
		{
			target_z = z_at_top;
		}

		z = target_z;
		if (dz > 0)
		{
			dz = -dz;
		}
		//Goes backwards first and reflects off of wall.  Could account for 6mrad here as well, consider later.
		y += -dy/dz*barDepth;
		//perform proper reflection - factor of 2 for the reflection, and an additional 1 for the corrected wedge.
		float enhance_rot_fac = 3;
		rotate_2d(dy,dz,cos(enhance_rot_fac*wedgeFarAngle/57.3),-sin(enhance_rot_fac*wedgeFarAngle/57.3));

		dirc_point out_val;

		//branch here

		mm_index += warp_wedge(\
				x,\
				y,\
				z,\
				dx,\
				dy,\
				dz);

		//mm_index_r = mm_index;

		int dx_lr = sgn(dx);

		out_val.last_wall_x = -1;
		//rotate_2d(dy,dz,box_angle_off_cval,box_angle_off_sval);
		bar_box_interface(x,y,z,dx,dy,dz);

		x += saveBarWidth/2-3;

		//And this set of ifs is why we try and do one photon at a time normally
		if (z < 0)
		{
			out_val.wedge_before_interface = -1;
			warp_readout_box(out_val,particle_bar,mm_index,x,y,z,dx,dy,dz);
			out_val.wedge_before_interface = wedgeBeforeInterface;
			if (out_val.t > 0)
			{
				out_val.t += particle_t;
				out_val.updown = tmp_updown;
				out_val.init_phi = regPhi;
				if (dx_lr < 0)
				{
					left_vals.push_back(out_val);
				}
				else
				{
					right_vals.push_back(out_val);
				}
			}
		}
	}
	barWidth = saveBarWidth;
	wedgeWidthOff = saveWedgeWidthOff;
	quartzIndex = saveGeneralQuartzIndex;
	liquidIndex = saveGeneralLiquidIndex;

	build_system();
	midLineMode = false;
	return true;
}
float DircBaseSim::warp_ray(\
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

	bool direct_ray = true;

	if (dy > 0) {
		//Going up
		dely = barLength*(0.5-grace_room) - y;
		direct_ray = true;
	} else {
		//going down, flip dy
		dely = barLength*(1.5-grace_room) + y;
		dy = -dy;
		direct_ray = false;
	}

	rval += dely*dely;

	y = barLength*(.5-grace_room);

	if (fabs(dz) > cos_critical_angle ||\
			fabs(dx) > cos_critical_angle ||\
			dy*dy < 1e-4) {
		//If it's not totally internally reflected, assume it escapes
		//also assume failure if it isn't going up fast enough (absorbed)
		//Positive z origin signals fail
		//Do not check y - is not totally internally reflected on the bottom - there is an actuall mirror there.
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


	if (store_bounces == true)
	{
		x_bounces.push_back(nbouncesx);
		z_bounces.push_back(nbouncesz);
		if (direct_ray == true)
		{
			x_direct_bounces.push_back(nbouncesx);	
			z_direct_bounces.push_back(nbouncesz);	
		}
		else
		{
			x_indirect_bounces.push_back(nbouncesx);	
			z_indirect_bounces.push_back(nbouncesz);	
		}
	}



	//	float bar_front = -barDepth/2;
	//	float bar_front = -17.25/2;
	float bar_front = 0;
	if (nbouncesz % 2 == 1) {
		dz = -dz;
		z = bar_front-remainderz;
	} else {
		z = remainderz-barDepth;
	}

	return sqrt(rval)*quartzIndex;
}

float DircBaseSim::warp_wedge(
		float &x,
		float &y,
		float &z,
		float &dx,
		float &dy,
		float &dz) {
	// No Critical angle checking - shouldn't need it
	// I am abusing the x angle being the same in the bar and wedge here
	// may have to rotate after to the mirror
	// Starts y at 0

	wedge_bounces = 0;

	//must change after close wedge bounce
	float startz = z;
	float starty = y;
	float mm_index = 0;
	bool passed_interface = false;
	float dt = 0; //dummy variable representing parametric line
	float n_dot_v = 0; //dummy variable to save dot product
	float n_dot_v0 = 0;

	//deal with yz first

	//Can't ever not totally internal reflect	
	//force dz to be negative	
	//float wedgeCloseIncidentAngle = acos(-dx*wedgeClosePlaneNx - dy*wedgeClosePlaneNy + fabs(dz)*wedgeClosePlaneNz);
	//printf("%12.04f\n",57.3*wedgeCloseIncidentAngle);

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


	//have to make sure ty intersects IN the wedge - will correctly return a negative if it's below
	if (ty < wedgeHeight && ty > 0) {
		//reflect it off bottom wedge

		//Again, always true
		mm_index += dt*quartzIndex;
		if (!(x_wedge_coerce_check(x,y,z,dx,dy,dz,dt))) {
			//Goes out the window, fail and return
			//Shouldn't
			z = 1337;
			return -1;
		}
		//printf("ty: %12.04f z-bd/2: %12.04f\n",ty,z+barDepth/2);
		dy += 2*n_dot_v*wedgeClosePlaneNy;
		dz += 2*n_dot_v*wedgeClosePlaneNz;

		//start the bouncing here
		starty = y;
		startz = z;

		wedgeBeforeInterface = 1;

	} 
	else if (ty > upperWedgeBottom && ty < upperWedgeTop) 
	{
		//passed interface before reflection
		//get dt to propagate to interface
		dt = (quartzLiquidY + barLength/2 - y)/dy;
		mm_index += dt*quartzIndex;
		if (!(x_wedge_coerce_check(x,y,z,dx,dy,dz,dt))) 
		{
			//Goes out the window, fail and return
			z = 1337;
			return -1;
		}
		//printf("top_of_lowerwedge_postinter_dxdydz %12.04f %12.04f %12.04f\n",dx,dy,dz);
		//correct ordering below - interface is xz plane, so dx and dz go first
		optical_interface_z(quartzIndex,liquidIndex,dx,dz,dy);
		passed_interface = true;

		wedgeBeforeInterface = 0;

		//NOW intersect with upper wedge

		//Signs are weird here, it works, but I think they are twisted up with the initialization too
		n_dot_v = -(dy*upperWedgeClosePlaneNy + dz*upperWedgeClosePlaneNz + dx*upperWedgeClosePlaneNx);
		n_dot_v0 = -(y*upperWedgeClosePlaneNy + z*upperWedgeClosePlaneNz + x*upperWedgeClosePlaneNx);


		dt = -(wedgeClosePlaneD+n_dot_v0)/n_dot_v;
		//printf("uwcpd: %12.04f nv: %12.04f nv0: %12.04f dt: %12.04f\n",upperWedgeClosePlaneD,n_dot_v,n_dot_v0,dt);

		if (dt*dy + y < upperWedgeTop + barLength/2) 
		{
			mm_index += dt*liquidIndex;
			if (!(x_wedge_coerce_check(x,y,z,dx,dy,dz,dt))) {
				//Goes out the window, fail and return
				z = 1337;
				return -1;
			}

			if (y > wedgeHeight + barLength/2 && y < upperWedgeBottom + barLength/2)
			{
				//We let it go out the gap on the slanted side, but not on the left and right.
				z = 1337;
				return -1;
			}

			if (upperWedgeAngleStore == true)
			{
				upper_wedge_incident.push_back(57.3*acos(n_dot_v));
			}

			dx += 2*n_dot_v*upperWedgeClosePlaneNx;
			dy += 2*n_dot_v*upperWedgeClosePlaneNy;
			dz += 2*n_dot_v*upperWedgeClosePlaneNz;

			//printf("upper wedge ty: %12.04f z-bd/2: %12.04f uwb: %12.04f uwt: %12.04f\n",ty,z+barDepth/2,upperWedgeBottom,upperWedgeTop);
			//yet again, start after bouncing
			starty = y;
			startz = z;

		} 
		else 
		{
			//refracted such that it no longer bounces
			//Still have to check from here later
			starty = y;
			startz = z;
			dt = (upperWedgeTop + barLength/2 - y)/dy;
			mm_index += dt*liquidIndex;
			if (!(x_wedge_coerce_check(x,y,z,dx,dy,dz,dt))) 
			{
				//Goes out the window, fail and return
				z = 1337;
				return -1;
			}
		}

	} else if (ty < upperWedgeTop && ty > 0) {
		//out the window
		//printf("%12.04f\n",ty);
		z = 1337;
		return -1;
	}

	//Have now performed the maximum number of reflections (besides x direction)
	//Sort of - still have to account for z direction as well
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
		//printf("top_of_lowerwedge_preinter_dxdydz %12.04f %12.04f %12.04f\n",dx,dy,dz);
	}

	//Now we just finish
	//        printf("liquidIndex: %12.04f  quartzIndex: %12.04f\n",liquidIndex,quartzIndex);
	//	printf("%12.04f %12.04f %12.04f %12.04f %12.04f %12.04f\n",x,y,z,dx,dy,dz);

	if (wedge_bounces > 2)
	{	
		//		printf("Wedge Before Interface   :  bounces: %d  :  %d\n",wedgeBeforeInterface,wedge_bounces);
	}
	dt = (upperWedgeTop + barLength/2 - y)/dy;
	mm_index += dt*liquidIndex;
	if (!(x_wedge_coerce_check(x,y,z,dx,dy,dz,dt)) || dt > 1000) {
		//Cut out large dts - more than a m (~30 reflections) in the wedge - only seems to show up in very low dy situations from the LUT table.
		//I'm not sure if it can go out the window at this point, but just in case
		//Goes out the window, fail and return
		//y position when hitting the back wall

		z = 1337;
		return -1;
	}


	//and we're done.  Everything should be correct now with the ray at the top of the wedge
	//return the distance traveled times the index it traveled in

	//printf("DTDTDT: %12.04f\n",dt);

	if (z > 0)
	{
		if (dz < 0)
		{
			printf("WTH Over: Somethings wrong with Zs in the wedge %d\n",wedgeBeforeInterface);
			printf("wthxyz: %12.04f %12.04f %12.04f %12.04f %12.04f %12.04f dt: %12.04f\n",x,y,z,dx,dy,dz,dt);
		}
		//implement reflection off the vack of the walls
		float wally = starty + (y - starty)*(-startz/(z-startz));
		//account for window and gap
		if (wally > wedgeHeight + barLength/2)
		{
			//bounced off upper part of mirror
			if (wally < upperWedgeBottom + barLength/2)
			{
				//printf("%12.04f\n",wally);
				z = 1337;
				return -1;
			}
			z = -z;
			dz = -dz;
		}
		else
		{
			//float sz = z;
			float sdz = dz;
			float sdy = dy;
			//bounced off lower part - just an appoximation
			dz = -dz;
			//printf("%12.04f %12.04f\n",dz,dy);
			rotate_2d(dz,dy,cos(2*wedgeFarAngle/57.3),sin(2*wedgeFarAngle/57.3));
			//printf("  %12.04f %12.04f\n",dz,dy);
			//z = -z;
			//account for change in angle


			z = z*dz/sdz;
			y = y - (y-wally)*(sdy-dy)/sdy;
		}
	}
	return mm_index;
}

//Possibly inline these or something for speed, but right now, leave them for sanity
bool DircBaseSim::optical_interface_z(\
		float n1,\
		float n2,\
		float &dx,\
		float &dy,\
		float &dz) {
	// n1 is starting index of refraction, n2 is ending
	// Assume that it's going through the z plane
	// (this will change the ordering when called in our current geometry)

	//dz and dy are flipped, so this is really acos
	float n12rat = n1/n2;
	float n12rat2 = n12rat*n12rat;

	//takes care of trig
	dz = sqrt(1-n12rat2*(1-dz*dz));

	// 	if (dz != dz) return false;//total internal reflection

	//simpler expression than I expected
	dx *= n12rat;
	dy *= n12rat;

	return true;
}
bool DircBaseSim::x_wedge_coerce_check(\
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

	//This is bad and encapsulation breaking, but much faster and reduces code duplication
	if (midLineMode == true)
	{
		if (ty < wedgeHeight)
		{
			tdt = (wedgeHeight - ty)/dy;

			if (tdt > dt)
			{
				//does not go above wedge
				tdt = dt;
			}
			//only randomize if bouncing in infinitely thin wall
			//dx = midLineWedgeWallFlip*fabs(dx);
			//dx = midLineWedgeWallFlip*fabs(dx);
			//x starts at 0
			x = 0;
			midLineWedgeWallFlip *= -1;
		}
		else
		{
			//does not bounce
			tdt = dt;
			x += dx*tdt;
		}

		x += dx*(dt - tdt);
		y += dy*dt;
		z += dz*dt;
		//ignores possibility of bouncing out - only care about that below window anyway
		return true;
	}

	while (true) {
		//Get time of this bounce based on x direction
		if (dx < 0) {
			tdt = x + barWidth/2;
			tdt /= -dx;
			//coming from right wall
		} else {
			tdt = -x + barWidth/2 - wedgeWidthOff;
			tdt /= dx;
			//coming from left wall
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
			//printf("in wedge: x: %12.04f y: %12.04f z: %12.04f\n", x + tdt*dx, ty, z + barDepth/2);
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
		wedge_bounces++;
	}
	//Finish last leg of trip:
	x += dx*tdt;
	y = barLength/2 + ty + dy*tdt;//Note that y starts at zero (bar middle), but ty starts at bar top
	z += dt*dz;//Implemented for convenience - z direction irrelevant to window calculation
	//printf("%12.04f\n",z);

	return true;//not out the window for the whole trip
}

void DircBaseSim::plane_reflect(\
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
		float offang /*=0*/) 
{
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
float DircBaseSim::get_z_intercept(\
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
float DircBaseSim::get_intercept_plane(\
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
float DircBaseSim::sgn(float val) {
	//stolen from http://stackoverflow.com/questions/1903954/is-there-a-standard-sign-function-signum-sgn-in-c-c
	return (0 < val) - (val < 0);
}
void DircBaseSim::fill_bounces_vecs(\
                std::vector<int> &fxbounces,\
                std::vector<int> &fzbounces,\
                std::vector<int> &fxdirbounces,\
                std::vector<int> &fzdirbounces,\
                std::vector<int> &fxindirbounces,\
                std::vector<int> &fzindirbounces)
{
	fxbounces.clear();
	fzbounces.clear();
	fxdirbounces.clear();
	fzdirbounces.clear();
	fxindirbounces.clear();
	fzindirbounces.clear();


	for (unsigned int i = 0; i < x_bounces.size(); i++)
	{
		fxbounces.push_back(x_bounces[i]);
	}
	for (unsigned int i = 0; i < z_bounces.size(); i++)
	{
		fzbounces.push_back(z_bounces[i]);
	}
	for (unsigned int i = 0; i < x_direct_bounces.size(); i++)
	{
		fxdirbounces.push_back(x_direct_bounces[i]);
	}
	for (unsigned int i = 0; i < z_direct_bounces.size(); i++)
	{
		fzdirbounces.push_back(z_direct_bounces[i]);
	}
	for (unsigned int i = 0; i < x_indirect_bounces.size(); i++)
	{
		fxindirbounces.push_back(x_indirect_bounces[i]);
	}
	for (unsigned int i = 0; i < z_indirect_bounces.size(); i++)
	{
		fzindirbounces.push_back(z_indirect_bounces[i]);
	}

}

