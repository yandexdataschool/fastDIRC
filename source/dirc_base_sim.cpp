#include <vector>

#include "../include/dirc_base_sim.h"
#include "../include/dirc_point.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
//
#include <stdlib.h>
#include <math.h>

#include <TRandom3.h>



DircBaseSim::DircBaseSim(
		int rand_seed /*=4357*/,\
		double ibarLength/*=4900*/, \
		double ibarWidth/*=35*/, \
		double ibarDepth/*=17.25*/,
		double iupperWedgeTop/*=178.6*/) {

	barLength=ibarLength;
	barWidth=ibarWidth;
	barDepth=ibarDepth;

	wedgeWidthOff = 1.75;
	wedgeDepthOff = 10;
	wedgeDepthOff = 9.75;
	wedgeFarAngle = .006*57.3;
	wedgeCloseAngle = 30;
	wedgeWidth=barWidth - wedgeWidthOff;
	wedgeDepthHigh = 79;
	//wedgeDepthHigh = 79.53;
	wedgeHeight = 91;

	upperWedgeAngleStore = false;

	wedgeDepthHigh = wedgeDepthOff+barDepth+wedgeHeight*sin(wedgeCloseAngle/57.296);

	windowThickness = 9.6;

	upperWedgeGap = 20;

	upperWedgeTop = iupperWedgeTop;
	upperWedgeBottom = wedgeHeight + windowThickness + upperWedgeGap;
	upperWedgeHeight = upperWedgeTop - upperWedgeBottom + upperWedgeGap;
	upperWedgeDepthHigh = wedgeDepthHigh + (upperWedgeTop-wedgeHeight)*sin(wedgeCloseAngle/57.296);


	lowerWedgeExtensionZ = -wedgeDepthHigh - tan(wedgeCloseAngle/57.296)*(upperWedgeBottom - wedgeHeight);


	//wedgeCloseAngle -= .5;

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

	rand_gen = new TRandom3(rand_seed);

	box_angle_off_cval = 1;
	box_angle_off_sval = 0;


	num_transmittance = 36;
	min_transmittance = 300;
	max_transmittance = 660;
	sep_transmittance = (max_transmittance - min_transmittance)/(num_transmittance - 1);

	//Numbers from Baptise
	double tmp_quartz_transmittance[36] = {\
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

	midLineMode = false;
	midLineWedgeWallFlip = 1;

	store_bounces = false;
	x_bounces.clear();
	z_bounces.clear();
	x_direct_bounces.clear();
	z_direct_bounces.clear();
	x_indirect_bounces.clear();
	z_indirect_bounces.clear();

	useMoliere = true;
	moliereP = 1000000;//minimal scattering 

	build_system();
}
void DircBaseSim::set_store_bounces(bool isb)
{
	store_bounces = isb;
}
void DircBaseSim::rotate_2d(double &x, double &y, double cval, double sval) {
	//Standard rotatation allows precomputation of matrix elements
	//Only store one variable for speed
	//Sin should be actual sin (not negative sin)
	double tx = x;
	x = cval*tx - sval*y;
	y = sval*tx + cval*y;
}

void DircBaseSim::set_bar_box_angle(double ang) {
	//expect degrees
	//sets angle between readout box and bars - rotates angle coming out of bars
	box_angle_off_cval = cos(ang/57.3);
	box_angle_off_sval = sin(ang/57.3);
}
void DircBaseSim::set_bar_box_offsets(double x, double y, double z) {
	//expect degrees
	//sets angle between readout box and bars - rotates angle coming out of bars
	bar_box_xoff = x;
	bar_box_yoff = y;
	bar_box_zoff = z;
}
std::vector<double> DircBaseSim::get_dist_traveled() {
	return dist_traveled;
}
void DircBaseSim::set_liquid_index(double li) {
	//Sets the index of the upper quartz wedge
	liquidIndex = li;
	//printf("Liquid index set: %12.04f\n",liquidIndex);

}
void DircBaseSim::set_upper_wedge_angle_diff(double rads, double rads_y) {
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
void DircBaseSim::set_wedge_mirror_rand(double ispread) {
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
		double spread_ang = wedgeCloseAngle;
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
double DircBaseSim::get_quartz_n(double lambda) {
	double B1,B2,B3,C1,C2,C3;
	B1 = 0.6961663;             // B1
	B2 = 0.4079426;             // B2
	B3 = 0.8974794;             // B3
	C1 = 0.0046791;             // C1
	C2 = 0.0135121;             // C2
	C3 = 97.9340025;          // C3
	double lam2;
	double n_lam;
	lam2 = lambda*lambda/1000000;

	n_lam = sqrt(1 + B1*lam2/(lam2-C1) + B2*lam2/(lam2-C2) + B3*lam2/(lam2-C3));

	return n_lam;
}
void DircBaseSim::set_use_quartz_n_for_liquid(bool iu)
{
	use_quartz_n_for_liquid = iu;
}
double DircBaseSim::get_liquid_n(double lambda) 
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
		double ln1 = 1.3786 - 0.1*lambda/1000;
		double ln2 = 1.3529 - 0.0353*lambda/1000;
		double rval = std::max(ln1,ln2);
		return rval;
	}
}
bool DircBaseSim::quartz_transmission_mc(double R, double lambda) {
	int ind_transmittance;
	double above_ind;
	double tmp_transmittance;
	double tmp_transmittance_constant;
	double trans_prob;

	ind_transmittance = (lambda - min_transmittance)/sep_transmittance;

	above_ind = lambda - (min_transmittance + sep_transmittance*ind_transmittance);

	//Simple linear interp between values.  5th order poly fit looked like it worked too
	tmp_transmittance = quartz_transmittance[ind_transmittance]*(sep_transmittance-above_ind)/sep_transmittance\
			    + quartz_transmittance[ind_transmittance+1]*above_ind/sep_transmittance;

	tmp_transmittance_constant = log(tmp_transmittance);



	trans_prob = exp(R/1000.*tmp_transmittance_constant);


	return (rand_gen->Uniform(0,1) < trans_prob);

}
double DircBaseSim::get_beta(double E, double m) {
	double gam = E/m;
	if (gam > 5) {
		//Approximate and possibly save time
		return 1-1/(2*gam*gam);
	} else {
		return sqrt(1-1/(gam*gam));
	}

}
double DircBaseSim::get_bar_offset(int bar)
{
	return fabs(bar)/bar*((150-0.5*(barWidth))+bar*(barWidth+.015));
}
int DircBaseSim::get_bar_from_x(double x)
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
std::vector<double> DircBaseSim::get_upper_wedge_incident()
{
	return upper_wedge_incident;
}
void DircBaseSim::sim_rand_n_photons(\
		std::vector<dirc_point> &out_points,\
		int n_photons, \
		double ckov_theta /*= 47*/, \
		double particle_bar /*=0*/, \
		double particle_x /*= 0*/, \
		double particle_y /*= 0*/, \
		double particle_t /*=0*/,\
		double particle_theta /*= 0*/, \
		double particle_phi /*= 0*/,\
		double phi_theta_unc /*= .0015*57.3*/,\
		double ckov_theta_unc /* = .0055*57.3*/,\
		double beta /* = -1*/) {

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
}
void DircBaseSim::sim_reg_n_photons(\
		std::vector<dirc_point> &out_points,\
		int n_photons_phi, \
		int n_photons_z,\
		double ckov_theta /*= 47*/, \
		double particle_bar /*=0*/, \
		double particle_x /*= 0*/, \
		double particle_y /*= 0*/, \
		double particle_t /*=0*/,\
		double particle_theta /*= 0*/, \
		double particle_phi /*= 0*/,\
		double phi_theta_unc /*= 0*/,\
		double ckov_theta_unc /* = 0*/,\
		double beta /* = -1*/) {
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
void DircBaseSim::test_from_wedge_top(\
		std::vector<dirc_point> &ovals,\
		int n_photons, \
		double particle_bar /*= 1*/, \
		double particle_x /*= 0*/, \
		double phot_theta /*= 0*/, \
		double phot_phi /*= 0*/,\
		double theta_unc, /*= 0*/
		double phi_unc /* = 0*/,\
		double overall_theta) {

	ovals.clear();
	//Note that Theta and Phi are defined along the bar axis, not as they are elsewhere
	//negative bar number is - x axis?
	int numPhots = n_photons;

	double sourcez = -barDepth/2;
	sourcez = -wedgeDepthHigh/2;
	sourcez = 0;
	double sourcey = barLength/2 + upperWedgeTop + .1;
	double sourcex = particle_x;

	double temit, randPhi;

	double x,y,z,dx,dy,dz;

	double mm_index = 0;

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
void DircBaseSim::bar_box_interface(\
		double &x,\
		double &y,\
		double &z,\
		double &dx,\
		double &dy,\
		double &dz)
{
	x += bar_box_xoff;
	y += bar_box_yoff;
	z += bar_box_zoff;
	rotate_2d(dy,dz,box_angle_off_cval,box_angle_off_sval);
}
void DircBaseSim::sim_lut_points(\
		std::vector<dirc_point> &ovals,\
		std::vector<double> &phis,\
		std::vector<double> &thetas,\
		int n_photons, \
		double particle_bar /*= 0*/){

	ovals.clear();
	phis.clear();
	thetas.clear();


	double x,y,z,dx,dy,dz;

	double randPhi;
	double randTheta;	

	double mm_index = 0;

//	double c_mm_ns = 300;

	double saveGeneralQuartzIndex = quartzIndex;
	double saveGeneralLiquidIndex = liquidIndex;
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
void DircBaseSim::fill_rand_phi(\
		std::vector<dirc_point> &ovals,\
		int n_photons, \
		double ckov_theta /*= 47*/, \
		double particle_bar /*= 0*/, \
		double particle_x /*= 0*/, \
		double particle_y /*= 0*/, \
		double particle_t /*= 0*/, \
		double particle_theta /*= 0*/, \
		double particle_phi /*= 0*/,\
		double phi_theta_unc, /*= .0015*57.3*/
		double ckov_theta_unc /* = .0055*57.3*/,\
		double beta/* = -1*/) {


	double emitAngle = ckov_theta;
	double particleTheta = particle_theta + rand_gen->Gaus(0,phi_theta_unc);
	double particlePhi = particle_phi + rand_gen->Gaus(0,phi_theta_unc);
	int numPhots = n_photons/cos(particle_theta/57.3);

	double randPhi;

	double temit, rand_add;
	double wavelength = 400;

	double x,y,z,dx,dy,dz;

	double cos_ptheta = cos(particleTheta/57.3);
	double sin_ptheta = sin(particleTheta/57.3);
	double cos_pphi = cos(particlePhi/57.3);
	double sin_pphi = sin(particlePhi/57.3);

	double mm_index = 0;

	int tmp_updown = 0;

	double saveGeneralQuartzIndex = quartzIndex;
	double saveGeneralLiquidIndex = liquidIndex;

	std::vector<dirc_base_sim_tracking_step> track_steps;
	double dist_traveled = -1;
	//hardcode step length to 1mm for now
	double step_length = 1;
	//note: These are general tracking vectors - you can implement your own MC with them.  Each refers to position and direction at the start of the step
	if (useMoliere == true)
	{
		fill_moliere_tracking_steps(\
				track_steps,\
				dist_traveled,\
				step_length,\
				particleTheta/57.3,\
				particlePhi/57.3,\
				particle_x,\
				particle_y,\
				-barDepth);
	}
	else
	{
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
		step_length = barDepth/cos_ptheta;

		
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

	}

	double track_loc = -1;

	//define numbers to linearly interp
	int low_ind = 0;
	double above_ind = 0;

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
		ovals.push_back(out_val);
	}
	quartzIndex = saveGeneralQuartzIndex;
	liquidIndex = saveGeneralLiquidIndex;
}
std::vector<std::pair<double,double> > DircBaseSim::get_refraction_rand_phi(\
		std::vector<double> &before_interface,\
		std::vector<double> &after_interface,\
		std::vector<double> &pmt_incidence,\
		int n_photons, \
		double ckov_theta /*= 47*/, \
		double particle_x /*= 0*/, \
		double particle_y /*= 0*/, \
		double particle_theta /*= 0*/, \
		double particle_phi /*= 0*/,\
		double phi_theta_unc /*= .0015*57.3*/,\
		double ckov_theta_unc /* = .0055*57.3*/,\
		double beta/* = -1*/) {
	//returns theta versus cerenkov phi_theta_unc

	std::vector<std::pair<double,double> > rval;

	before_interface.clear();
	after_interface.clear();
	refraction_before.clear();
	refraction_after.clear();
	pmt_incidence.clear();

	double emitAngle = ckov_theta;
	double particleTheta = particle_theta + rand_gen->Gaus(0,phi_theta_unc);
	double particlePhi = particle_phi + rand_gen->Gaus(0,phi_theta_unc);

	int numPhots = n_photons/cos(particle_theta/57.3);

	// 	double sourcez = -sDepth;
	double sourcey = particle_y-barDepth*tan(particleTheta/57.3);
	double sourcex = particle_x;
	double tsy = sourcey;
	double tsx = sourcex;
	sourcey = tsy*cos(particlePhi/57.3)-tsx*sin(particlePhi/57.3);
	sourcex = tsy*sin(particlePhi/57.3)+tsx*cos(particlePhi/57.3);


	double sourceOff,randPhi;

	double temit, rand_add;
	double wavelength = 400;

	double x,y,z,dx,dy,dz;

	double cos_ptheta = cos(particle_theta/57.3);
	double sin_ptheta = sin(particle_theta/57.3);
	double cos_pphi = cos(particle_phi/57.3);
	double sin_pphi = sin(particle_phi/57.3);

	double mm_index = 0;
	double c_mm_ns = 300;

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
		std::pair<double, double> add_val;

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
		std::vector<dirc_point> &ovals,
		int n_photons_phi,
		int n_photons_z,
		double ckov_theta /*= 47*/,
		double particle_bar /*= 0*/,
		double particle_x /*= 0*/,
		double particle_y /*= 0*/,
		double particle_t /*= 0*/,
		double particle_theta /*= 0*/,
		double particle_phi /*= 0*/,
		double phi_theta_unc /*= 0*/,
		double ckov_theta_unc /* = 0*/,
		double beta /* = -1*/)
{
	double sDepth = .95*barDepth;
	double emitAngle = ckov_theta;

	int tmp_updown = 0;

	double sourceOff,regPhi;

	double temit;
	double rand_add;
	double wavelength = 400;

	double x,y,z,dx,dy,dz;

	double cos_ptheta = cos(particle_theta/57.3);
	double sin_ptheta = sin(particle_theta/57.3);
	double cos_pphi = cos(particle_phi/57.3);
	double sin_pphi = sin(particle_phi/57.3);

	double mm_index = 0;

	double sin_emit;
	double cos_emit;
	double sin_regphi;
	double cos_regphi;

	int adj_n_photons_phi = n_photons_phi/cos(particle_theta/57.3);
	double saveGeneralQuartzIndex = quartzIndex;
	double saveGeneralLiquidIndex = liquidIndex;

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
			//rotate_2d(dy,dz,box_angle_off_cval,box_angle_off_sval);
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
			ovals.push_back(out_val);
		}
	}
	quartzIndex = saveGeneralQuartzIndex;
	liquidIndex = saveGeneralLiquidIndex;

}
bool DircBaseSim::track_single_photon(\
		dirc_point &out_val,\
		double emit_theta,\
		double emit_phi,\
		double particle_theta,\
		double particle_phi,\
		double particle_x,\
		double particle_y,\
		double particle_z,\
		double particle_t,\
		int particle_bar)
{
	double x,y,z,dx,dy,dz;
	double mm_index;

	double temit = emit_theta;
	double regPhi = emit_phi;

	double cos_emit, sin_emit;
	double cos_regphi, sin_regphi;
	int tmp_updown = 0;

	double cos_ptheta = cos(particle_theta/57.3);
	double sin_ptheta = sin(particle_theta/57.3);
	double cos_pphi = cos(particle_phi/57.3);
	double sin_pphi = sin(particle_phi/57.3);

	double sourceOff = particle_z + barDepth;

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
		double emit_theta,\
		double emit_phi,\
		double particle_theta,\
		double particle_phi,\
		double particle_x,\
		double particle_y,\
		double particle_z,\
		double particle_t,\
		int particle_bar,\
		double z_at_top /*=1*/)
{
	//I'm being bad and breaking encapsulation here, but it's for the greater good
	//Think of the children
	double saveBarWidth = barWidth;
	double saveWedgeWidthOff = wedgeWidthOff;

	midLineMode = true;
	midLineMode = false;
	double width_fraction = 1;
	barWidth *= width_fraction;
	wedgeWidthOff *= width_fraction;
	wedgeWidthOff = 0;
	build_system();

	double x,y,z,dx,dy,dz;
	double mm_index;

	double temit = emit_theta;
	double regPhi = emit_phi;

	double cos_emit, sin_emit;
	double cos_regphi, sin_regphi;
	int tmp_updown = 0;

	double cos_ptheta = cos(particle_theta/57.3);
	double sin_ptheta = sin(particle_theta/57.3);
	double cos_pphi = cos(particle_phi/57.3);
	double sin_pphi = sin(particle_phi/57.3);

	double sourceOff = particle_z + barDepth;

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
	double target_z = -17.25/2;

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
		double particle_beta,\
		double emit_phi,\
		double particle_theta,\
		double particle_phi,\
		double particle_x,\
		double particle_y,\
		double particle_z,\
		double particle_t,\
		int particle_bar)
{
	double x,y,z,dx,dy,dz;
	double mm_index;

	double wavelength = 400;

	double temit = get_cerenkov_angle_rand(particle_beta,0,wavelength);
	double regPhi = emit_phi;

	double cos_emit, sin_emit;
	double cos_regphi, sin_regphi;
	int tmp_updown = 0;

	double cos_ptheta = cos(particle_theta/57.3);
	double sin_ptheta = sin(particle_theta/57.3);
	double cos_pphi = cos(particle_phi/57.3);
	double sin_pphi = sin(particle_phi/57.3);

	double sourceOff = particle_z + barDepth;

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
		double emit_theta,\
		double particle_theta,\
		double particle_phi,\
		double particle_x,\
		double particle_y,\
		double particle_z,\
		double particle_t,\
		int particle_bar,\
		double z_at_top /*=1*/)
{
	//I'm being bad and breaking encapsulation here, but it's for the greater good
	//Think of the children
	double saveBarWidth = barWidth;
	double saveWedgeWidthOff = wedgeWidthOff;

	left_vals.clear();
	right_vals.clear();

	midLineMode = true;
	double width_fraction = .001;
	barWidth *= width_fraction;
	wedgeWidthOff *= width_fraction;
	wedgeWidthOff = 0;
	build_system();

	std::vector<double> mustache_phi;

	double x,y,z,dx,dy,dz;
	//righthand side versions - right hand not needed, just reuse
	//double x_r,y_r,z_r,dx_r,dy_r,dz_r;
	double mm_index;
//	double mm_index_r;

	double temit = emit_theta;
	double regPhi;

	double cos_emit, sin_emit;
	double cos_regphi, sin_regphi;
	int tmp_updown = 0;

	double cos_ptheta = cos(particle_theta/57.3);
	double sin_ptheta = sin(particle_theta/57.3);
	double cos_pphi = cos(particle_phi/57.3);
	double sin_pphi = sin(particle_phi/57.3);

	double sourceOff = particle_z + barDepth;

	double phi_inc = (2*3.14159265359)/points_per_side;
	regPhi = 0;
	double saveGeneralQuartzIndex = quartzIndex;
	double saveGeneralLiquidIndex = liquidIndex;
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


		double target_z = -17.25/2;

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
		double target_z = 0;

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
		double enhance_rot_fac = 3;
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
double DircBaseSim::warp_ray(\
		double &x,\
		double &y,\
		double &z,\
		double &dx,\
		double &dy,\
		double &dz,\
		double cos_critical_angle) {
	//Implemented to avoid total internal reflection computations
	//Expects ray to be "prerotated" by particle theta and phi - just propagates the vectors
	//Uses bar geometry.
	//Modifys the ray.
	//Be careful about x,y,and z - can't straight rotate that
	//returns distance traveled (mm) times index


	double rval = 0; //stores total distance;

	double delx, dely,delz;
	// 	double xzR;
	//Start the rays in the quartz for simulation reasons
	double grace_room = 0;

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

	double remainderx;
	// 	double remaindery;
	double lrx;
	// 	double lrz = 1;
	double remainderz = 0;
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



	//	double bar_front = -barDepth/2;
	//	double bar_front = -17.25/2;
	double bar_front = 0;
	if (nbouncesz % 2 == 1) {
		dz = -dz;
		z = bar_front-remainderz;
	} else {
		z = remainderz-barDepth;
	}

	return sqrt(rval)*quartzIndex;
}
double DircBaseSim::warp_wedge(
		double &x,
		double &y,
		double &z,
		double &dx,
		double &dy,
		double &dz) {
	//No Critical angle checking - shouldn't need it
	//I am abusing the x angle being the same in the bar and wedge here
	//may have to rotate after to the mirror
	//Starts y at 0

	wedge_bounces = 0;

	//must change after close wedge bounce
	double startz = z;
	double starty = y;
	double mm_index = 0;
	bool passed_interface = false;
	// 	double x0 = x;
	double dt = 0; //dummy variable representing parametric line
	double n_dot_v = 0; //dummy variable to save dot product
	double n_dot_v0 = 0;

	//deal with yz first

	//Can't ever not totally internal reflect	
	//force dz to be negative	
	//double wedgeCloseIncidentAngle = acos(-dx*wedgeClosePlaneNx - dy*wedgeClosePlaneNy + fabs(dz)*wedgeClosePlaneNz);
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
	double ty = dt*dy + y - barLength/2;


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
		double wally = starty + (y - starty)*(-startz/(z-startz));
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
			//double sz = z;
			double sdz = dz;
			double sdy = dy;
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
		double n1,\
		double n2,\
		double &dx,\
		double &dy,\
		double &dz) {
	// n1 is starting index of refraction, n2 is ending
	// Assume that it's going through the z plane
	// (this will change the ordering when called in our current geometry)

	//dz and dy are flipped, so this is really acos
	double n12rat = n1/n2;
	double n12rat2 = n12rat*n12rat;

	//takes care of trig
	dz = sqrt(1-n12rat2*(1-dz*dz));

	// 	if (dz != dz) return false;//total internal reflection

	//simpler expression than I expected
	dx *= n12rat;
	dy *= n12rat;

	return true;
}
bool DircBaseSim::x_wedge_coerce_check(\
		double &x,\
		double &y,\
		double &z,\
		double &dx,\
		double &dy,\
		double &dz,\
		double dt) {
	//already have dt, be sure to add it in when calling this

	//assumes x starts between the sides of the wedge
	//assumes dy/dx constant over the whole distance
	//also performs the index plane crossing

	//changes x and dx, correctly propagating them
	//propagates y and z as well
	//returns false if it goes out the window on the x step
	//Will not check for y and z directions (assumes dt corresponds to an intersection/reflection point)

	double cdt = 0;
	double tdt = 0;
	double ty = y - barLength/2;

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
void DircBaseSim::set_moliere_p(double ip)
{
	moliereP = ip;
}
void DircBaseSim::set_use_moliere(bool ium)
{
	useMoliere = ium;
}
double DircBaseSim::generate_cos_moliere_angle(\
		double rad_length)
{
	//Algorithm from Macro by Mike Williams 07/26/16
	//Magic numbers taken from limits of scattering
//	TF1 g1("g1","exp(-[0]*(1-x))",-1,1);
//	TF1 g2("g2","pow([0]-x,-2)",-1,1); 

	// .17 corresponds to rad_length in this limit
	//double theta0=(13.6/P)*sqrt(0.17)*(1+0.038*log(0.17)); 
	rad_length = .17;
	double theta0=(13.6/moliereP)*sqrt(rad_length)*(1+0.038*log(rad_length)); 
	double a = pow(theta0,-2);
	double u0 = 1-3*pow(theta0,2);
	double b = 1-pow(theta0,2);
//	g1.SetParameter(0,a);
//	g2.SetParameter(0,b);

	double c1 = a;
	double c2 = 1/(1/(3*pow(theta0,2))-1/2.);

	double g10 = exp(-a*(1-u0))*c1;
	double g20 = pow(b-u0,-2)*c2;
	double p = g20/(g10+g20);

	double max_g = 0;

	//binning chosen based on original ROOT histograms
	//actually undercounts forward scattering for this formulation
	//done this way to be able to simulate this in nonzero time
	double u = 1-.005/100000;
	double g1val = 0, g2val = 0;
	if(u >= u0) g1val = c1*exp(-a*(1-u));
	if(u <= u0) g2val = c2*pow(b-u,-2);
	double g = p*g1val + (1-p)*g2val;
	if(g < 0) g=0;
	//printf("%12.04f\n",g);
	max_g = std::max(g,max_g);


	while (true)
	{
		double u = rand_gen->Uniform(.99,1);
		double g1val = 0, g2val = 0;
		if(u >= u0) g1val = c1*exp(-a*(1-u));
		if(u <= u0) g2val = c2*pow(b-u,-2);
		double g = p*g1val + (1-p)*g2val;
		if(g < 0) g=0;
		//printf("%12.04f\n",g);
		if (rand_gen->Uniform(0,max_g) < g)
		{
			return u;
		}
	}
}
void DircBaseSim::fill_moliere_tracking_steps(\
		std::vector<dirc_base_sim_tracking_step> &rsteps,\
		double &travel_distance,\
		double step_length,\
		double start_theta,\
		double start_phi,\
		double start_x,\
		double start_y,\
		double start_z)
{
		rsteps.clear();
		travel_distance = 0;

		double max_travel = fabs(barDepth*2/cos(start_theta));
		
		double cur_x = start_x;
		double cur_y = start_y;
		double cur_z = start_z;
		double cur_theta = start_theta;
		double cur_phi = start_phi;

		//printf("Enter ctheta: %12.04f\n",cur_theta);
	
		double cur_sintheta = sin(cur_theta);
		double cur_costheta = cos(cur_theta);
		double cur_sinphi = sin(cur_phi);
		double cur_cosphi = cos(cur_phi);

		double sin_scatter_theta;
		double cos_scatter_theta;
		double scatter_phi;
		double sin_scatter_phi;
		double cos_scatter_phi;

		double sdx,sdy,sdz;


		while (cur_z < 0 && travel_distance < max_travel)
		{
			dirc_base_sim_tracking_step add_step;
			add_step.x = cur_x;
			add_step.y = cur_y;
			add_step.z = cur_z;
			add_step.sin_theta = cur_sintheta;
			add_step.cos_theta = cur_costheta;
			add_step.sin_phi = cur_sinphi;
			add_step.cos_phi = cur_cosphi;
		
	
			rsteps.push_back(add_step);

			cos_scatter_theta = generate_cos_moliere_angle(step_length);
			sin_scatter_theta = sqrt(1-cos_scatter_theta*cos_scatter_theta);

			scatter_phi = rand_gen->Uniform(0,2*3.14159265);

			cos_scatter_phi = cos(scatter_phi);
			//sin_scatter_phi = sqrt(1-cos_scatter_phi*cos_scatter_phi);			
			//compute together to save time in the future - not currently the limiting factor
			sin_scatter_phi = sin(scatter_phi);			

			sdx = step_length*sin_scatter_theta*cos_scatter_phi;
			sdy = step_length*sin_scatter_theta*sin_scatter_phi;
			sdz = step_length*cos_scatter_theta;

			sdx = 0;
			sdy = 0;
			sdz = 1;

			rotate_2d(sdz,sdy,cur_costheta,cur_sintheta);
			rotate_2d(sdy,sdx,cur_cosphi,cur_sinphi);

			//printf("inmol: %12.04f %12.04f %12.04f step_length: %12.04f cos_theta: %12.04f sdxyz:  %12.04f %12.04f %12.04f\n",cur_x,cur_y,cur_z,step_length,cos_scatter_theta,sdx,sdy,sdz);
			//printf("cst,cct,csp,ccp: %12.04f %12.04f %12.04f %12.04f step_length: %12.04f cos_theta: %12.04f sdxyz:  %12.04f %12.04f %12.04f\n",cur_sintheta,cur_costheta,cur_sinphi,cur_cosphi,step_length,cos_scatter_theta,sdx,sdy,sdz);
		
			//exits bar at this step, add it again
			if (cur_z + sdz > 0)
			{
				//Adjust to ensure the step 
				sdx = -sdx*cur_z/sdz*1.0001;
				sdy = -sdy*cur_z/sdz*1.0001;
				sdz = -sdz*cur_z/sdz*1.0001;
				travel_distance += -step_length*cur_z/sdz;
			}
			else
			{
				travel_distance += step_length;
			}
			
			cur_x += sdx;
			cur_y += sdy;
			cur_z += sdz;
			cur_theta = acos(sdz);
			//ordering below is correct - has to do with the rotate_2d coordinates above, which are left as they are for consistency with production
			//still produces the same phi (0 along x axis)
			cur_phi = atan2(sdx,sdy);

			cur_sintheta = sin(cur_theta);
			cur_costheta = cos(cur_theta);
			cur_sinphi = sin(cur_phi);
			cur_cosphi = cos(cur_phi);

			//same as condition above, finish with last step.
			if (cur_z > 0)
			{
				dirc_base_sim_tracking_step add_step;
				add_step.x = cur_x;
				add_step.y = cur_y;
				add_step.z = cur_z;
				add_step.sin_theta = cur_sintheta;
				add_step.cos_theta = cur_costheta;
				add_step.sin_phi = cur_sinphi;
				add_step.cos_phi = cur_cosphi;
				
				rsteps.push_back(add_step);
			}
		}
}

void DircBaseSim::plane_reflect(\
		double Nx,\
		double Ny,\
		double Nz,\
		double D,\
		double &x,\
		double &y,\
		double &z,\
		double &dx,\
		double &dy,\
		double &dz,\
		double &dt,\
		double offang /*=0*/) 
{
	//Propagate to the intercept and reflects off a plane
	//Could use this in the Wedge, but would loose some speed

	double n_dot_v = (Nx*dx + Ny*dy + Nz*dz);
	double n_dot_v0 = (Nx*x + Ny*y + Nz*z);

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
double DircBaseSim::get_z_intercept(\
		double Nx,\
		double Ny,\
		double Nz,\
		double D,\
		double x,\
		double y,\
		double z,\
		double dx,\
		double dy,\
		double dz) {
	//finds z intercept (without modifying x y and z)
	//Could use this in the Wedge, but would loose some speed

	double n_dot_v = (Nx*dx + Ny*dy + Nz*dz);
	double n_dot_v0 = (Nx*x + Ny*y + Nz*z);
	double dt;//optimized by math compression?
	dt = (D - n_dot_v0)/n_dot_v;

	return z + dt*dz;
}
double DircBaseSim::get_intercept_plane(\
		double Nx,\
		double Ny,\
		double Nz,\
		double D,\
		double &x,\
		double &y,\
		double &z,\
		double dx,\
		double dy,\
		double dz) {
	//Just intersects a plane (modifying x,y,z)
	//Could use this in the Wedge, but would loose some speed
	//Non returning verion for speed when needed?

	double n_dot_v = (Nx*dx + Ny*dy + Nz*dz);
	double n_dot_v0 = (Nx*x + Ny*y + Nz*z);
	double dt;//optimized by math compression?
	dt = (D - n_dot_v0)/n_dot_v;

	x += dx*dt;
	y += dy*dt;
	z += dz*dt;

	return dt;
}
double DircBaseSim::sgn(double val) {
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

