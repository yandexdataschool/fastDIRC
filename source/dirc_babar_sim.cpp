#include <vector>

#include "../include/dirc_babar_sim.h"
#include "../include/dirc_point.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
//
#include <stdlib.h>
#include <math.h>

#include <TRandom3.h>



DircBaBarSim::DircBaBarSim(
		int rand_seed /*=4357*/,\
		double isens_r/*=540.66*/, \
		double isens_subtend_angle/*=52.4*/, \
		double ibar_length/*=4900*/,\
		double ibar_width/*=35*/,\
		double ibar_depth/*=17.25*/,
		double iupper_wedge_top/*=178.6*/) : 
			DircBaseSim(
				rand_seed,\
				ibar_length,\
				ibar_width,\
				ibar_depth) {

//Attempts to simulate the babar dirc - I tried to take the parameters from the NIM, but there may be errors, this has not been reviewed extensively

	printf("BarLWD: %12.04f %12.04f %12.04f\n",barLength,barWidth,barDepth);
	sens_subtend_angle = isens_subtend_angle;
	sens_r = isens_r;
	storeOpticalAngles = false;

	//only used for checking collision

	sensCylMinZ = -sens_r*sin(sens_subtend_angle/57.3);
	sensCylY = barLength/2 + 1174.0 - sens_r;
	sensCylZ = 0;

	sidemirror_xl = -1000000;
	sidemirror_xr = 1000000;
	sidemirror_reflectivity = .9;

	//Take average
	quartzIndex = 1.47;
	liquidIndex = 1.47;
	quartzLiquidY = upperWedgeBottom;

	num_QE = 31;
	min_QE = 300;
	max_QE = 600;
	sep_QE = (max_QE - min_QE)/(num_QE - 1);

	double t_QE[31] = {\
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
	//This removes the upper wedge
	upperWedgeTop = wedgeHeight + windowThickness;
}
double DircBaBarSim::get_sens_r()
{
	return sens_r;
}
double DircBaBarSim::get_sens_subtend_angle()
{
	return sens_subtend_angle;
}
const double DircBaBarSim::get_cerenkov_angle_rand(double beta, double additional_spread, double &wavelength) {
        //May be slow enough to consider approximating in distribution generation
        double out_ang = 0;
        double tmp_lam = 0;
        double tmp_QE_val;
        double above_ind;
        int ind_QE;
        double n_lam;

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

void DircBaBarSim::build_readout_box()
{
	//NA
	//build_system();
	//still rebuild wedge and bars when this is called:
}
void DircBaBarSim::set_pmt_plane_zs(double imin, double imax)
{
	//NA
}
void DircBaBarSim::set_large_mirror_zs(double imin, double imax)
{
	//NA
}
void DircBaBarSim::set_sidemirror_reflectivity(double isr) 
{
	sidemirror_reflectivity = isr;
}
void DircBaBarSim::set_foc_mirror_r(double ifoc_r) {
	//Does not apply
}
void DircBaBarSim::set_focmirror_nonuniformity(double nonuni_deg) {
	//Does not apply
}
void DircBaBarSim::sidemirror_reflect_points(std::vector<dirc_point> &points) {
	double tmpx = 0;
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
void DircBaBarSim::sidemirror_reflect_point(dirc_point &point) {
	double tmpx = 0;
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
void DircBaBarSim::set_sidemirror(double ixr, double ixl) {
	sidemirror_xr = ixr;
	sidemirror_xl = ixl;
}
void DircBaBarSim::set_three_seg_mirror(bool itsm) {
	//Not applicable
	//Ideally, the inheritance structure would be rewritten such that 
	//these do not need to be implemented, but it's hard to have them 
	//accesable for geometry changes AND have a flexible geometry
}
void DircBaBarSim::set_pmt_offset(double r) {
	//positive r increases path length
	build_readout_box();
}
void DircBaBarSim::set_liquid_absorbtion(double iabs) {
	liquidAbsorbtion = iabs;
}
std::vector<double> DircBaBarSim::get_dist_traveled() {
	return dist_traveled;
}
void DircBaBarSim::set_store_traveled(bool sst/* = true*/) {
	store_traveled = sst;
}
//Liquid Index is the same as the upper quartz wedge - the volumes are connected

//Note that these Not Applicable functions are implemented for ease of integrating with dircfit.cpp
//They do not need to be implemented if you are writing your own geometry and driver
void DircBaBarSim::set_focus_mirror_angle(double ang,double yang, double zang) {
	//no focus mirror angle
}
void DircBaBarSim::set_pmt_angle(double ang) {
	//no PMT angles
}
void DircBaBarSim::set_store_optical_angles(bool istore)
{
	storeOpticalAngles = istore;
}
std::vector<double> DircBaBarSim::get_focus_photon_angles()
{
	//Not Applicable
	std::vector<double> rval;
	return rval;
}
std::vector<double> DircBaBarSim::get_large_flat_photon_angles()
{
	//Not Applicable
	std::vector<double> rval;
	return rval;
}
std::vector<double> DircBaBarSim::get_side_photon_angles()
{
	//Not Applicable
	std::vector<double> rval;
	return rval;
}
void DircBaBarSim::warp_sens_cyl(
	dirc_point &out_val,\
	double &mm_index,\
	double &x,\
	double &y,\
	double &z,\
	double dx,\
	double dy,\
	double dz)
{
	//Assumes that we start at the top of the wedge
	
	//intersects and reflects the sensitive cylinder
        //Pretending the cylinder has no Nx for intersection purposes
        //Should be a valid approximation, and saves a ton of speed

        double dydz_norm = sqrt(dz*dz+dy*dy);
        double dy_norm = dy/dydz_norm;
        double dz_norm = dz/dydz_norm;

        //there's gotta be a faster way to do all of this
        //different than plane intercept D.  Took this algorithm from internet (wolfram Mathworld)
        //Parametric intercept sucks
        double D = (z - sensCylZ)*dy_norm - (y - sensCylY)*dz_norm;
        double detD = sqrt(sens_r*sens_r - D*D);

        //dy > 0 always (or should be.  Remove sgn and fabs function calls after confirming)
        double zrel = D*dy_norm + sgn(dy_norm)*dz_norm*detD;
        double yrel = -D*dz_norm + fabs(dy_norm)*detD;

	//printf("D and detD: %12.04f %12.04f dy_norm and dz_norm: %12.04f %12.04f\n",D,detD,dy_norm,dz_norm);
        double newy = yrel + sensCylY;
        double newz = zrel + sensCylZ;

	if (newz > 0)
	{
		newz = 2*sensCylZ - newz;
	}
	if (newz < sensCylMinZ)
	{
		//double plane_theta = atan2(-zrel,yrel);
		//printf("minZ: %12.04f\n",sensCylMinZ);
		//printf("failed: z and y: %12.04f %12.04f dz and dy: %12.04f %12.04f zrel and yrel: %12.04f %12.04f atan: %12.04f\n",z,y,dz,dy,-zrel,yrel,plane_theta);
		z = 1337;
		out_val.t = -1337;
		return;
	}



        double ydiff = newy - y;
        double zdiff = newz - z;

	double mm_traveled = sqrt((ydiff*ydiff+zdiff*zdiff))/dydz_norm;// total mm
	double plane_theta = atan2(-zrel,yrel);
	//printf("z and y: %12.04f %12.04f dz and dy: %12.04f %12.04f zrel and yrel: %12.04f %12.04f atan: %12.04f\n",z,y,dz,dy,-zrel,yrel,plane_theta);
	out_val.x = x + dx*mm_traveled;
	out_val.y = sens_r*plane_theta;
	out_val.t = mm_index + liquidIndex*mm_traveled;
}
void DircBaBarSim::warp_readout_box(
	dirc_point &out_val,\
	int particle_bar,\
	double &mm_index,\
	double &x,\
	double &y,\
	double &z,\
	double &dx,\
	double &dy,\
	double &dz)
{
	double c_mm_ns = 300;
	if (z > 0) {
		out_val.t = -1337;
		return;
	}
	warp_sens_cyl(\
		out_val,\
		mm_index,\
		x,\
		y,\
		z,\
		dx,\
		dy,\
		dz);
	//check absorbtion
	if (!(absorbtion_mc(dx,dy))) {
		out_val.t = -1337;
		return;
	}

	if (out_val.t > 0)
	{
		//printf("wsc: %12.04f %12.04f %12.04f\n",out_val.x,out_val.y,out_val.t);
	}
	out_val.t = mm_index/(c_mm_ns);
	out_val.x += get_bar_offset(particle_bar);
	
	//Must reflect after offset....
	sidemirror_reflect_point(out_val);
}
bool DircBaBarSim::absorbtion_mc(double dx, double dy) {
	//True if the particle makes it
	//expects vector pointing after bounce
	//Magic number corresponding to minimal distancel traveled
	double min_dist = 650+56;

	//approximating here - assume yz distance is independent of where on mirror it hit
	//measuring side plots with a ruler seems to mean this is good to ~5-10%
	double approx_dist = min_dist*sqrt(dx*dx+dy*dy)/dy;

	if (store_traveled == true) {
		dist_traveled.push_back(approx_dist);
	}

	double prob_transmitted = exp(-liquidAbsorbtion*approx_dist);
	if (rand_gen->Uniform(0,1) < prob_transmitted) {
		return true;
	} else {
		return false;
	}
}

