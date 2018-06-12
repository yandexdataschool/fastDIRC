#include "dirc_point.h"
#include "../include/dirc_babar_digitizer.h"
#include <vector>
DircBaBarDigitizer::DircBaBarDigitizer(\
		float iminx,\
		float imaxx,\
		float iminy,\
		float imaxy,\
		float ipmt_r,\
		float it_unc,\
		float it_bin_size):DircDigitizer()
{
//Assumes a set of Hexagonally packed circles and rounds to the center that was hit
//Rejects things outside the circle
	minx = iminx;
	maxx = imaxx;
	pmt_r = ipmt_r;
	miny = iminy;
	maxy = imaxy;

	pmt_sep_x = 2*pmt_r;
	pmt_sep_y = 2*sqrt(3)*pmt_r;
	
	second_grid_minx = minx + pmt_sep_x/2;
	second_grid_miny = miny + pmt_sep_y/2;

	t_unc = it_unc;
	t_bin_size = it_bin_size;
	
	dig_rand = new TRandom3();
}
void DircBaBarDigitizer::digitize_point(dirc_point &pt)
{
		//overflow/underflow?
	float xout,yout,tout;
	float x = pt.x;
	float y = pt.y;
	int xdig_1 = 1.0000001*round((x - minx)/pmt_sep_x);
	int ydig_1 = 1.0000001*round((y - miny)/pmt_sep_y);
	int xdig_2 = 1.0000001*round((x - second_grid_minx)/pmt_sep_x);
	int ydig_2 = 1.0000001*round((y - second_grid_miny)/pmt_sep_y);

	float xd1 = pmt_sep_x*xdig_1 + minx;
	float xd2 = pmt_sep_x*xdig_2 + second_grid_minx;
	float yd1 = pmt_sep_y*ydig_1 + miny;
	float yd2 = pmt_sep_y*ydig_2 + second_grid_miny;
	
	float dxd1 = x - xd1;
	float dyd1 = y - yd1;
	float dxd2 = x - xd2;
	float dyd2 = y - yd2;

	float dist1 = sqrt(dxd1*dxd1 + dyd1*dyd1);
	float dist2 = sqrt(dxd2*dxd2 + dyd2*dyd2);
	
	float mindist = std::min(dist1,dist2);

	//printf("x,ydig1: %5d %5d x,y,dig2: %5d %5d x,yd1: %12.04f %12.04f x,yd2: %12.04f %12.04f mindist: %12.04f %5d\n",xdig_1,ydig_1,xdig_2,ydig_2,xd1,yd1,xd2,yd2,mindist,ttat);
	
	if (mindist < pmt_r)
	{
		if (dist1 < dist2)
		{
			xout = xd1;
			yout = yd1;
		}
		else
		{
			xout = xd2;
			yout = yd2;
		}
	}
	else
	{
		xout = minx;
		yout = miny;
		tout = -1;
	}
	
	if (x < minx)
	{
		xout = minx - resx/2;
		tout = -1;
	}
	else if (x > maxx)
	{
		xout = maxx + resx/2;
		tout = -1;
	}
	
	if (y < miny)
	{
		yout = miny - resy/2;
		tout = -1;
	}
	else if (y > maxy)
	{
		yout = maxy + resy/2;
		tout = -1;
	}
	
	pt.x = xout;
	pt.y = yout;
	pt.t += dig_rand->Gaus(0,t_unc);
	if (fabs(t_bin_size) > t_unc/100)
	{
		//Don't bother with binning if it's small
		//perhaps slow - a hard cut or a constant variable could be better here
		int tmp_t_bin = pt.t/(t_bin_size);
		pt.t = tmp_t_bin*t_bin_size + t_bin_size/2;
	}
	if (tout < 0) pt.t = tout;
}
void DircBaBarDigitizer::digitize_points(std::vector<dirc_point> &points)
{
	for (unsigned int i = 0; i < points.size(); i++)
	{
		//Wasting time Make it faster later TODO
		//printf("bef xyz: %12.04f %12.04f %12.04f minx,y: %12.04f %12.04f\n",points[i].x,points[i].y,points[i].t,minx,miny);
		digitize_point(points[i]);
		//printf("pts xyz: %12.04f %12.04f %12.04f\n",points[i].x,points[i].y,points[i].t);
		if (points[i].t < 0)
		{
			//Sinful way of dealing with the deletion, sorry
			points.erase(points.begin()+i);
			i--;
		}
	}
}

	
