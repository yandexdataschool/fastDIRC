#include "dirc_point.h"
#include "../include/dirc_digitizer.h"
#include <vector>
DircDigitizer::DircDigitizer(\
		float iminx,\
		float imaxx,\
		float iresx,\
		float iminy,\
		float imaxy,\
		float iresy,\
		float it_unc,\
		float it_bin_size)
{
	minx = iminx;
	maxx = imaxx;
	resx = iresx;
	miny = iminy;
	maxy = imaxy;
	resy = iresy;
	t_unc = it_unc;
	t_bin_size = it_bin_size;
	
	dig_rand = new TRandom3();
}
void DircDigitizer::digitize_point(dirc_point &pt)
{
		//overflow/underflow?
	float x = pt.x;
	float y = pt.y;
	int xdig = 1.0000001*round((x - minx)/resx);
	int ydig = 1.0000001*round((y - miny)/resy);
	
	float xout = resx*xdig + minx;
	float yout = resy*ydig + miny;
	
	if (x < minx)
	{
		xout = minx - resx/2;
	}
	else if (x > maxx)
	{
		xout = maxx + resx/2;
	}
	
	if (y < miny)
	{
		yout = miny - resy/2;
	}
	else if (y > maxy)
	{
		yout = maxy + resy/2;
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
}
void DircDigitizer::digitize_points(std::vector<dirc_point> &points)
{
	for (unsigned int i = 0; i < points.size(); i++)
	{
		//Wasting time Make it faster later TODO
		digitize_point(points[i]);
	}
}

	
