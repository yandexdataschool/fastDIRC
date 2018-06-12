#include "dirc_point.h"
#include "../include/dirc_gluex_lut_enum.h"
#include <vector>
DircGluexLUTEnum::DircGluexLUTEnum(\
		float iminx,\
		float imaxx,\
		float iresx,\
		float iminy,\
		float imaxy,\
		float iresy):DircLUTEnum()
{
	minx = iminx;
	maxx = imaxx;
	resx = iresx;
	miny = iminy;
	maxy = imaxy;
	resy = iresy;
	x_n = (maxx-minx)/resx + 1;
	y_n = (maxy-miny)/resy + 1;
	max_n = x_n*y_n;
}
int DircGluexLUTEnum::return_enum(dirc_point &pt)
{
	//overflow/underflow?
	float x = pt.x;
	float y = pt.y;
	int xdig = (x - minx)/resx;
	int ydig = (y - miny)/resy;
	
	if (x < minx || x > maxx || y < miny || y > maxy)
	{
		return -1;
	}

	//Independent integer for each point
	return x_n*ydig + xdig;
}
	
