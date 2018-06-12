#include "dirc_point.h"
#include "dirc_digitizer.h"
#include <vector>

//Root is infecting more files :(
#include <TRandom3.h>

#ifndef DIRC_RECT_DIGITIZER
#define DIRC_RECT_DIGITIZER
class DircRectDigitizer : public DircDigitizer
{
private:
	float minx,maxx,miny,maxy;
	float resx, resy;
	float t_unc;
	float t_bin_size;
	
	TRandom3 dig_rand;
	
public:
	DircRectDigitizer(\
		float iminx,\
		float imaxx,\
		float iresx,\
		float iminy,\
		float imaxy,\
		float iresy,\
		float it_unc,\
		float it_bin_size);

	void digitize_point(dirc_point &pt);
	void digitize_points(std::vector<dirc_point> &points);
};
#endif
