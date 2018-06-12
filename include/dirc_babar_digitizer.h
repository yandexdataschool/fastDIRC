#include "dirc_point.h"
#include "dirc_digitizer.h"
#include <vector>

//Root is infecting more files :(
#include <TRandom3.h>

#ifndef DIRC_BABAR_DIGITIZER
#define DIRC_BABAR_DIGITIZER
class DircBaBarDigitizer : public DircDigitizer
{
private:
	float minx,maxx,miny,maxy;
	float pmt_r;
	float pmt_sep_x,pmt_sep_y;
	float second_grid_minx,second_grid_miny;
	float resx, resy;
	float t_unc;
	float t_bin_size;
	
	TRandom3 *dig_rand;
	
public:
	DircBaBarDigitizer(\
		float iminx,\
		float imaxx,\
		float iminy,\
		float imaxy,\
		float ipmt_r,\
		float it_unc,\
		float it_bin_size);

	void digitize_point(dirc_point &pt);
	//return negative time if point is lost
	void digitize_points(std::vector<dirc_point> &points);
};
#endif
