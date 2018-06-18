#ifndef DIRC_RECT_DIGITIZER
#define DIRC_RECT_DIGITIZER
#include <vector>
#include <TRandom3.h>
#include "dirc_point.h"

class DircRectDigitizer {
private:
	float minx,maxx,miny,maxy;
	float resx, resy;
	float t_unc;
	float t_bin_size;
	std::unique_ptr<TRandom3> dig_rand;
public:
	DircRectDigitizer(
		float iminx,
		float imaxx,
		float iresx,
		float iminy,
		float imaxy,
		float iresy,
		float it_unc,
		float it_bin_size);

	void digitize_point(dirc_point &pt);
	void digitize_points(std::vector<dirc_point> &points);
};
#endif
