#include "dirc_point.h"
#include "dirc_lut_enum.h"

#ifndef DIRC_GLUEX_LUT_ENUM
#define DIRC_GLUEX_LUT_ENUM
class DircGluexLUTEnum: public DircLUTEnum
{
private:
	float minx,maxx,miny,maxy,resx,resy;
	int x_n, y_n;
public:
	DircGluexLUTEnum(float iminx,\
                float imaxx,\
                float iresx,\
                float iminy,\
                float imaxy,\
                float iresy);

	int return_enum(dirc_point &pt);
};
#endif
