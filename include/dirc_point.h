#ifndef DIRC_POINT_H
#define DIRC_POINT_H
struct dirc_point
{
	float x;
	float y;
	float t;
	int updown;
	int last_wall_x;
	int wedge_before_interface;
	float weight;
	float init_phi; //internal validation only
};
#endif
