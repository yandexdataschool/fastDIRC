#include "dirc_point.h"
#include <utility>
#include <vector> 

#ifndef DIRC_LUT
#define DIRC_LUT
struct lut_entry
{
	float phi;
	float theta;
	float time;
};

class DircLUT
{
	private:

		DircLUTEnum* pt_to_ind;

		std::vector<std::vector<lut_entry> > lu_table;

		//Needed?  Private?
		std::vector<lut_entry> get_base_phi_theta(dirc_point pt);
	public:
		DircLUT(DircLUTEnum* ipt_to_ind);

		void add_table_pt(dirc_point pt, float phi, float theta);

		void get_ckov_theta_all(std::vector<float> &rval, std::vector<float> &ret_dt, std::vector<dirc_point> pt,float inc_phi, float inc_theta, float inc_y);
		//Needed?  Private?
		void get_base_phi_theta_all(std::vector<lut_entry> &rval, std::vector<dirc_point> pts);
	/*
		void get_ckov_theta_single_oval_cut(std::vector<float> &rval, \
				std::vector<float> &ret_dt, \
				std::vector<dirc_point> pts, \
				float inc_phi, \
				float inc_theta, \
				float inc_y, \
				float center_ang, \
				float center_ang_spread,\
				float center_time_spread);
	*/
		void get_chromcorr_m_b_single_oval_cut(
				std::vector<dirc_point> pts, \
				float inc_phi, \
				float inc_theta, \
				float inc_y, \
				float center_ang, \
				float center_ang_spread_sq,\
				float time_spread_sq,\
				float expected_angle,\
				float &m_direct,\
				float &b_direct,\
				float &m_indirect,\
				float &b_indirect);
		void get_ckov_theta_single_oval_cut(
				std::vector<float> &rval, \
				std::vector<float> &ret_dt_dl, \
				std::vector<dirc_point> pts, \
				float inc_phi, \
				float inc_theta, \
				float inc_y, \
				float center_ang, \
				float center_ang_spread_sq,\
				float time_spread_sq,\
				float m_direct=0,\
				float b_direct=0,\
				float m_indirect=0,\
				float b_indirect=0);



};
#endif
