#ifndef DIRC_SPREAD_GAUSSIAN_H
#define DIRC_SPREAD_GAUSSIAN_H
#include <vector>
#include <nanoflann.hpp>
#include <memory>
#include "dirc_point.h"

// nanoflann
struct PointCloud {
    struct Point
    {
	float x,y,z;
    };
    std::vector<Point> pts;

    // Must return the number of data points
    inline size_t kdtree_get_point_count() const { return pts.size(); }

    // Returns the dim'th component of the idx'th point in the class:
    // Since this is inlined and the "dim" argument is typically an immediate value, the
    //  "if/else's" are actually solved at compile time.
    inline float kdtree_get_pt(const size_t idx, int dim) const
    {
	if (dim == 0) return pts[idx].x;
	else if (dim == 1) return pts[idx].y;
	else return pts[idx].z;
    }

    // Optional bounding-box computation: return false to default to a
    //   standard bbox computation loop.  Return true if the BBOX was
    //   already computed by the class and returned in "bb" so it can
    //   be avoided to redo it again.  Look at bb.size() to find out
    //   the expected dimensionality (e.g. 2 or 3 for point clouds)
    template <class BBOX>
    bool kdtree_get_bbox(BBOX& /* bb */) const { return false; }
};

typedef nanoflann::KDTreeSingleIndexAdaptor<
     nanoflann::L2_Simple_Adaptor<float, PointCloud>,
     PointCloud,
     3 /* dim */
     > kd_tree;

class DircSpreadGaussian {
private:
	float x_unc, y_unc, t_unc;
	float spread_func_norm, spread_func_norm_inv;
	float lin_slope, r_trans, sigma2, sigma2inv,max_val;
	float get_weight(dirc_point inpoint);
	const float min_probability = 1e-5;
	const unsigned int kd_max_leaf = 1000;
	float support_cutoff_radius2;
	PointCloud support_points;
	std::unique_ptr<kd_tree> support_index;
	float const_ll_component;
public:
	DircSpreadGaussian(
		float isigma,
		const std::vector<dirc_point>& isupport,
		float x_unc,
		float y_unc,
		float t_unc);

	// Zero-cutoff is the caller's responsibility
	const inline float radius_spread_function(const float r2) const {
	    return exp(-r2*sigma2inv);
	};
	const float get_log_likelihood(const std::vector<dirc_point>& inpoints) const;
};
#endif
