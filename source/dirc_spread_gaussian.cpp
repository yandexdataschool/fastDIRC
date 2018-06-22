#include "../include/dirc_spread_gaussian.h"

DircSpreadGaussian::DircSpreadGaussian(
	float isigma,
	const std::vector<dirc_point>& isupport,
	float ix_unc,
	float iy_unc,
	float it_unc): x_unc(ix_unc), y_unc(iy_unc), t_unc(it_unc) {
	sigma2 = isigma*isigma;
	sigma2inv = 1/sigma2;
	spread_func_norm = 1;
	spread_func_norm_inv=1/spread_func_norm;
	support_cutoff_radius2 = 5*sigma2;
	support_points.pts.reserve(isupport.size());
	for (auto& point: isupport) {
	    support_points.pts.push_back({
                        point.x / x_unc,
			point.y / y_unc,
			point.t / t_unc});
	}
	support_index = std::make_unique<kd_tree>(3, support_points,
            nanoflann::KDTreeSingleIndexAdaptorParams(
						      kd_max_leaf));
	support_index->buildIndex();
}

const float DircSpreadGaussian::get_log_likelihood(const std::vector<dirc_point>& inpoints) const {
    float rval = 0;
    // TODO(kazeevn) math for the below threshold cases is fishy

    // We expected to see no hits, and saw none. Fine!
    if ((inpoints.size() == 0) && (support_index->m_size == 0)) {
	return log(1.);
    }

    // We expected hits, but got none
    // The honest value = -inf
    if (inpoints.size() == 0) {
	return log(min_probability);
    }

    for (auto& point: inpoints) {
	float tprob = 0;
	const std::array<float, 3> scaled_coords({
		point.x / x_unc,
		    point.y / y_unc,
		    point.t / t_unc
		    });
	std::vector<std::pair<size_t, float> > ret_matches;
	support_index->radiusSearch(
				    scaled_coords.data(), support_cutoff_radius2, 
				    ret_matches, nanoflann::SearchParams());
		 
	for (auto& support_point: ret_matches) {
	    tprob += radius_spread_function(std::get<1>(support_point));
	}
	// TODO deal with normalization....
	rval += log(tprob * spread_func_norm_inv / support_index->m_size + min_probability);
    }
    rval -= log(inpoints.size());
    return rval;
}

