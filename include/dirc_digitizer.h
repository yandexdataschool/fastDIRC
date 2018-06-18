#ifndef DIRC_DIGITIZER
#define DIRC_DIGITIZER
#include <vector>
#include <memory>
#include <TRandom3.h>
#include "dirc_point.h"

class DircDigitizer {
public:
    virtual void digitize_point(dirc_point &pt) = 0;
    virtual void digitize_points(std::vector<dirc_point> &points) = 0;
};
#endif
