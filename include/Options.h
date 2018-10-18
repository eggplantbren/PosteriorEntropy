#ifndef TransitInfo_Options_h
#define TransitInfo_Options_h

#include <cstdlib>

namespace TransitInfo
{

// Number of data points
static constexpr size_t N = 101;

// Time range of data points
static constexpr double t_min = 0.0;
static constexpr double t_max = 1.0;
static constexpr double t_range = t_max - t_min;
static constexpr double dt = t_range/(N - 1);

} // namespace TransitInfo

#endif

