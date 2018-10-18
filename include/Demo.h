#ifndef TransitInfo_Demo_h
#define TransitInfo_Demo_h

#include "Options.h"
#include "InfoNest/cpp/RNG.h"

namespace TransitInfo
{

// A simple demo planet model specification.
class Demo
{
    private:

        // Just parameterise centre time, depth and width of the transit.
        double tc;
        double depth;
        double width;
        double sigma; // Noise sd.

    public:

        void from_prior(InfoNest::RNG& rng);
        std::vector<double> simulate_data(InfoNest::RNG& rng) const;
};


/****************** IMPLEMENTATION FOLLOWS ********************/

void Demo::from_prior(InfoNest::RNG& rng)
{
    // Uniform(t_min, t_max)
    tc = t_min + t_range*rng.rand();

    // Exponential(0.1)
    depth = -0.1*log(rng.rand());

    // Lognormal(median=0.01, sd_ln=1)
    width = exp(log(0.01) + rng.randn());

    // Lognormal(median=0.01, sd_ln=1)
    sigma = exp(log(0.01) + rng.randn());
}


std::vector<double> Demo::simulate_data(InfoNest::RNG& rng) const
{
    std::vector<double> ys(N);

    // Beginning and end of transit
    double t_start = tc - 0.5*width;
    double t_end   = tc + 0.5*width;

    double t;
    for(size_t i=0; i<N; ++i)
    {
        // Timestamp of the data point
        t = t_min + i*dt;

        if(t > t_start && t < t_end)
            ys[i] = -depth;
        else
            ys[i] = 0.0;

        // Add the noise
        ys[i] += sigma*rng.randn();
    }

    return ys;
}

} // namespace


#endif

