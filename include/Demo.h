#ifndef TransitInfo_Demo_h
#define TransitInfo_Demo_h

#include <ostream>
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

        // Noise-free model curve
        std::vector<double> mus;
        void compute_mus();

    public:

        Demo();
        void from_prior(InfoNest::RNG& rng);
        std::vector<double> simulate_data(InfoNest::RNG& rng) const;
        double log_likelihood(const std::vector<double>& ys) const;
        friend std::ostream& operator << (std::ostream& out,
                                          const Demo& demo);
};


/****************** IMPLEMENTATION FOLLOWS ********************/

Demo::Demo()
:mus(N)
{

}

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

    // Compute noise-free curve
    compute_mus();
}

void Demo::compute_mus()
{
    // Beginning and end of transit
    double t_start = tc - 0.5*width;
    double t_end   = tc + 0.5*width;

    double t;
    for(size_t i=0; i<N; ++i)
    {
        // Timestamp of the data point
        t = t_min + i*dt;

        if(t > t_start && t < t_end)
            mus[i] = -depth;
        else
            mus[i] = 0.0;
    }
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


double Demo::log_likelihood(const std::vector<double>& ys) const
{
    double var = sigma*sigma;
    double tau = 1.0/var;
    double C = -0.5*log(2*M_PI*var);

    double logL = 0.0;
    for(size_t i=0; i<N; ++i)
        logL += C - 0.5*pow(ys[i] - mus[i], 2)*tau;

    return logL;
}

std::ostream& operator << (std::ostream& out,
                                          const Demo& demo)
{
    out << demo.tc << ' ' << demo.depth << ' ';
    out << demo.width << ' ' << demo.sigma;
    return out;
}

std::ostream& operator << (std::ostream& out,
                                          const std::vector<double>& ys)
{
    for(size_t i=0; i<ys.size(); ++i)
        out << ys[i] << ' ';
    return out;
}



} // namespace


#endif

