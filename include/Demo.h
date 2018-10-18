#ifndef PosteriorEntropy_Demo_h
#define PosteriorEntropy_Demo_h

#include <ostream>
#include "Options.h"
#include "InfoNest/cpp/RNG.h"
#include "InfoNest/cpp/Utils.h"

namespace PosteriorEntropy
{

// A simple demo planet model specification.
class Demo
{
    private:

        // Number of data points
        static constexpr size_t N = 101;

        // Time range of data points
        static constexpr double t_min = 0.0;
        static constexpr double t_max = 1.0;
        static constexpr double t_range = t_max - t_min;
        static constexpr double dt = t_range/(N - 1);

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
        double perturb(InfoNest::RNG& rng);
        std::vector<double> simulate_data(InfoNest::RNG& rng) const;
        double log_likelihood(const std::vector<double>& ys) const;
        friend std::ostream& operator << (std::ostream& out,
                                          const Demo& demo);
        friend double distance(const Demo& x, const Demo& y);
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

double Demo::perturb(InfoNest::RNG& rng)
{
    double logH = 0.0;
    int which = rng.rand_int(4);

    if(which == 0)
    {
        tc += t_range*rng.randh();
        InfoNest::wrap(tc, t_min, t_max);
    }
    else if(which == 1)
    {
        depth = 1.0 - exp(-depth/0.1);
        depth += rng.randh();
        InfoNest::wrap(depth, 0.0, 1.0);
        depth = -0.1*log(1.0 - depth);
    }
    else if(which == 2)
    {
        width = log(width);
        logH -= -0.5*pow(width - log(0.01), 2);
        width += rng.randh();
        logH += -0.5*pow(width - log(0.01), 2);
        width = exp(width);
    }
    else if(which == 3)
    {
        sigma = log(sigma);
        logH -= -0.5*pow(sigma - log(0.01), 2);
        sigma += rng.randh();
        logH += -0.5*pow(sigma - log(0.01), 2);
        sigma = exp(sigma);
    }

    // Compute noise-free curve
    if(which != 3)
        compute_mus();

    return logH;
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

std::ostream& operator << (std::ostream& out, const std::vector<double>& ys)
{
    for(size_t i=0; i<ys.size(); ++i)
        out << ys[i] << ' ';
    return out;
}

double distance(const Demo& x, const Demo& y)
{
    return std::abs(log(x.width) - log(y.width));
}


} // namespace


#endif

