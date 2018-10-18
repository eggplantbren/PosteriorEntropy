#ifndef PosteriorEntropy_Sinewave_h
#define PosteriorEntropy_Sinewave_h

#include <ostream>
#include "InfoNest/cpp/RNG.h"
#include "InfoNest/cpp/Utils.h"

namespace PosteriorEntropy
{

// A sinusoidal situation
class Sinewave
{
    private:

        // Number of data points
        static constexpr size_t N = 21;

        // Time range of data points
        static constexpr double t_min = 0.0;
        static constexpr double t_max = 1.0;
        static constexpr double t_range = t_max - t_min;
        static constexpr double dt = t_range/(N - 1);

        // Width of observations
        static constexpr double h = 0.01;

        // Noise SD
        static constexpr double sigma = 1.0;

    private:

        // Just parameterise by amplitude, period, phase
        double amplitude, period, phase;

        // Noise-free model curve
        std::vector<double> mus;
        void compute_mus();

    public:

        Sinewave();
        void from_prior(InfoNest::RNG& rng);
        double perturb(InfoNest::RNG& rng);
        std::vector<double> simulate_data(InfoNest::RNG& rng) const;
        double log_likelihood(const std::vector<double>& ys) const;
        friend std::ostream& operator << (std::ostream& out,
                                          const Sinewave& demo);
        friend double distance(const Sinewave& x, const Sinewave& y);
};


/****************** IMPLEMENTATION FOLLOWS ********************/

Sinewave::Sinewave()
:mus(N)
{

}

void Sinewave::from_prior(InfoNest::RNG& rng)
{
    // Naive priors
    amplitude = 100.0*rng.rand();
    period = 0.05 + 0.95*rng.rand();
    phase = 2.0*M_PI*rng.rand();

    // Compute noise-free curve
    compute_mus();
}

double Sinewave::perturb(InfoNest::RNG& rng)
{
    double logH = 0.0;
    int which = rng.rand_int(4);

    if(which == 0)
    {
        amplitude += 100.0*rng.randh();
        InfoNest::wrap(amplitude, 0.0, 100.0);
    }
    else if(which == 1)
    {
        period += 0.95*rng.randh();
        InfoNest::wrap(period, 0.05, 1.0);
    }
    else
    {
        phase += 2.0*M_PI*rng.randh();
        InfoNest::wrap(phase, 0.0, 2.0*M_PI);
    }

    // Compute noise-free curve
    compute_mus();

    return logH;
}

void Sinewave::compute_mus()
{
    double t, a, b;
    InfoNest::RNG temp_rng(123);
    double hh = h;

    for(size_t i=0; i<N; ++i)
    {
        // Timestamp of the data point
        t = t_min + i*dt;

        // Left and right edge of observations
        hh = 2*h*temp_rng.rand();
        a = t - 0.5*hh;
        b = t + 0.5*hh;

        mus[i] = amplitude*(-period*cos(phase + 2*M_PI*b/period)/(2*M_PI)
                 + period*cos(2*M_PI*a/period + phase)/(2*M_PI)) / h;
    }
}


std::vector<double> Sinewave::simulate_data(InfoNest::RNG& rng) const
{
    std::vector<double> ys(N);

    for(size_t i=0; i<N; ++i)
    {
        // Add the noise
        ys[i] = mus[i] + sigma*rng.randn();
    }

    return ys;
}


double Sinewave::log_likelihood(const std::vector<double>& ys) const
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
                                          const Sinewave& s)
{
    out << s.amplitude << ' ' << s.period << ' ' << s.phase << ' ';
    return out;
}

std::ostream& operator << (std::ostream& out, const std::vector<double>& ys)
{
    for(size_t i=0; i<ys.size(); ++i)
        out << ys[i] << ' ';
    return out;
}

double distance(const Sinewave& x, const Sinewave& y)
{
    return std::abs(x.period - y.period);
}

} // namespace


#endif

