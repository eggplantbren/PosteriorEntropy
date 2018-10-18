#ifndef PosteriorEntropy_MassDemo_h
#define PosteriorEntropy_MassDemo_h

#include <ostream>
#include "InfoNest/cpp/RNG.h"
#include "InfoNest/cpp/Utils.h"

namespace PosteriorEntropy
{

// Total mass of a source measured with noise. Is it better to be model-based,
// even if the model is wrong?
class MassDemo
{
    private:

        // Number of data points
        static constexpr size_t N = 1001;

        // Time range of data points
        static constexpr double x_min = 0.0;
        static constexpr double x_max = 1.0;
        static constexpr double x_range = x_max - x_min;
        static constexpr double dx = x_range/(N - 1);

        // Noise level
        static constexpr double sigma = 1.0;

    private:

        // Central position, amplitude, width
        double xc, A, width;
        
        // Noise-free model curve
        std::vector<double> mus;
        void compute_mus();

    public:

        MassDemo();
        void from_prior(InfoNest::RNG& rng);
        double perturb(InfoNest::RNG& rng);
        std::vector<double> simulate_data(InfoNest::RNG& rng) const;
        double log_likelihood(const std::vector<double>& ys) const;
        friend std::ostream& operator << (std::ostream& out,
                                          const MassDemo& demo);
        friend double distance(const MassDemo& x, const MassDemo& y);
};


/****************** IMPLEMENTATION FOLLOWS ********************/

MassDemo::MassDemo()
:mus(N)
{

}

void MassDemo::from_prior(InfoNest::RNG& rng)
{
    xc = x_min + 0.25*x_range + 0.5*x_range*rng.rand();
    A = exp(3.0*rng.randn());
    width = 0.1*x_range*rng.rand();

    // Compute noise-free curve
    compute_mus();
}

double MassDemo::perturb(InfoNest::RNG& rng)
{
    double logH = 0.0;
    int which = rng.rand_int(3);

    if(which == 0)
    {
        xc += 0.5*x_range*rng.randh();
        InfoNest::wrap(xc, x_min + 0.25*x_range, x_min + 0.75*x_range);
    }
    else if(which == 1)
    {
        A = log(A);
        logH -= -0.5*pow(A/3.0, 2);
        A += 3.0*rng.randh();
        logH += -0.5*pow(A/3.0, 2);
        A = exp(A);
    }
    else
    {
        width += 0.1*x_range*rng.randh();
        InfoNest::wrap(width, 0.0, 0.1*x_range);
    }

    // Compute noise-free curve
    compute_mus();

    return logH;
}

void MassDemo::compute_mus()
{
    double x;
    for(size_t i=0; i<N; ++i)
    {
        x = x_min + i*dx;
        mus[i] = A*exp(-0.5*pow((x - xc)/width, 2));
    }
}


std::vector<double> MassDemo::simulate_data(InfoNest::RNG& rng) const
{
    std::vector<double> ys(N);

    for(size_t i=0; i<N; ++i)
    {
        // Add the noise
        ys[i] = mus[i] + sigma*rng.randn();
    }

    return ys;
}


double MassDemo::log_likelihood(const std::vector<double>& ys) const
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
                                          const MassDemo& demo)
{
    out << demo.xc << ' ' << demo.A << ' ';
    out << demo.width;
    return out;
}

std::ostream& operator << (std::ostream& out, const std::vector<double>& ys)
{
    for(size_t i=0; i<ys.size(); ++i)
        out << ys[i] << ' ';
    return out;
}

double distance(const MassDemo& x, const MassDemo& y)
{
    double mass1 = 0.0;
    double mass2 = 0.0;
    for(double mu: x.mus)
        mass1 += mu*MassDemo::dx;
    for(double mu: y.mus)
        mass2 += mu*MassDemo::dx;

    return std::abs(log(mass1) - log(mass2));
}


} // namespace


#endif

