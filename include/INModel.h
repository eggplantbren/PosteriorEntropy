#ifndef PosteriorEntropy_INModel_h
#define PosteriorEntropy_INModel_h

#include "InfoNest/cpp/RNG.h"

namespace PosteriorEntropy
{

// T is planet parameterisation
// D is data type (e.g. std::vector<double>)

// The following must be implemented:

//      void T::from_prior(InfoNest::RNG&)
//      D    T::simulate_data(InfoNest::RNG&) const
//      double T::log_likelihood(const D&) const
//      double T::perturb(InfoNest::RNG&)
//      std::ostream& operator << (std::ostream&, const T&)
//      std::ostream& operator << (std::ostream&, const D&)
//      double distance(const T&, const T&)
template <typename T, typename D>
class INModel
{
    private:

        // A point in the joint space.
        T params;
        D data;

        // Log likelihood
        double logL;

    public:

        INModel();

        // Generate from the distribution
        void generate(InfoNest::RNG& rng);

        // Metropolis proposal
        double perturb(InfoNest::RNG& rng);

        // Printing to stream
        void print(std::ostream& out) const;

        // Parameter distance function
        static double parameter_distance
                            (const INModel& x, const INModel& y);

};




/****************** IMPLEMENTATION FOLLOWS ********************/

template <typename T, typename D>
INModel<T, D>::INModel()
{

}


template <typename T, typename D>
void INModel<T, D>::generate(InfoNest::RNG& rng)
{
    params.from_prior(rng);
    data = params.simulate_data(rng);
    logL = params.log_likelihood(data);
}

template <typename T, typename D>
double INModel<T, D>::perturb(InfoNest::RNG& rng)
{
    double logH = -logL;

    logH += params.perturb(rng);
    logL = params.log_likelihood(data);
    logH += logL;

    return logH;
}

// Satisfies clang
std::ostream& operator << (std::ostream& out, const std::vector<double>& ys);

template <typename T, typename D>
void INModel<T, D>::print(std::ostream& out) const
{
    out << params << ' ' << data;
}

template <typename T, typename D>
double INModel<T, D>::parameter_distance
                            (const INModel<T, D>& x, const INModel<T, D>& y)
{
    return distance(x.params, y.params);
}

} // namespace


#endif

