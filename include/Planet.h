#ifndef TransitInfo_Planet_h
#define TransitInfo_Planet_h

#include "InfoNest/cpp/RNG.h"

namespace TransitInfo
{

// T is planet parameterisation
// D is data type (e.g. std::vector<double>)
// The following must be implemented:
//      void T::from_prior(InfoNest::RNG&)
//      D    T::simulate_data(InfoNest::RNG&) const
template <typename T, typename D>
class Planet
{
    private:

        // A point in the joint space.
        T params;
        D data;

    public:

        Planet();

        // Generate from the distribution
        void generate(InfoNest::RNG& rng);

        // Metropolis proposal
        double perturb(InfoNest::RNG& rng);

        // Printing to stream
        void print(std::ostream& out) const;

        // Parameter distance function
        static double parameter_distance
                            (const Planet& x, const Planet& y);

};




/****************** IMPLEMENTATION FOLLOWS ********************/

template <typename T, typename D>
Planet<T, D>::Planet()
{

}


template <typename T, typename D>
void Planet<T, D>::generate(InfoNest::RNG& rng)
{
    params.from_prior(rng);
    data = params.simulate_data(rng);
}


} // namespace


#endif

