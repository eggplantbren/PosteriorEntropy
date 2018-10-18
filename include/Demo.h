#ifndef TransitInfo_Demo_h
#define TransitInfo_Demo_h

#include "InfoNest/cpp/RNG.h"

namespace TransitInfo
{

// A simple demo planet model specification.
class Demo
{
    private:

        // Just parameterise depth and width of the transit.
        double depth;
        double width;

    public:

        void from_prior(InfoNest::RNG& rng);
};



} // namespace


#endif

