#include <iostream>
#include "Demo.h"

using namespace InfoNest;
using namespace TransitInfo;

int main()
{
    RNG rng;
    Demo params;
    params.from_prior(rng);
    auto data = params.simulate_data(rng);
    std::cout << "Parameters: " << params << std::endl;
    std::cout << "Data: " << data << std::endl;

    return 0;
}

