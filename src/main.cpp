#include <iostream>
#include "Planet.h"
#include "Demo.h"
#include "InfoNest/cpp/Execute.hpp"

using namespace InfoNest;
using namespace TransitInfo;

int main()
{
    // The InfoNest model
    using INModel = Planet<Demo, std::vector<double>>;

    // Create random number generators
    // The first one is used to generate reference points
    // and the second is used for the Nested Sampling
    unsigned long seed0 = 0;
    unsigned long seed1 = time(0);
    InfoNest::RNG rng0(seed0);
    InfoNest::RNG rng1(seed1);

    // Define run parameters
    constexpr double depth         = 30.0;
    constexpr size_t num_reps      = 1000;
    constexpr size_t num_particles = 10;
    constexpr size_t mcmc_steps    = 1000;

    // Do the run.
    InfoNest::execute<INModel>(rng0, rng1, depth, num_reps, num_particles,
                               mcmc_steps, INModel::parameter_distance,
                               Mode::conditional_entropy, mcmc_steps);

    return 0;
}

