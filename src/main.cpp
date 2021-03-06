#include <iostream>
#include "INModel.h"
#include "TransitDemo.h"
#include "InfoNest/cpp/Execute.hpp"

using namespace InfoNest;
using namespace PosteriorEntropy;

int main()
{
    // The InfoNest model
    using TheModel = INModel<TransitDemo, std::vector<double>>;

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
    InfoNest::execute<TheModel>(rng0, rng1, depth, num_reps, num_particles,
                                mcmc_steps, TheModel::parameter_distance,
                                Mode::conditional_entropy, mcmc_steps);

    return 0;
}

