#pragma once

#include <vector>

std::vector<std::vector<long long> > solveTransportation(
    std::vector<long long> const& capacities,
    std::vector<long long> const& demands,
    std::vector<std::vector<float> > const& costs);