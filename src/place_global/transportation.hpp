#pragma once

#include <vector>

std::vector<std::vector<long long> > solveTransportation(
    const std::vector<long long>& capacities,
    const std::vector<long long>& demands,
    const std::vector<std::vector<float> >& costs);