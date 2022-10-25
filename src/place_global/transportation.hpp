#pragma once

#include <vector>

namespace coloquinte {
/**
 * @brief Solve a generic transportation problem
 *
 * @param capacities Capacities for the regions
 * @param demands Demands for the cells
 * @param costs Cost of allocating a region to a source
 * @return The quantity allocated between each region and source
 */
std::vector<std::vector<long long> > solveTransportation(
    const std::vector<long long>& capacities,
    const std::vector<long long>& demands,
    const std::vector<std::vector<float> >& costs);

typedef std::pair<int, long long> t1D_elt;

/**
 * @brief Solve a special case of the transportation problem, where the sources
 * and sinks are all on a line
 *
 * @param sources by position and capacity
 * @param sinks by position and capacity
 * @return An "absolute position" of the sources in the sinks
 */
std::vector<long long> solveTransportation1D(std::vector<t1D_elt> sources,
                                             std::vector<t1D_elt> sinks);

}  // namespace coloquinte