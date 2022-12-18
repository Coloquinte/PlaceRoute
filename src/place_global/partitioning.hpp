#pragma once

#include "coloquinte.hpp"
#include "place_global/density_grid.hpp"

namespace coloquinte {
class Partitioner {
 public:
  Partitioner fromIspdCircuit(const Circuit &circuit, float sizeFactor = 10.0);

 private:
  Partitioner(const Circuit &circuit, HierarchicalDensityPlacement leg)
      : circuit_(circuit), placement_(leg) {}

  /**
   * @brief Reoptimize the partitioning across multiple bins
   */
  void reoptimize(const std::vector<std::pair<int, int> > &bins);

  /**
   * @brief Refine the placement grid
   */
  void refine();

  /**
   * @brief Improve the current placement with one round of partitioning
   */
  void improve();

 private:
  const Circuit &circuit_;
  HierarchicalDensityPlacement placement_;
};
}  // namespace coloquinte