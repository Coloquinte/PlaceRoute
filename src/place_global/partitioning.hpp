#pragma once

#include "coloquinte.hpp"
#include "place_global/density_grid.hpp"

namespace coloquinte {
class Partitioner : public HierarchicalDensityPlacement {
 public:
  static Partitioner fromIspdCircuit(const Circuit &circuit, float sizeFactor = 10.0);

  /**
   * @brief Run the whole partitioning algorithm
   */
  void run();

  /**
   * @brief Refine the placement grid
   */
  void refine();

  /**
   * @brief Improve the current placement with one round of partitioning
   */
  void improve();

 private:
  Partitioner(const Circuit &circuit, HierarchicalDensityPlacement leg)
      : circuit_(circuit), HierarchicalDensityPlacement(leg) {}

  /**
   * @brief Reoptimize the partitioning across multiple bins
   */
  void reoptimize(const std::vector<std::pair<int, int> > &bins);

  /**
   * @brief Improve neighbouring bin pairs in the x direction
   */
  void improveXNeighbours(bool sameParent = true);

  /**
   * @brief Improve neighbouring bin pairs in the y direction
   */
  void improveYNeighbours(bool sameParent = true);

 private:
  const Circuit &circuit_;
};
}  // namespace coloquinte