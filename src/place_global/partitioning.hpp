#pragma once

#include "coloquinte.hpp"
#include "place_global/density_grid.hpp"

class Partitioner {
 public:
  Partitioner fromIspdCircuit(const Circuit &circuit, float sizeFactor=10.0);

 private:
  /**
   * @brief Optimize the partitioning between two bins, vertically
   */
  void partitioningX(int x, int y);
  
  /**
   * @brief Optimize the partitioning between two bins, horizontally
   */
  void partitioningY(int x, int y);

  /**
   * @brief Optimize the partitioning between four adjacent bins
   */
  void quadripartitioning(int x, int y);

  /**
   * @brief Refine the placement grid
   */
  void refine();

  /**
   * @brief Improve the current placement with one round of partitioning
   */
  void improve();

 private:
  HierarchicalDensityPlacement placement_;

  // Compressed sparse representation
  std::vector<int> netCellLimits_;
  std::vector<int> netCells_;
  std::vector<int> netPinLimits_;
  std::vector<int> netPinX_;
  std::vector<int> netPinY_;
};