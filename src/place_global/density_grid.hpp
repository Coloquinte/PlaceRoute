#pragma once

#include <cassert>

#include "coloquinte.hpp"

/**
 * @brief Representation of the density available for placement on a grid
 */
class DensityGrid {
 public:
  /**
   * @brief Represent a rectangular group of bins in the grid; max coords are
   * non-inclusive
   */
  struct BinGroup {
    int minXCoord;
    int maxXCoord;
    int minYCoord;
    int maxYCoord;
  };

  /**
   * @brief Initialize a grid with a simple placement area
   */
  DensityGrid(Rectangle placementArea);

  /**
   * @brief Initialize a grid with placement regions (usually rows) and
   * obstacles
   */
  DensityGrid(std::vector<Rectangle> regions,
              std::vector<Rectangle> obstacles = std::vector<Rectangle>());

  /**
   * @brief Initialize a grid from a circuit, and initialize the bins from the
   * standard cell height
   */
  static DensityGrid fromIspdCircuit(const Circuit &circuit,
                                     float sizeFactor = 10.0);

  /**
   * @brief Get the total number of bins
   */
  int nbBins() const { return nbBinsX() * nbBinsY(); }

  /**
   * @brief Get the number of bins in x direction
   */
  int nbBinsX() const { return binX_.size(); }

  /**
   * @brief Get the number of bins in y direction
   */
  int nbBinsY() const { return binY_.size(); }

  /**
   * @brief Get the sum of the capacities of the bins
   */
  long long totalCapacity() const;

  /**
   * @brief Get the capacity of a given bin
   */
  long long binCapacity(int x, int y) const { return binCapacity_[x][y]; }
  long long &binCapacity(int x, int y) { return binCapacity_[x][y]; }

  /**
   * @brief Get the center x coordinate of a given bin
   */
  float binX(int x) const {
    assert(x < nbBinsX());
    return binX_[x];
  }

  /**
   * @brief Get the center y coordinate of a given bin
   */
  float binY(int y) const {
    assert(y < nbBinsY());
    return binY_[y];
  }

  /**
   * @brief Get the x coordinate of the limit between bins at positions x and
   * x+1
   */
  float binLimitX(int x) const {
    assert(x <= nbBinsX());
    return binLimitX_[x];
  }

  /**
   * @brief Get the y coordinate of the limit between bins at positions y and
   * y+1
   */
  float binLimitY(int y) const {
    assert(y <= nbBinsY());
    return binLimitY_[y];
  }

  /**
   * @brief Compute the x coordinate of the center of mass of a group of bins
   */
  float groupCenterX(BinGroup g) const;

  /**
   * @brief Compute the y coordinate of the center of mass of a group of bins
   */
  float groupCenterY(BinGroup g) const;

  /**
   * @brief Update the density grid to this exact number of bins
   */
  void updateBinsToNumber(int binsX, int binsY);

  /**
   * @brief Update the density grid so the bins are of the given dimension or
   * smaller
   */
  void updateBinsToSize(int maxBinSize) {
    updateBinsToSize(maxBinSize, maxBinSize);
  }

  /**
   * @brief Update the density grid so the bins are of the given dimensions or
   * smaller
   */
  void updateBinsToSize(int maxXSize, int maxYSize);

  /**
   * @brief Check the consistency of the datastructure
   */
  void check() const;

 private:
  /**
   * @brief Compute the boundaries of the placement area
   */
  Rectangle computePlacementArea() const;

  /**
   * @brief Obtain cleaned-up placement regions by removing any overlap between
   * regions and obstacles
   */
  // TODO
  std::vector<Rectangle> computeActualRegions() const;

  /**
   * @brief Compute the bin capacity from the cleaned-up placement regions
   */
  // TODO
  std::vector<std::vector<long long> > computeBinCapacity() const;

 private:
  /**
   * @brief Boundaries of the placement area
   */
  Rectangle placementArea_;

  /**
   * @brief Regions where placement is possible
   */
  std::vector<Rectangle> regions_;

  /**
   * @brief Obstacles, where nothing can be placed
   */
  std::vector<Rectangle> obstacles_;

  /**
   * @brief Cleaned-up placement regions, with obstacles removed and no overlap
   */
  std::vector<Rectangle> actualRegions_;

  /**
   * Center x coordinates of the bins
   */
  std::vector<float> binX_;

  /**
   * @brief Center y coordinates of the bins
   */
  std::vector<float> binY_;

  /**
   * @brief Boundary x coordinates of the bins
   */
  std::vector<int> binLimitX_;

  /**
   * @brief Boundary y coordinates of the bins
   */
  std::vector<int> binLimitY_;

  /**
   * @brief Placement capacity of the bins (i.e. area available for placement)
   */
  std::vector<std::vector<long long> > binCapacity_;
};