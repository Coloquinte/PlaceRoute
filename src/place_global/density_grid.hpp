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
  DensityGrid(int binSize, Rectangle placementArea);

  /**
   * @brief Initialize a grid with placement regions (usually rows) and
   * obstacles
   */
  DensityGrid(int binSize, std::vector<Rectangle> regions,
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
   * @brief Compute the total capacity of a group of bins
   */
  long long binCapacity(BinGroup g) const;

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
   * @brief Check the consistency of the datastructure
   */
  void check() const;

 private:
  /**
   * @brief Initialize a grid from limits and capacity
   */
  DensityGrid(std::vector<int> xLimits, std::vector<int> yLimits,
              std::vector<std::vector<long long> > binCapacity);

  /**
   * @brief Compute the centers of the bins
   */
  void updateBinCenters();

  /**
   * @brief Compute the bin capacity from their sizes
   */
  void updateBinCapacity();

  /**
   * @brief Compute the bin capacity from the cleaned-up placement regions
   */
  void updateBinCapacity(const std::vector<Rectangle> &regions);

  /**
   * @brief Update the density grid to this exact number of bins
   */
  void updateBinsToNumber(int binsX, int binsY);

  /**
   * @brief Update the density grid so the bins are of the given dimension or
   * smaller
   */
  void updateBinsToSize(int maxSize);

  /**
   * @brief Compute the boundaries of the placement area
   */
  static Rectangle computePlacementArea(const std::vector<Rectangle> &regions);

  /**
   * @brief Obtain cleaned-up placement regions by removing any overlap between
   * regions and obstacles
   */
  static std::vector<Rectangle> computeActualRegions(
      const std::vector<Rectangle> &regions,
      const std::vector<Rectangle> &obstacles);

 private:
  /**
   * @brief Boundaries of the placement area
   */
  Rectangle placementArea_;

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

  friend class HierarchicalDensityPlacement;
};

/**
 * @brief Represent the state of a hierarchical legalizer or partitioner, with a
 * density grid that can be refined and coarsened at will.
 *
 */
class HierarchicalDensityPlacement {
 public:
  /**
   * @brief Initialize the datastructure from the grid and the cell demands
   * (single bin)
   */
  HierarchicalDensityPlacement(DensityGrid grid, std::vector<int> cellDemand);

  /**
   * @brief Initialize from a circuit
   */
  static HierarchicalDensityPlacement fromIspdCircuit(const Circuit &circuit,
                                                      float sizeFactor = 10.0);

  /**
   * @brief Access the underlying grid
   */
  const DensityGrid &grid() const { return grid_; }
  DensityGrid &grid() { return grid_; }

  /**
   * @brief Get the total number of bins
   */
  int nbBins() const { return nbBinsX() * nbBinsY(); }

  /**
   * @brief Get the number of bins in the x direction with the current view
   */
  int nbBinsX() const { return nbBinsX(levelX_); };

  /**
   * @brief Get the number of bins in the y direction with the current view
   */
  int nbBinsY() const { return nbBinsY(levelY_); };

  /**
   * @brief Get the number of cells
   */
  int nbCells() const { return cellDemand_.size(); }

  /**
   * @brief Get the current coarsening level in the x direction. 0 is the fully
   * refined grid.
   */
  int levelX() const { return levelX_; }

  /**
   * @brief Get the current coarsening level in the y direction. 0 is the fully
   * refined grid.
   */
  int levelY() const { return levelY_; }

  /**
   * @brief Get the total number of coarsening levels in the x direction
   */
  int nbLevelX() const { return xLimits_.size(); }

  /**
   * @brief Get the total number of coarsening levels in the y direction
   */
  int nbLevelY() const { return yLimits_.size(); }

  /**
   * @brief Get the demand for a given cell
   */
  int cellDemand(int c) const {
    assert(c < nbCells());
    return cellDemand_[c];
  }

  /**
   * @brief Get the x coordinate of the boundary between bins with the current
   * view
   */
  int binLimitX(int x) const {
    assert(x <= nbBinsX());
    return grid_.binLimitX(xLimits_[levelX_][x]);
  }

  /**
   * @brief Get the y coordinate of the boundary between bins with the current
   * view
   */
  int binLimitY(int y) const {
    assert(y <= nbBinsY());
    return grid_.binLimitY(yLimits_[levelY_][y]);
  }

  /**
   * @brief Index of the parent bin at the upper level of coarsening (x
   * direction). This is useful to decide which bins to reoptimize after a
   * refinement step.
   */
  int parentX(int x) const { return parentX(levelX_, x); }

  /**
   * @brief Index of the parent bin at the upper level of coarsening (y
   * direction). This is useful to decide which bins to reoptimize after a
   * refinement step.
   */
  int parentY(int y) const { return parentY(levelY_, y); }

  /**
   * @brief Get the bin index corresponding to this x coordinate
   */
  int findBinByX(int x) const;

  /**
   * @brief Get the bin index corresponding to this y coordinate
   */
  int findBinByY(int y) const;

  /**
   * @brief Get the capacity of a given bin in the current view
   */
  long long binCapacity(int x, int y) const {
    return grid_.binCapacity(getGroup(x, y));
  }
  /**
   * @brief Get the x center of a given bin in the current view
   */
  float binX(int x, int y) const { return grid_.groupCenterX(getGroup(x, y)); }

  /**
   * @brief Get the x center of a given bin in the current view
   */
  float binY(int x, int y) const { return grid_.groupCenterY(getGroup(x, y)); }

  /**
   * @brief Get the usage of a given bin in the current view
   */
  long long binUsage(int x, int y) const;

  /**
   * @brief Return the x coordinates for the cells (center of the bin)
   */
  std::vector<float> simpleCoordX() const;

  /**
   * @brief Return the y coordinates for the cells (center of the bin)
   */
  std::vector<float> simpleCoordY() const;

  /**
   * @brief Return the x coordinates for the cells (spread in the bin according
   * to the target coordinates)
   */
  std::vector<float> spreadCoordX(const std::vector<float> &target) const;

  /**
   * @brief Return the y coordinates for the cells (spread in the bin according
   * to the target coordinates)
   */
  std::vector<float> spreadCoordY(const std::vector<float> &target) const;

  /**
   * @brief Get the sum of the demands of the cells
   */
  long long totalDemand() const;

  /**
   * @brief Compute the total overflowing area of this placement
   */
  long long totalOverflow() const;

  /**
   * @brief Compute the overflowing area of this placement, relative to the
   * total demand
   */
  float overflowRatio() const;

  /**
   * @brief Return the cells allocated to a given bin in the current view
   */
  const std::vector<int> &binCells(int x, int y) const {
    return binCells_[x][y];
  }

  /**
   * @brief Update the cells in the given bin
   */
  void setBinCells(int x, int y, std::vector<int> cells);

  /**
   * @brief Return the x index of the bin where the cell is located
   */
  int cellBinX(int c) const { return cellBinX_[c]; }

  /**
   * @brief Return the y index of the bin where the cell is located
   */
  int cellBinY(int c) const { return cellBinY_[c]; }

  /**
   * @brief Refine vertically (more bins in the x direction). The cells are
   * assigned to one side without rebalancing.
   */
  void refineX();

  /**
   * @brief Refine horizontally (more bins in the y direction). The cells are
   * assigned to one side without rebalancing.
   */
  void refineY();

  /**
   * @brief Coarsen vertically (less bins in the x direction).
   */
  void coarsenX();

  /**
   * @brief Coarsen horizontally (less bins in the y direction).
   */
  void coarsenY();

  /**
   * @brief Check the consistency of the datastructure
   */
  void check() const;

 private:
  DensityGrid::BinGroup getGroup(int x, int y) const {
    return DensityGrid::BinGroup(
        {xLimits_[levelX_][x], xLimits_[levelX_][x + 1], yLimits_[levelY_][y],
         yLimits_[levelY_][y + 1]});
  }

  /**
   * @brief Construct the datastructures describing the hierarchy
   */
  void setupHierarchy();

  /**
   * @brief Update the cell locations
   */
  void updateCellToBin();

  /**
   * @brief Get the number of bins in the x direction at the given level
   */
  int nbBinsX(int lvl) const { return xLimits_[lvl].size() - 1; };

  /**
   * @brief Get the number of bins in the y direction at the given level
   */
  int nbBinsY(int lvl) const { return yLimits_[lvl].size() - 1; };

  /**
   * @brief Index of the parent bin at the upper level of coarsening (x
   * direction).
   */
  int parentX(int lvl, int x) const { return parentX_[lvl][x]; }

  /**
   * @brief Index of the parent bin at the upper level of coarsening (y
   * direction).
   */
  int parentY(int lvl, int y) const { return parentY_[lvl][y]; }

 private:
  /**
   * @brief Underlying density grid
   */
  DensityGrid grid_;

  /**
   * @brief Level of coarsening in the x direction
   */
  int levelX_;

  /**
   * @brief Level of coarsening in the y direction
   */
  int levelY_;

  /**
   * @brief Index of the bin boundaries at each coarsening level in the x
   * direction
   */
  std::vector<std::vector<int> > xLimits_;

  /**
   * @brief Index of the bin boundaries at each coarsening level in the y
   * direction
   */
  std::vector<std::vector<int> > yLimits_;

  /**
   * @brief Index of the parent in the x direction
   */
  std::vector<std::vector<int> > parentX_;

  /**
   * @brief Index of the parent in the y direction
   */
  std::vector<std::vector<int> > parentY_;

  /**
   * @brief Demand of the cells
   */
  std::vector<int> cellDemand_;

  /**
   * @brief Allocation of the cells to the bins
   */
  std::vector<std::vector<std::vector<int> > > binCells_;

  std::vector<int> cellBinX_;
  std::vector<int> cellBinY_;
};
