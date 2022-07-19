#pragma once

#include "coloquinte.hpp"

#include <cassert>

/**
 * Representation of the density available for placement on a grid
 */
class DensityGrid {
  public:
    /**
     * Represent a rectangular group of bins in the grid; max coords are non-inclusive
     */
    struct BinGroup {
        int minXCoord;
        int maxXCoord;
        int minYCoord;
        int maxYCoord;
    };

    /**
     * Initialize a grid with a simple placement area
     */
    DensityGrid(Rectangle placementArea);

    /**
     * Initialize a grid with placement regions (usually rows) and obstacles
     */
    DensityGrid(std::vector<Rectangle> regions, std::vector<Rectangle> obstacles=std::vector<Rectangle>());

    /**
     * Initialize a grid from a circuit, and initialize the bins from the standard cell height
     */
    static DensityGrid fromIspdCircuit(const Circuit &circuit, float sizeFactor=10.0);

    /**
     * Get the total number of bins
     */
    int nbBins() const { return nbBinsX() * nbBinsY(); }

    /**
     * Get the number of bins in x direction
     */
    int nbBinsX() const { return binX_.size(); }

    /**
     * Get the number of bins in y direction
     */
    int nbBinsY() const { return binY_.size(); }

    /**
     * Get the total capacity of the grid
     */
    long long totalCapacity() const;

    /**
     * Get the capacity of a given bin
     */
    long long binCapacity(int x, int y) const { return binCapacity_[x][y]; }
    long long &binCapacity(int x, int y) { return binCapacity_[x][y]; }

    /**
     * Get the center x coordinate of a given bin
     */
    float binX(int x) const { assert (x < nbBinsX()); return binX_[x]; }

    /**
     * Get the center y coordinate of a given bin
     */
    float binY(int y) const { assert (y < nbBinsY()); return binY_[y]; }

    /**
     * Get the x coordinate of the limit between bins at positions x and x+1
     */
    float binLimitX(int x) const { assert (x <= nbBinsX()); return binLimitX_[x]; }

    /**
     * Get the y coordinate of the limit between bins at positions y and y+1
     */
    float binLimitY(int y) const { assert (y <= nbBinsY()); return binLimitY_[y]; }

    /**
     * Compute the x coordinate of the center of mass of a group of bins
     */
    float groupCenterX(BinGroup g) const;

    /**
     * Compute the y coordinate of the center of mass of a group of bins
     */
    float groupCenterY(BinGroup g) const;

    /**
     * Update the density grid to this exact number of bins
     */
    void updateBinsToNumber(int binsX, int binsY);

    /**
     * Update the density grid so the bins are of the given dimension or smaller
     */
    void updateBinsToSize(int maxBinSize) { updateBinsToSize(maxBinSize, maxBinSize); }

    /**
     * Update the density grid so the bins are of the given dimensions or smaller
     */
    void updateBinsToSize(int maxXSize, int maxYSize);

    /**
     * Check the consistency of the datastructure
     */
    void check() const;

  private:
    /**
     * Compute the boundaries of the placement area
     */
    Rectangle computePlacementArea() const;

    /**
     * Obtain cleaned-up placement regions by removing any overlap between regions and obstacles
     */
    // TODO
    std::vector<Rectangle> computeActualRegions() const;

    /**
     * Compute the bin capacity from the cleaned-up placement regions
     */
    // TODO
    std::vector<std::vector<long long> > computeBinCapacity() const;

  private:
    /**
     * Boundaries of the placement area
     */
    Rectangle placementArea_;

    /**
     * Regions where placement is possible
     */
    std::vector<Rectangle> regions_;

    /**
     * Obstacles, where nothing can be placed
     */
    std::vector<Rectangle> obstacles_;

    /**
     * Cleaned-up placement regions, with obstacles removed and no overlap
     */
    std::vector<Rectangle> actualRegions_;

    /**
     * Center x coordinates of the bins
     */
    std::vector<float> binX_;

    /**
     * Center y coordinates of the bins
     */
    std::vector<float> binY_;

    /**
     * Boundary x coordinates of the bins
     */
    std::vector<int> binLimitX_;

    /**
     * Boundary y coordinates of the bins
     */
    std::vector<int> binLimitY_;

    /**
     * Placement capacity of the bins (i.e. area available for placement)
     */
    std::vector<std::vector<long long> > binCapacity_;
};
