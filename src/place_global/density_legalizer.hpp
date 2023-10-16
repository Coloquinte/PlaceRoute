#pragma once

#include "coloquinte.hpp"
#include "place_global/density_grid.hpp"

namespace coloquinte {
/**
 * Representation of an almost legalized placement, where density constraints
 * are met per bin
 */
class DensityLegalizer : public HierarchicalDensityPlacement {
 public:
  /**
   * @brief Parameters to perform rough legalization
   */
  struct Parameters {
    // Number of improvement steps at each grid coarsening level
    int nbSteps;
    LegalizationModel costModel;
    int lineReoptSize;
    int lineReoptOverlap;
    int diagReoptSize;
    int diagReoptOverlap;
    int squareReoptSize;
    int squareReoptOverlap;
    double quadraticPenaltyFactor;
    double coarseningLimit;
    bool unidimensionalTransport;

    Parameters();
  };

  /**
   * @brief Initialize the datastructure from the grid and the cell demands
   * (single bin)
   */
  DensityLegalizer(DensityGrid grid, std::vector<int> cellDemand,
                   const Parameters &params = Parameters());

  /**
   * @brief Initialize from an existing hierarchical placement
   */
  explicit DensityLegalizer(HierarchicalDensityPlacement pl,
                            const Parameters &params = Parameters());

  /**
   * @brief Initialize from a circuit
   */
  static DensityLegalizer fromIspdCircuit(const Circuit &circuit,
                                          float sizeFactor, float sideMargin);

  /**
   * @brief Access the parameters
   */
  const Parameters &params() const { return params_; }

  /**
   * @brief Set the parameters
   */
  void setParams(const Parameters &params) { params_ = params; }

  /**
   * @brief Target x position for the cell
   */
  float cellTargetX(int c) const { return cellTargetX_[c]; }

  /**
   * @brief Target y position for the cell
   */
  float cellTargetY(int c) const { return cellTargetY_[c]; }

  /**
   * @brief Update the x coordinates to target
   */
  void updateCellTargetX(const std::vector<float> &cellTargetX) {
    assert((int) cellTargetX.size() == nbCells());
    cellTargetX_ = cellTargetX;
  }
  /**
   * @brief Update the y coordinates to target
   */
  void updateCellTargetY(const std::vector<float> &cellTargetY) {
    assert((int) cellTargetY.size() == nbCells());
    cellTargetY_ = cellTargetY;
  }

  /**
   * @brief Return the mean distance
   */
  float meanDistance() const;

  /**
   * @brief Return the root-mean-square distance
   */
  float rmsDistance() const;

  /**
   * @brief Return the maximum distance
   */
  float maxDistance() const;

  /**
   * @brief Run the whole legalization process
   */
  void run();

  /**
   * @brief Run the coarsening part of the legalization process
   */
  void runCoarsening();

  /**
   * @brief Run the refinement part of the legalization process
   */
  void runRefinement();

  /**
   * @brief Do one refinement step
   */
  void refine();

  /**
   * @brief Improve the solution at the current refinement level
   */
  void improve();

  /**
   * @brief Check the consistency of the datastructure
   */
  void check() const;

  /**
   * @brief Report on stdout
   */
  void report(bool verbose = false) const;

 private:
  /**
   * @brief Improve neighbouring bin pairs in the x direction
   */
  void improveXNeighbours(bool sameParent = true);

  /**
   * @brief Improve neighbouring bin pairs in the y direction
   */
  void improveYNeighbours(bool sameParent = true);

  /**
   * @brief Improve neighbouring bin squares
   */
  void improveSquareNeighbours(bool sameParentX = true, bool sameParentY = true);

  /**
   * @brief Improve groups of side-adjacent bins
   */
  void improveXY();

  /**
   * @brief Improve all bins in x or y direction at once
   */
  void improveUnidimensionalTransport();

  /**
   * @brief Improve all bins in x direction
   */
  void improveXTransport();

  /**
   * @brief Improve all bins in x direction
   */
  void improveYTransport();

  /**
   * @brief Improve groups of diagonally adjacent bins
   */
  void improveDiagonals();

  /**
   * @brief Improve squares of adjacent bins
   */
  void improveSquare();

  /**
   * @brief Generic improvement of rectangles applied over the grid
   */
  void improveRectangles(int width, int height, int strideX, int strideY,
                         int startX, int startY);

  /**
   * @brief Generic improvement of diagonal rectangles applied over the grid
   */
  void improveDiagonalRectangles(int xmySize, int xpySize, int strideX,
                                 int strideY, int startX, int startY);

  /**
   * @brief Redo the bisection for two bins
   */
  void rebisect(int x1, int y1, int x2, int y2);

  /**
   * @brief Improve a single rectangle of the grid
   */
  void improveRectangle(int i, int j, int width, int height);

  /**
   * @brief Redo the distribution using a transportation algorithm
   */
  void reoptimize(const std::vector<std::pair<int, int> > &bins);

  // Bisection algorithm helpers
  std::vector<std::pair<float, int> > computeCellCosts(
      float cx1, float cy1, float cx2, float cy2,
      const std::vector<int> &cells) const;
  int findIdealSplitPos(
      const std::vector<std::pair<float, int> > &cellCosts) const;
  int findConstrainedSplitPos(
      const std::vector<std::pair<float, int> > &cellCosts, int targetPos,
      long long capa1, long long capa2) const;
  std::pair<std::vector<int>, std::vector<int> > doSplit(
      const std::vector<std::pair<float, int> > &cellCosts, int ind) const;

  /**
   * @brief Return the distance given the parameters (cost model  + penalty)
   */
  float distance(float x, float y) const;

  /**
   * @brief Return all distances
   */
  std::vector<float> allDistances() const;

 private:
  Parameters params_;

  std::vector<float> cellTargetX_;
  std::vector<float> cellTargetY_;
};
}  // namespace coloquinte