#pragma once

#include "coloquinte.hpp"
#include "place_global/density_grid.hpp"

/**
 * Representation of an almost legalized placement, where density constraints
 * are met per bin
 */
class DensityLegalizer : public HierarchicalDensityPlacement {
 public:
  /**
   * @brief Initialize the datastructure from the grid and the cell demands
   * (single bin)
   */
  DensityLegalizer(DensityGrid grid, std::vector<int> cellDemand);

  /**
   * @brief Initialize from an existing hierarchical placement
   */
  DensityLegalizer(HierarchicalDensityPlacement pl);

  /**
   * @brief Initialize from a circuit
   */
  static DensityLegalizer fromIspdCircuit(const Circuit &circuit,
                                          float sizeFactor = 10.0);

  /**
   * @brief Target x position for the cell
   */
  float cellTargetX(int c) const { return cellTargetX_[c]; }

  /**
   * @brief Target y position for the cell
   */
  float cellTargetY(int c) const { return cellTargetY_[c]; }

  /**
   * @brief Get the cost model used
   */
  LegalizationModel costModel() const { return costModel_; }

  /**
   * @brief Set the cost model to use
   */
  void setCostModel(LegalizationModel model) { costModel_ = model; }

  /**
   * @brief Update the demands of the cells
   */
  void updateCellDemand(std::vector<int> cellDemand) {
    assert(cellDemand.size() == nbCells());
    cellDemand_ = cellDemand;
  }

  /**
   * @brief Update the x coordinates to target
   */
  void updateCellTargetX(std::vector<float> cellTargetX) {
    assert(cellTargetX.size() == nbCells());
    cellTargetX_ = cellTargetX;
  }
  /**
   * @brief Update the y coordinates to target
   */
  void updateCellTargetY(std::vector<float> cellTargetY) {
    assert(cellTargetY.size() == nbCells());
    cellTargetY_ = cellTargetY;
  }

  /**
   * @brief Return the mean displacement with the current cost model
   */
  float quality() const { return meanDistance(costModel_); }

  /**
   * @brief Return the mean displacement with the given cost model
   */
  float meanDistance(LegalizationModel model) const;

  /**
   * @brief Return the root-mean-square displacement with the given cost model
   */
  float rmsDistance(LegalizationModel model) const;

  /**
   * @brief Return the maximum displacement with the given cost model
   */
  float maxDistance(LegalizationModel model) const;

  /**
   * @brief Run the whole legalization process
   */
  void run();

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
  void report(bool verbose=false) const;

 private:
  /**
   * @brief Improve all neighbouring bins in the x direction
   */
  void improveX(bool sameParent);

  /**
   * @brief Improve all neighbouring bins in the y direction
   */
  void improveY(bool sameParent);

  /**
   * @brief Redo the bisection for two bins
   */
  void rebisect(int x1, int y1, int x2, int y2);

  // Bisection algorithm helpers
  std::vector<std::pair<float, int> > computeCellCosts(
      float cx1, float cy1, float cx2, float cy2, std::vector<int> cells) const;
  int findIdealSplitPos(
      const std::vector<std::pair<float, int> > &cellCosts) const;
  int findConstrainedSplitPos(
      const std::vector<std::pair<float, int> > &cellCosts, int targetPos,
      long long capa1, long long capa2) const;
  std::pair<std::vector<int>, std::vector<int> > doSplit(
      const std::vector<std::pair<float, int> > &cellCosts, int ind) const;

  /**
   * @brief Return all distances with a given cost model
   */
  std::vector<float> allDistances(LegalizationModel model) const;

 private:
  LegalizationModel costModel_;
  std::vector<float> cellTargetX_;
  std::vector<float> cellTargetY_;
};
