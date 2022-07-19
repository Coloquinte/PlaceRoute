#pragma once

#include "coloquinte.hpp"

/**
 * Algorithms to obtain a legal placement for standard cells
 *
 */
class Legalizer {
 public:
  /**
   * Initialize the datastructure from a circuit
   */
  static Legalizer fromIspdCircuit(const Circuit &circuit);

  /**
   * Export the placement obtained to the circuit datastructure
   */
  void exportPlacement(Circuit &circuit);

  /**
   * Initialize the datastructure
   *      @param rows: Available rows for placement; must all be the right
   * height for the cells
   *      @param width: Width of the cells when placed in a row
   *      @param targetX: Target x coordinate for legalization
   *      @param targetY: Target y coordinate for legalization
   */
  Legalizer(const std::vector<Rectangle> &rows, const std::vector<int> &width,
            const std::vector<int> &targetX, const std::vector<int> &targetY);

  /**
   * Return the number of rows
   */
  int nbRows() const { return rows_.size(); }

  /**
   * Return the number of cells
   */
  int nbCells() const { return cellWidth_.size(); }

  /**
   * Return the width of the cells
   */
  const std::vector<int> &cellWidth() const { return cellWidth_; }

  /**
   * Return the target x coordinates for legalization
   */
  const std::vector<int> &cellTargetX() const { return cellTargetX_; }

  /**
   * Return the target y coordinates for legalization
   */
  const std::vector<int> &cellTargetY() const { return cellTargetY_; }

  /**
   * Return the cost model used by the algorithm
   */
  LegalizationModel costModel() const { return costModel_; }

  /**
   * Set the cost model used by the algorithm
   */
  void setCostModel(LegalizationModel m) { costModel_ = m; }

  /**
   * Run the algorithm
   */
  void run();

  /**
   * Compute the x coordinates after legalization
   */
  std::vector<int> cellLegalX() const;

  /**
   * Compute the y coordinates after legalization
   */
  std::vector<int> cellLegalY() const;

  /**
   * Check consistency of the datastructure
   */
  void check() const;

 private:
  /**
   * Place a single cell optimally
   * Return true if successful
   */
  bool placeCellOptimally(int cell);

  /**
   * Simulate placing a single cell in a given row
   * Return a pair: true if successful and the X coordinate
   */
  std::pair<bool, int> placeCellOptimally(int cell, int row) const;

  /**
   * Compute the ordering of the cells
   */
  std::vector<int> computeCellOrder(float weightX, float weightWidth,
                                    float weightY) const;

  /**
   * Returns true if the cell is to be ignored by legalization
   */
  bool isIgnored(int cell) const { return cellWidth_[cell] == -1; }

  /**
   * Materialize the placement of a cell
   */
  void doPlacement(int cell, int row, int x);

  /**
   * Undo the placement of a cell
   */
  void undoPlacement(int cell);

 private:
  // Placement data
  LegalizationModel costModel_;
  std::vector<Rectangle> rows_;
  std::vector<int> cellWidth_;
  std::vector<int> cellTargetX_;
  std::vector<int> cellTargetY_;

  // Placement status
  std::vector<std::vector<int> > rowToCells_;
  std::vector<std::vector<int> > rowToX_;
  std::vector<int> cellToRow_;
  std::vector<int> cellToX_;
  std::vector<int> cellToY_;
};
