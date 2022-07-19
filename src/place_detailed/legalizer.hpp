#pragma once

#include "coloquinte.hpp"

/**
 * @brief Algorithms to obtain a legal placement for standard cells
 *
 */
class Legalizer {
 public:
  /**
   * @brief Initialize the datastructure from a circuit
   */
  static Legalizer fromIspdCircuit(const Circuit &circuit);

  /**
   * @brief Export the placement obtained to the circuit datastructure
   */
  void exportPlacement(Circuit &circuit);

  /**
   * @brief Initialize the datastructure
   *      @param rows: Available rows for placement; must all be the right
   * height for the cells
   *      @param width: Width of the cells when placed in a row
   *      @param targetX: Target x coordinate for legalization
   *      @param targetY: Target y coordinate for legalization
   */
  Legalizer(const std::vector<Rectangle> &rows, const std::vector<int> &width,
            const std::vector<int> &targetX, const std::vector<int> &targetY);

  /**
   * @brief Return the number of rows
   */
  int nbRows() const { return rows_.size(); }

  /**
   * @brief Return the number of cells
   */
  int nbCells() const { return cellWidth_.size(); }

  /**
   * @brief Return the width of the cells
   */
  const std::vector<int> &cellWidth() const { return cellWidth_; }

  /**
   * @brief Return the target x coordinates for legalization
   */
  const std::vector<int> &cellTargetX() const { return cellTargetX_; }

  /**
   * @brief Return the target y coordinates for legalization
   */
  const std::vector<int> &cellTargetY() const { return cellTargetY_; }

  /**
   * @brief Return the cost model used by the algorithm
   */
  LegalizationModel costModel() const { return costModel_; }

  /**
   * @brief Set the cost model used by the algorithm
   */
  void setCostModel(LegalizationModel m) { costModel_ = m; }

  /**
   * @brief Run the algorithm
   */
  void run();

  /**
   * @brief Compute the x coordinates after legalization
   */
  std::vector<int> cellLegalX() const;

  /**
   * @brief Compute the y coordinates after legalization
   */
  std::vector<int> cellLegalY() const;

  /**
   * @brief Check consistency of the datastructure
   */
  void check() const;

 private:
  /**
   * @brief Place a single cell optimally
   * Return true if successful
   */
  bool placeCellOptimally(int cell);

  /**
   * @brief Simulate placing a single cell in a given row
   * Return a pair: true if successful and the X coordinate
   */
  std::pair<bool, int> placeCellOptimally(int cell, int row) const;

  /**
   * @brief Compute the ordering of the cells
   */
  std::vector<int> computeCellOrder(float weightX, float weightWidth,
                                    float weightY) const;

  /**
   * @brief Returns true if the cell is to be ignored by legalization
   */
  bool isIgnored(int cell) const { return cellWidth_[cell] == -1; }

  /**
   * @brief Materialize the placement of a cell
   */
  void doPlacement(int cell, int row, int x);

  /**
   * @brief Undo the placement of a cell
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
