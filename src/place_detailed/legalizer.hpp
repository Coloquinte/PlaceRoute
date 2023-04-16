#pragma once

#include "coloquinte.hpp"
#include "place_detailed/row_legalizer.hpp"

namespace coloquinte {
/**
 * @brief Algorithms to obtain a legal placement for standard cells
 *
 * TODO: Create a separate template class to handle more cost models with
 * Abacus-like legalization
 * TODO: Create a separate class to handle backtracking during legalization for
 * hard cases
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
   *      @param width: Width of the cells
   *      @param height: Height of the cells
   *      @param polarity: Polarity of the cells with respect to rows
   *      @param targetX: Target x coordinate for legalization
   *      @param targetY: Target y coordinate for legalization
   */
  Legalizer(const std::vector<Row> &rows, const std::vector<int> &width,
            const std::vector<int> &height,
            const std::vector<CellRowPolarity> &polarity,
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
   * @brief Return the sum of the widths of the cells
   */
  long long totalCellWidth() const;

  /**
   * @brief Run the algorithm
   */
  void run(const ColoquinteParameters &parameters);

  /**
   * @brief Compute the x coordinates after legalization
   */
  std::vector<int> cellLegalX() const;

  /**
   * @brief Compute the y coordinates after legalization
   */
  std::vector<int> cellLegalY() const;

  /**
   * @brief Compute the cell orientation after legalization
   */
  std::vector<CellOrientation> cellLegalOrientation() const;

  /**
   * @brief Check consistency of the datastructure
   */
  void check() const;

  /**
   * @brief Report on the datastructure on stdout
   */
  void report(bool verbose = false) const;

 private:
  /**
   * @brief Return all distances with a given cost model
   */
  std::vector<float> allDistances(LegalizationModel model) const;

  /**
   * @brief Place a single cell optimally
   * Return true if successful
   */
  bool placeCellOptimally(int cell, LegalizationModel costModel);

  /**
   * @brief Simulate placing a single cell in a given row, pushing other cells
   * as needed Return a pair: true if successful and the added distance
   */
  std::pair<bool, long long> placeCellOptimally(int cell, int row);

  /**
   * @brief Compute the ordering of the cells
   *
   * @param weightX Weight allocated to the target x coordinate; main ordering
   * component, should usually be 1
   * @param weightWidth Weight allocated to the width; allows ordering by left
   * side (0), center (0.5) or right side (1)
   * @param weightY Weight allocated to the target y coordinate; allows diagonal
   * ordering, should usually be close to 0
   *
   * @return An ordering of the cell indices
   */
  std::vector<int> computeCellOrder(float weightX, float weightWidth,
                                    float weightY) const;

  /**
   * @brief Returns true if the cell is already placed by the algorithm
   */
  bool isPlaced(int cell) const { return cellToRow_[cell] != -1; }

  /**
   * @brief Find the row that is closest to the target position
   */
  int closestRow(int y) const;

 private:
  // Placement data
  std::vector<Row> rows_;
  std::vector<int> cellWidth_;
  std::vector<int> cellHeight_;
  std::vector<CellRowPolarity> cellRowPolarity_;
  std::vector<int> cellTargetX_;
  std::vector<int> cellTargetY_;

  // Placement status
  std::vector<std::vector<int> > rowToCells_;
  std::vector<RowLegalizer> rowLegalizers_;
  std::vector<int> cellToRow_;
  std::vector<int> cellToX_;
  std::vector<int> cellToY_;
};
}  // namespace coloquinte