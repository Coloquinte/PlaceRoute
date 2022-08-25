#pragma once

#include "coloquinte.hpp"
#include "place_detailed/detailed_placement.hpp"
#include "place_detailed/incr_net_model.hpp"

namespace coloquinte {
/**
 * @brief Main class for detailed placement
 */
class DetailedPlacer {
 public:
  /**
   * @brief Run detailed placement on the circuit representation
   *
   * @param circuit The circuit to be modified
   * @param effort Effort level, between 0 and 9
   */
  static void place(Circuit &circuit, int effort);

  /**
   * @brief Initialize the datastructure
   */
  explicit DetailedPlacer(const Circuit &circuit);

  /**
   * @brief Check the consistency of the datastructure
   */
  void check() const;

  /**
   * @brief Run a simple optimization using only cell swapping
   *
   * @param nbRows Number of neighbouring rows to look at for each row
   * @param nbNeighbours Number of closest neighbours to consider for each cell
   */
  void runSwaps(int nbRows, int nbNeighbours);

  /**
   * @brief Run a simple optimization using only cell shifting (no reordering)
   *
   * @param nbRows Number of neighbouring rows considered
   */
  void runShifts(int nbRows);

  /**
   * @brief Run the cell swapping optimization within a row
   *
   * @param row Row to look for swaps
   * @param nbNeighbours Number of closest neighbours to consider for each cell
   */
  void runSwapsOneRow(int row, int nbNeighbours);

  /**
   * @brief Run the cell swapping optimization on the two rows
   *
   * @param r1 First row to look for swaps
   * @param r2 Second row to look for swaps
   * @param nbNeighbours Number of closest neighbours to consider for each cell
   */
  void runSwapsTwoRows(int r1, int r2, int nbNeighbours);

  /**
   * @brief Attempt to swap the two cells; keep the modification if it improves
   * the result strictly
   */
  bool trySwap(int c1, int c2);

  /**
   * @brief Attempt to insert the cell here; keep the modification if it
   * improves the result strictly
   */
  bool tryInsert(int c, int row, int pred);

  /**
   * @brief Change the cell coordinates to optimize the wirelength without
   * reordering them
   */
  void optimizeShift(const std::vector<int> &cells);

  /**
   * @brief Return the current objective value
   */
  long long value() const { return xtopo_.value() + ytopo_.value(); }

 private:
  /**
   * @brief Update the cell position in the objective after a change in the
   * placement
   */
  void updateCellPos(int c);

 private:
  DetailedPlacement placement_;
  IncrNetModel xtopo_;
  IncrNetModel ytopo_;
};
}  // namespace coloquinte
