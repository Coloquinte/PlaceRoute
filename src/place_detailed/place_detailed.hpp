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

 private:
  /**
   * @brief Initialize the datastructure
   */
  explicit DetailedPlacer(const Circuit &circuit);

  /**
   * @brief Run a simple optimization using only cell swapping
   */
  void runSwaps();

  /**
   * @brief Attempt to swap the two cells; keep the modification if it improves
   * the result
   */
  bool trySwap(int c1, int c2);

  /**
   * @brief Attempt to insert the cell here; keep the modification if it
   * improves the result
   */
  bool tryInsert(int c, int row, int pred);

  /**
   * @brief Attempt to swap the two cells and shift them; keep the modification
   * if it improves the result
   */
  bool trySwapShift(int c1, int c2);

  /**
   * @brief Change the cell coordinates to optimize the wirelength without
   * reordering them
   */
  void optimizeShift(const std::vector<int> &cells);

  /**
   * @brief Update the cell position in the objective after a change in the
   * placement
   */
  void updateCellPos(int c);

  /**
   * @brief Return the current objective value
   */
  long long value() const { return xtopo_.value() + ytopo_.value(); }

 private:
  DetailedPlacement placement_;
  IncrNetModel xtopo_;
  IncrNetModel ytopo_;
};
}  // namespace coloquinte
