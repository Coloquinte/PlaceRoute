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
  static void place(Circuit &circuit, int effort) {
    place(circuit, DetailedPlacerParameters(effort));
  }

  /**
   * @brief Run detailed placement on the circuit representation
   *
   * @param circuit The circuit to be modified
   * @param params Placement parameters
   */
  static void place(Circuit &circuit, const DetailedPlacerParameters &params);

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
   * @brief Run a simple optimization using only cell insertion
   *
   * @param nbRows Number of neighbouring rows to look at for each row
   * @param nbNeighbours Number of closest neighbours to consider for each cell
   */
  void runInserts(int nbRows, int nbNeighbours);

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
   * @brief Run the cell insertion optimization within a row
   *
   * @param row Row to look for insertions
   * @param nbNeighbours Number of closest neighbours to consider for each cell
   */
  void runInsertsOneRow(int row, int nbNeighbours);

  /**
   * @brief Run the cell swapping optimization on the two rows
   *
   * @param r1 First row to look for swaps
   * @param r2 Second row to look for swaps
   * @param nbNeighbours Number of closest neighbours to consider for each cell
   */
  void runSwapsTwoRows(int r1, int r2, int nbNeighbours);

  /**
   * @brief Run the cell insertion optimization on the two rows
   *
   * @param r1 Row to look for cells
   * @param r2 Row to insert into
   * @param nbNeighbours Number of closest neighbours to consider for each cell
   */
  void runInsertsTwoRows(int r1, int r2, int nbNeighbours);

  /**
   * @brief Run the cell swapping optimization on the two rows
   *
   * After a swap is found, retry to find other improvements
   *
   * @param r1 First row to look for swaps
   * @param r2 Second row to look for swaps
   * @param nbNeighbours Number of closest neighbours to consider for each cell
   */
  void runSwapsTwoRowsAmplify(int r1, int r2, int nbNeighbours);

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
   * @brief Given two ordered rows, obtain the index of the closest cell in the
   * second row for each cell in the firt row
   */
  std::vector<int> computeClosestIndexInRow(
      const std::vector<int> &row1Cells,
      const std::vector<int> &row2Cells) const;

  /**
   * @brief Swap the two cells
   */
  void doSwap(int c1, int c2);

  /**
   * @brief Insert the cell here
   */
  void doInsert(int c, int row, int pred);

  /**
   * @brief Try to swap two cells
   *
   * @return True on a change
   */
  bool trySwap(int c1, int c2) { return bestSwap(c1, {c2}); }

  /**
   * @brief Insert the cell at a given position
   *
   * @return True on a change
   */
  bool tryInsert(int c, int row, int pred) {
    return bestInsert(c, row, {pred});
  }

  /**
   * @brief Perform the best swap out of many candidates
   *
   * @return True on a change
   */
  bool bestSwap(int c, const std::vector<int> &candidates);

  /**
   * @brief Perform the best insert out of many candidates
   *
   * @return True on a change
   */
  bool bestInsert(int c, int row, const std::vector<int> &candidates);

  /**
   * @brief Perfom the best improvement swap between a cell and another row;
   * update the inputs if necessary on a swap
   *
   * @return True on a change
   */
  bool bestSwapUpdate(int &c, int &from, int nbNeighbours);

  /**
   * @brief Return the feasibility of the swap, and the new value if feasible
   */
  std::pair<bool, long long> valueOnSwap(int c1, int c2);

  /**
   * @brief Return the feasibility of the insertion, and the new value if
   * feasible
   */
  std::pair<bool, long long> valueOnInsert(int c, int row, int pred);

  /**
   * @brief Return the first cell in a row that has x larger or equal to the
   * target cell
   */
  int findCellAfter(int target, int fromCell) const;

  /**
   * @brief Return the first cell in a row that has x smaller or equal to the
   * target cell
   */
  int findCellBefore(int target, int fromCell) const;

 private:
  /**
   * @brief Update the cell position in the objective after a change in the
   * placement
   */
  void updateCellPos(int c);

  /**
   * @brief Update the cell position in the objective (not the placement)
   */
  void updateCellPos(int c, Point pos);

 private:
  DetailedPlacement placement_;
  IncrNetModel xtopo_;
  IncrNetModel ytopo_;
};
}  // namespace coloquinte
