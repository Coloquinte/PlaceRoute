#pragma once

#include <queue>
#include <vector>

namespace coloquinte {
/**
 * @brief Obtain the positions minimizing total weighted displacement along a
 * row.
 *
 * It is an ordered single row problem/fixed order single machine scheduling
 * problem, solved by the specialized cascading descent algorithm.
 *
 * The cost model is linear in the distance to the target position, weighted by
 * the width of the cells
 **/
class RowLegalizer {
 public:
  /// Initialize
  RowLegalizer(int b, int e) : begin_(b), end_(e), cumWidth_(1, 0) {}

  /**
   * @brief Return the space already used in the row
   */
  int usedSpace() const { return cumWidth_.back(); }

  /**
   * @brief Return the space remaining in the row
   */
  int remainingSpace() const { return end_ - begin_ - usedSpace(); }

  /**
   * @brief Return the current last position
   */
  int lastAvailablePos() const { return constrainingPos_.back() + usedSpace(); }

  /**
   * @brief Return the cost of pushing a new cell, without updating the
   * datastructure
   */
  long long getCost(int width, int targetPos);

  /**
   * @brief Update the datastructure with a new cell and return the cost
   */
  long long push(int width, int targetPos);

  /**
   * @brief Return the placement of each cell in the datastructure
   */
  std::vector<int> getPlacement() const;

  /**
   * @brief Check the consistency
   */
  void check() const;

  /**
   * @brief Remove all cells from the datastructure
   */
  void clear();

 private:
  /**
   * @brief Representation of the events in the cascading descent algorithm
   */
  struct Bound {
    /// Will be the target absolute position of the cell
    int absolutePos;
    /// Will be proportional to the width of the cell
    int weight;

    bool operator<(Bound const o) const {
      return absolutePos < o.absolutePos ||
             (absolutePos == o.absolutePos && weight < o.weight);
    }
    Bound(int w, int absPos) : absolutePos(absPos), weight(w) {}
  };

  /**
   * @brief Return the number of elements
   */
  int nbElements() const { return cumWidth_.size() - 1; }

  /**
   * @brief Return the width of this element
   */
  int width(int ind) const { return cumWidth_[ind + 1] - cumWidth_[ind]; }

  /// Get the cost of pushing a cell on the row
  long long getDisplacement(int width, int targetPos, bool update);

  /// Leftmost coordinate of the region
  int begin_;

  /// Rightmost coordinate of the region
  int end_;

  /// Where the cells constrain the positions of preceding cells
  std::vector<int> constrainingPos_;

  /// Cumulative width of the cells
  std::vector<int> cumWidth_;

  /// Priority queue for the cascading descent algorithm
  std::priority_queue<Bound> bounds;
};
}  // namespace coloquinte