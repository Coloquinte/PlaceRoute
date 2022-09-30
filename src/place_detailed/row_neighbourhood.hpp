#pragma once

#include "coloquinte.hpp"

namespace coloquinte {
/**
 * @brief Given standard cell rows, access neighbouring rows easily
 */
class RowNeighbourhood {
 public:
  /**
   * @brief Initialize the datastructure
   */
  explicit RowNeighbourhood(const std::vector<Rectangle> &rows,
                            int nbNeighbourRows = 1);

  /**
   * @brief Return the number of rows
   */
  int nbRows() const { return rowsBelow_.size(); }

  /**
   * @brief Test if r1 is below r2
   */
  static bool isBelow(Rectangle r1, Rectangle r2);

  /**
   * @brief Test if r1 is above r2
   */
  static bool isAbove(Rectangle r1, Rectangle r2);

  /**
   * @brief Test if r1 is left of r2
   */
  static bool isLeft(Rectangle r1, Rectangle r2);

  /**
   * @brief Test if r1 is right of r2
   */
  static bool isRight(Rectangle r1, Rectangle r2);

  /**
   * @brief Return a few rows below, closest ones first
   */
  const std::vector<int> &rowsBelow(int row) const { return rowsBelow_[row]; }

  /**
   * @brief Return a few rows above, closest ones first
   */
  const std::vector<int> &rowsAbove(int row) const { return rowsAbove_[row]; }

  /**
   * @brief Return a few rows on the left, closest ones first
   */
  const std::vector<int> &rowsLeft(int row) const { return rowsLeft_[row]; }

  /**
   * @brief Return a few rows on the right, closest ones first
   */
  const std::vector<int> &rowsRight(int row) const { return rowsRight_[row]; }

  /**
   * @brief Check consistency of the datastructure
   */
  void check() const;

 private:
  /**
   * @brief Simple unoptimized setup with quadratic complexity
   *
   * TODO: Implement a setup method with lower complexity
   */
  void simpleSetup(const std::vector<Rectangle> &rows, int nbNeighbourRows);

  /**
   * @brief Build the structure indicating the rows below
   */
  static std::vector<std::vector<int> > rowsBelow(
      const std::vector<Rectangle> &rows, int nbNeighbourRows);

  /**
   * @brief Build the structure indicating the rows above
   */
  static std::vector<std::vector<int> > rowsAbove(
      const std::vector<Rectangle> &rows, int nbNeighbourRows);

  /**
   * @brief Build the structure indicating the rows on the sides (assumes rows
   * above/below are built already)
   */
  void buildRowsSides(const std::vector<Rectangle> &rows, int nbNeighbourRows);

  std::vector<int> buildLeftFrom(Rectangle row, const std::vector<Rectangle> &rows, int ind);
  std::vector<int> buildRightFrom(Rectangle row, const std::vector<Rectangle> &rows, int ind);

 private:
  std::vector<std::vector<int> > rowsBelow_;
  std::vector<std::vector<int> > rowsAbove_;
  std::vector<std::vector<int> > rowsLeft_;
  std::vector<std::vector<int> > rowsRight_;
};
}  // Namespace coloquinte