#pragma once

#include "coloquinte.hpp"

namespace coloquinte {
/**
 * @brief Representation of a detailed placement of standard cells
 *
 */
class DetailedPlacement {
 public:
  /**
   * @brief Initialize the datastructure from a circuit
   */
  static DetailedPlacement fromIspdCircuit(const Circuit &circuit);

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
  DetailedPlacement(const std::vector<Rectangle> &rows,
                 const std::vector<int> &width, const std::vector<int> &targetX,
                 const std::vector<int> &targetY);

  /**
   * @brief Return the number of rows
   */
  int nbRows() const { return rows_.size(); }

  /**
   * @brief Return the number of cells
   */
  int nbCells() const { return cellWidth_.size(); }

  /**
   * @brief Returns true if the cell is to be ignored by the detailed placement
   */
  bool isIgnored(int cell) const { return cellWidth_[cell] == -1; }

  /**
   * @brief Return the width of the cell
   */
  int cellWidth(int cell) const { return cellWidth_[cell]; }

  /**
   * @brief Return the row the cell is currently allocated to, -1 if it is not
   * placed
   */
  int cellRow(int c) const {
    assert(c >= 0 && c < nbCells());
    return cellRow_[c];
  }
  /**
   * @brief Return the predecessor of the cell in its row, -1 if it is the first
   */
  int cellPred(int c) const {
    assert(c >= 0 && c < nbCells());
    return cellPred_[c];
  }

  /**
   * @brief Return the successor of the cell in its row, -1 if it is the last
   */
  int cellNext(int c) const {
    assert(c >= 0 && c < nbCells());
    return cellNext_[c];
  }

  /**
   * @brief Return the first cell in the row, -1 if there is none
   */
  int rowFirstCell(int row) const {
    assert(row >= 0 && row < nbRows());
    return rowFirstCell_[row];
  }

  /**
   * @brief Return the last cell in the row, -1 if there is none
   */
  int rowLastCell(int row) const {
    assert(row >= 0 && row < nbRows());
    return rowLastCell_[row];
  }

  /**
   * @brief Return the x position of the cell
   */
  int cellX(int c) const {
    assert(c >= 0 && c < nbCells());
    return cellX_[c];
  }

  /**
   * @brief Return the y position of the cell
   */
  int cellY(int c) const {
    assert(c >= 0 && c < nbCells());
    return cellY_[c];
  }

  /**
   * @brief Place a single cell
   */
  void place(int c, int row, int pred, int next, int x);

  /**
   * @brief Remove the placement of the cell
   */
  void remove(int c);

  /**
   * @brief Run the algorithm
   */
  void run();

  /**
   * @brief Check consistency of the datastructure
   */
  void check() const;

 private:
 private:
  std::vector<Rectangle> rows_;
  std::vector<int> rowFirstCell_;
  std::vector<int> rowLastCell_;
  std::vector<int> cellWidth_;
  std::vector<int> cellPred_;
  std::vector<int> cellNext_;
  std::vector<int> cellRow_;
  std::vector<int> cellX_;
  std::vector<int> cellY_;
};
}