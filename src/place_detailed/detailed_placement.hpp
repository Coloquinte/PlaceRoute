#pragma once

#include "coloquinte.hpp"

namespace coloquinte {
/**
 * @brief Representation of the detailed placement state of a standard cells
 * design
 *
 * Each row is represented as a linked lists, so that accessing nearby cells is
 * cheap.
 */
class DetailedPlacement {
 public:
  /**
   * @brief Initialize the datastructure from a circuit
   */
  static DetailedPlacement fromIspdCircuit(const Circuit &circuit);

  /**
   * @brief Initialize the datastructure from a region of a circuit
   */
  static DetailedPlacement fromIspdCircuit(const Circuit &circuit,
                                           const Rectangle &region);

  /**
   * @brief Export the placement obtained to the circuit datastructure
   */
  void exportPlacement(Circuit &circuit);

  /**
   * @brief Create the datastructure
   *      @param rows: Available rows for placement; must all be the right
   * height for the cells
   *      @param width: Width of the cells when placed in a row
   *      @param posX: X coordinate (must be legal)
   *      @param posY: Y coordinate (must be legal)
   */
  static DetailedPlacement fromPos(const std::vector<Row> &rows,
                                   const std::vector<int> &width,
                                   const std::vector<int> &posX,
                                   const std::vector<int> &posY);

  /**
   * @brief Initialize the datastructure
   *      @param rows: Available rows for placement; must all be the right
   * height for the cells
   *      @param width: Width of the cells when placed in a row
   *      @param posX: X coordinate (must be legal)
   *      @param posY: Y coordinate (must be legal)
   *      @param orient: Orientation of the cells (must be legal)
   *      @param pol: Polarity of the cells with respect to the rows
   *      @param cellIndex: Cell index in the original circuit
   */
  DetailedPlacement(const std::vector<Row> &rows, const std::vector<int> &width,
                    const std::vector<int> &posX, const std::vector<int> &posY,
                    const std::vector<CellOrientation> &orient,
                    const std::vector<CellRowPolarity> &pol,
                    const std::vector<int> &cellIndex);

  /**
   * @brief Return the number of rows
   */
  int nbRows() const { return rows_.size(); }

  /**
   * @brief Return the number of cells
   */
  int nbCells() const { return cellWidth_.size(); }

  /**
   * @brief Return all row geometries
   */
  const std::vector<Row> &rows() const { return rows_; }

  /**
   * @brief Returns true if the cell is to be ignored by the detailed placement
   */
  bool isIgnored(int cell) const { return cellWidth_[cell] == -1; }

  /**
   * @brief Returns true if the cell is assigned a placement
   */
  bool isPlaced(int cell) const { return cellRow_[cell] != -1; }

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
   * @brief Return all cells in the row
   */
  std::vector<int> rowCells(int row) const;

  /**
   * @brief Return all cells in the s, sorted by position
   */
  std::vector<int> rowCells(const std::vector<int> &rows) const;

  /**
   * @brief Return all cells in the row strictly between two boundary cells
   */
  std::vector<int> cellsBetween(int row, int cellBefore, int cellAfter) const;

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
   * @brief Return the position of the cell
   */
  Point cellPos(int c) const {
    assert(c >= 0 && c < nbCells());
    return Point(cellX_[c], cellY_[c]);
  }

  /**
   * @brief Return the orientation of the cell
   */
  CellOrientation cellOrientation(int c) const {
    assert(c >= 0 && c < nbCells());
    return cellOrientation_[c];
  }

  /**
   * @brief Return the polarity of the cell
   */
  CellRowPolarity cellRowPolarity(int c) const {
    assert(c >= 0 && c < nbCells());
    return cellRowPolarity_[c];
  }

  /**
   * @brief Return the y position of the row
   */
  int rowY(int r) const {
    assert(r >= 0 && r < nbRows());
    return rows_[r].minY;
  }

  /**
   * @brief Return indices of the cells in the original circuit
   */
  const std::vector<int> &cellIndex() const { return cellIndex_; }

  /**
   * @brief Return the x boundary before the cell (beginning of row or end of
   * previous cell)
   */
  int boundaryBefore(int c) const;

  /**
   * @brief Return the x boundary before the cell (beginning of row or end of
   * previous cell), or the end of the row is cell is -1
   */
  int boundaryBefore(int row, int c) const;

  /**
   * @brief Return the x boundary after the cell (end of row or beginning of
   * next cell)
   */
  int boundaryAfter(int c) const;

  /**
   * @brief Return the x boundary after the cell (end of row or beginning of
   * next cell), or the beginning of the row is cell is -1
   */
  int boundaryAfter(int row, int c) const;

  /**
   * @brief Return the beginning position of the placement site after this cell
   * on this row
   */
  int siteBegin(int row, int pred) const;

  /**
   * @brief Return the ending position of the placement site after this cell on
   * this row
   */
  int siteEnd(int row, int pred) const;

  /**
   * @brief Return true if it is possible to place the cell here
   */
  bool canPlace(int c, int row, int pred, int x) const;

  /**
   * @brief Return true if it is possible to insert the cell with this
   * predecessor
   */
  bool canInsert(int c, int row, int pred) const;

  /**
   * @brief Return the position after inserting the cell here
   */
  Point positionOnInsert(int c, int row, int pred) const;

  /**
   * @brief Return true if it is possible to swap the two cells
   */
  bool canSwap(int c1, int c2) const;

  /**
   * @brief Return the positions after swapping the two cells
   */
  std::pair<Point, Point> positionsOnSwap(int c1, int c2) const;

  /**
   * @brief Do placement of a single cell
   */
  void place(int c, int row, int pred, int x);

  /**
   * @brief Undo placement of a single cell
   */
  void unplace(int c);

  /**
   * @brief Insert the cell with this predecessor
   */
  void insert(int c, int row, int pred);

  /**
   * @brief Insert the cell at this position
   */
  void insertAt(int c, int row, int pred, int x);

  /**
   * @brief Swap the two cells
   */
  void swap(int c1, int c2);

  /**
   * @brief Swap the two cells with the given x position
   */
  void swapAt(int c1, int c2, int x1, int x2);

  /**
   * @brief Run the algorithm
   */
  void run();

  /**
   * @brief Check consistency of the datastructure
   */
  void check() const;

 private:
  std::vector<Row> rows_;
  std::vector<int> rowFirstCell_;
  std::vector<int> rowLastCell_;
  std::vector<int> cellWidth_;
  std::vector<int> cellPred_;
  std::vector<int> cellNext_;
  std::vector<int> cellRow_;
  std::vector<int> cellX_;
  std::vector<int> cellY_;
  std::vector<CellOrientation> cellOrientation_;
  std::vector<CellRowPolarity> cellRowPolarity_;
  std::vector<int> cellIndex_;

  friend class DetailedPlacer;
};
}  // namespace coloquinte
