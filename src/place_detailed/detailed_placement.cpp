#include "place_detailed/detailed_placement.hpp"

#include <algorithm>
#include <stdexcept>

namespace coloquinte {

DetailedPlacement fromIspdCircuit(const Circuit &circuit) {
  // Represent fixed cells with -1 width so they are not considered
  std::vector<int> widths = circuit.cellWidths;
  for (int i = 0; i < circuit.nbCells(); ++i) {
    if (circuit.cellFixed[i]) {
      widths[i] = -1;
    }
  }
  return DetailedPlacement(circuit.computeRows(), widths, circuit.cellX,
                           circuit.cellY);
}

DetailedPlacement::DetailedPlacement(const std::vector<Rectangle> &rows,
                                     const std::vector<int> &width,
                                     const std::vector<int> &posX,
                                     const std::vector<int> &posY) {
  assert(posX.size() == width.size());
  assert(posY.size() == width.size());
  rows_ = rows;
  cellX_ = posX;
  cellY_ = posY;
  cellWidth_ = width;

  std::sort(rows_.begin(), rows_.end(), [](Rectangle a, Rectangle b) -> bool {
    return a.minY < b.minY || (a.minY == b.minY && a.minX < b.minX);
  });

  // Find the cells allocated to each row
  std::vector<std::vector<int> > rowToCells(nbRows());
  for (int i = 0; i < nbCells(); ++i) {
    if (isIgnored(i)) continue;
    int x = posX[i];
    int y = posY[i];
    // Find the first row starting after the cell
    auto it = std::upper_bound(
        rows_.begin(), rows_.end(), Rectangle(x, x, y, y),
        [](Rectangle a, Rectangle b) {
          return a.minY < b.minY || (a.minY == b.minY && a.minX < b.minX);
        });
    if (it == rows_.begin()) {
      throw std::runtime_error("No row found for the cell");
    }
    int row = it - rows_.begin() - 1;
    Rectangle rect = rows_[row];
    if (rect.minY != y) {
      throw std::runtime_error("Found row doesn't have the right y");
    }
    if (rect.minX > x) {
      throw std::runtime_error("Found row starts after the cell");
    }
    if (rect.maxX < x + width[i]) {
      throw std::runtime_error("Found row ends before the cell");
    }
    rowToCells[row].push_back(i);
  }

  // Now sort the cells
  for (std::vector<int> &cells : rowToCells) {
    std::sort(cells.begin(), cells.end(),
              [&posX](int c1, int c2) { return posX[c1] < posX[c2]; });
  }

  cellPred_.assign(nbCells(), -1);
  cellNext_.assign(nbCells(), -1);
  cellRow_.assign(nbCells(), -1);
  rowFirstCell_.assign(nbRows(), -1);
  rowLastCell_.assign(nbRows(), -1);
  // Now setup the cells in the rows
  for (int row = 0; row < nbRows(); ++row) {
    if (rowToCells[row].empty()) continue;
    for (int c : rowToCells[row]) {
      cellRow_[c] = row;
    }
    for (int i = 0; i + 1 < rowToCells[row].size(); ++i) {
      int c1 = rowToCells[row][i];
      int c2 = rowToCells[row][i + 1];
      cellNext_[c1] = c2;
      cellPred_[c2] = c1;
      if (cellX_[c1] + cellWidth_[c1] > cellX_[c2]) {
        throw std::runtime_error("Overlap between cells");
      }
    }
    rowFirstCell_[row] = rowToCells[row].front();
    rowLastCell_[row] = rowToCells[row].back();
  }
}

void DetailedPlacement::check() const {
  if (rows_.size() != nbRows()) throw std::runtime_error("Row size mismatch");
  if (rowFirstCell_.size() != nbRows())
    throw std::runtime_error("Row size mismatch");
  if (rowLastCell_.size() != nbRows())
    throw std::runtime_error("Row size mismatch");
  if (cellWidth_.size() != nbCells())
    throw std::runtime_error("Cell size mismatch");
  if (cellPred_.size() != nbCells())
    throw std::runtime_error("Cell size mismatch");
  if (cellNext_.size() != nbCells())
    throw std::runtime_error("Cell size mismatch");
  if (cellRow_.size() != nbCells())
    throw std::runtime_error("Cell size mismatch");
  if (cellX_.size() != nbCells())
    throw std::runtime_error("Cell size mismatch");
  if (cellY_.size() != nbCells())
    throw std::runtime_error("Cell size mismatch");
  for (int i = 0; i < nbRows(); ++i) {
    int fc = rowFirstCell(i);
    int lc = rowLastCell(i);
    if ((lc == -1) != (fc == -1)) {
      throw std::runtime_error("Inconcistency between first and last cell");
    }
    if (fc == -1) continue;
    if (cellRow(fc) != i) {
      throw std::runtime_error("Inconsistency in the first row cell");
    }
    if (cellPred(fc) != -1) {
      throw std::runtime_error("Inconsistency in the first row cell");
    }
    if (cellRow(lc) != i) {
      throw std::runtime_error("Inconsistency in the last row cell");
    }
    if (cellNext(lc) != -1) {
      throw std::runtime_error("Inconsistency in the last row cell");
    }
  }
  for (int i = 0; i < nbCells(); ++i) {
    int pc = cellPred(i);
    int nc = cellNext(i);
    int row = cellRow(i);
    if (row < -1 || row >= nbRows()) {
      throw std::runtime_error("Invalid row number");
    }
    if (row == -1) {
      if (pc != -1 || nc != -1) {
        throw std::runtime_error(
            "Non-placed cell should have no predecessor/successor");
      }
      continue;
    }
    if (pc != -1) {
      if (cellRow(pc) != row) {
        throw std::runtime_error("Row inconsistency with the predecessor");
      }
      if (cellX(pc) + cellWidth(pc) > cellX(i)) {
        throw std::runtime_error("Overlap with the predecessor");
      }
    } else {
      if (rowFirstCell(row) != i) {
        throw std::runtime_error("Inconsistent first row cell");
      }
      if (cellX(i) < rows_[row].minX) {
        throw std::runtime_error("Element is out of the row");
      }
    }
    if (nc != -1) {
      if (cellRow(nc) != row) {
        throw std::runtime_error("Row inconsistency with the successor");
      }
      if (cellX(i) + cellWidth(i) > cellX(nc)) {
        throw std::runtime_error("Overlap with the successor");
      }
    } else {
      if (rowLastCell(row) != i) {
        throw std::runtime_error("Inconsistent last row cell");
      }
      if (cellX(i) + cellWidth(i) > rows_[row].maxX) {
        throw std::runtime_error("Element is out of the row");
      }
    }
  }
}
}  // namespace coloquinte