#include "place_detailed/detailed_placement.hpp"

#include <algorithm>
#include <stdexcept>

namespace coloquinte {

DetailedPlacement DetailedPlacement::fromIspdCircuit(const Circuit &circuit) {
  // Represent fixed cells with -1 width so they are not considered
  std::vector<int> widths = circuit.cellWidth_;
  for (int i = 0; i < circuit.nbCells(); ++i) {
    if (circuit.cellIsFixed_[i]) {
      widths[i] = -1;
    }
  }
  std::vector<int> cellIndex;
  cellIndex.reserve(circuit.nbCells());

  for (int i = 0; i < circuit.nbCells(); ++i) {
    cellIndex.push_back(i);
  }
  // TODO: check that the orientations of the cells are compatible with the rows
  return DetailedPlacement(circuit.computeRows(), widths, circuit.cellX_,
                           circuit.cellY_, cellIndex);
}

DetailedPlacement DetailedPlacement::fromIspdCircuit(const Circuit &circuit,
                                                     const Rectangle &region) {
  // Compute the cells in the placement region
  std::vector<int> cellIndex;
  std::vector<Rectangle> obstacles;
  for (int c = 0; c < circuit.nbCells(); ++c) {
    if (circuit.isFixed(c)) {
      continue;
    }
    Rectangle pl = circuit.placement(c);
    if (region.contains(pl)) {
      cellIndex.push_back(c);
    } else if (region.intersects(pl)) {
      // This will cover some of the rows: we need to reduce their size
      obstacles.push_back(pl);
    }
  }

  // Compute the rows in the placement region
  std::vector<Rectangle> rows;
  for (Rectangle row : circuit.computeRows(obstacles)) {
    // Height must be completely contained in the region
    if (row.minY < region.minY) {
      continue;
    }
    if (row.maxY > region.maxY) {
      continue;
    }
    // Laterally we just need to have an intersection
    if (row.minX >= region.maxX) {
      continue;
    }
    if (row.maxX <= region.minX) {
      continue;
    }
    Rectangle actual(std::max(row.minX, region.minX),
                     std::min(row.maxX, region.maxX), row.minY, row.maxY);
    rows.push_back(actual);
  }

  // Setup the widths and positions
  std::vector<int> widths(cellIndex.size());
  std::vector<int> cellX(cellIndex.size());
  std::vector<int> cellY(cellIndex.size());
  for (int i = 0; i < cellIndex.size(); ++i) {
    int c = cellIndex[i];
    widths[i] = circuit.cellWidth()[c];
    cellX[i] = circuit.cellX()[c];
    cellY[i] = circuit.cellY()[c];
  }
  // Additional fixed cell
  cellIndex.push_back(-1);
  widths.push_back(-1);
  cellX.push_back(0);
  cellY.push_back(0);

  return DetailedPlacement(rows, widths, cellX, cellY, cellIndex);
}

void DetailedPlacement::exportPlacement(Circuit &circuit) {
  for (int i = 0; i < nbCells(); ++i) {
    int cell = cellIndex_[i];
    if (cell < 0) {
      continue;
    }
    if (circuit.isFixed(cell)) {
      continue;
    }
    circuit.cellX_[cell] = cellX(i);
    circuit.cellY_[cell] = cellY(i);
  }
}

DetailedPlacement::DetailedPlacement(const std::vector<Rectangle> &rows,
                                     const std::vector<int> &width,
                                     const std::vector<int> &posX,
                                     const std::vector<int> &posY,
                                     const std::vector<int> &cellIndex) {
  assert(posX.size() == width.size());
  assert(posY.size() == width.size());
  assert(cellIndex.size() == width.size());
  rows_ = rows;
  cellX_ = posX;
  cellY_ = posY;
  cellWidth_ = width;
  cellIndex_ = cellIndex;

  std::sort(rows_.begin(), rows_.end(), [](Rectangle a, Rectangle b) -> bool {
    return a.minY < b.minY || (a.minY == b.minY && a.minX < b.minX);
  });

  // Find the cells allocated to each row
  std::vector<std::vector<int> > rowToCells(nbRows());
  for (int i = 0; i < nbCells(); ++i) {
    if (isIgnored(i)) {
      continue;
    }
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
    if (rowToCells[row].empty()) {
      continue;
    }
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

  check();
}

std::vector<int> DetailedPlacement::rowCells(int row) const {
  std::vector<int> ret;
  for (int c = rowFirstCell(row); c != -1; c = cellNext(c)) {
    ret.push_back(c);
  }
  return ret;
}

std::vector<int> DetailedPlacement::rowCells(
    const std::vector<int> &rows) const {
  std::vector<std::pair<int, int> > sortedCells;
  for (int row : rows) {
    for (int c : rowCells(row)) {
      sortedCells.emplace_back(cellX(c), c);
    }
  }
  std::stable_sort(sortedCells.begin(), sortedCells.end());
  std::vector<int> cells;
  cells.reserve(sortedCells.size());

  for (auto p : sortedCells) {
    cells.push_back(p.second);
  }
  return cells;
}

std::vector<int> DetailedPlacement::cellsBetween(int row, int cellBefore,
                                                 int cellAfter) const {
  if (cellBefore != -1 && cellRow(cellBefore) != row) {
    throw std::runtime_error("First cell row is inconsistent");
  }
  if (cellAfter != -1 && cellRow(cellAfter) != row) {
    throw std::runtime_error("Last cell row is inconsistent");
  }
  int firstCell = cellBefore == -1 ? rowFirstCell(row) : cellNext(cellBefore);
  std::vector<int> ret;
  for (int c = firstCell; c != cellAfter; c = cellNext(c)) {
    if (c == -1) {
      throw std::runtime_error("End range cell was not encountered");
    }
    ret.push_back(c);
  }
  return ret;
}

int DetailedPlacement::boundaryBefore(int c) const {
  assert(isPlaced(c));
  int pred = cellPred(c);
  if (pred == -1) {
    return rows_[cellRow(c)].minX;
  }
  return cellX(pred) + cellWidth(pred);
}

int DetailedPlacement::boundaryAfter(int c) const {
  assert(isPlaced(c));
  int next = cellNext(c);
  if (next == -1) {
    return rows_[cellRow(c)].maxX;
  }
  return cellX(next);
}

int DetailedPlacement::boundaryBefore(int row, int c) const {
  if (c == -1) {
    return rows_[row].maxX;
  }
  return boundaryBefore(c);
}

int DetailedPlacement::boundaryAfter(int row, int c) const {
  if (c == -1) {
    return rows_[row].minX;
  }
  return boundaryAfter(c);
}

int DetailedPlacement::siteBegin(int row, int pred) const {
  return pred == -1 ? rows_[row].minX : cellX(pred) + cellWidth(pred);
}

int DetailedPlacement::siteEnd(int row, int pred) const {
  int next = pred == -1 ? rowFirstCell(row) : cellNext(pred);
  return next == -1 ? rows_[row].maxX : cellX(next);
}

bool DetailedPlacement::canPlace(int c, int row, int pred, int x) const {
  if (isPlaced(c)) {
    throw std::runtime_error("Cannot attempt to place already placed cell");
  }
  return x >= siteBegin(row, pred) && x + cellWidth(c) <= siteEnd(row, pred);
}

bool DetailedPlacement::canInsert(int c, int row, int pred) const {
  if (!isPlaced(c)) {
    throw std::runtime_error(
        "Cannot attempt to insert a cell that is not placed yet");
  }
  if (c == pred) {
    // Do not insert after itself
    return false;
  }
  if (cellRow(c) == row && cellPred(c) == pred) {
    // Do not insert before itself
    return false;
  }
  return siteEnd(row, pred) - siteBegin(row, pred) >= cellWidth(c);
}

bool DetailedPlacement::canSwap(int c1, int c2) const {
  if (!isPlaced(c1) || !isPlaced(c2)) {
    throw std::runtime_error("Cannot swap cells that are not placed yet");
  }
  if (c1 == c2) {
    // Do not swap a cell with itself
    return false;
  }
  if (cellPred(c1) == c2 || cellPred(c2) == c1) {
    // We can always swap neighbours
    return true;
  }  // Otherwise check if there is enough space for both cells

  int b1 = boundaryBefore(c1);
  int b2 = boundaryBefore(c2);
  int e1 = boundaryAfter(c1);
  int e2 = boundaryAfter(c2);
  return e2 - b2 >= cellWidth(c1) && e1 - b1 >= cellWidth(c2);
}

void DetailedPlacement::place(int c, int row, int pred, int x) {
  if (!canPlace(c, row, pred, x)) {
    throw std::runtime_error("Cannot place the cell");
  }
  cellRow_[c] = row;
  int next = pred == -1 ? rowFirstCell(row) : cellNext(pred);
  if (pred == -1) {
    rowFirstCell_[row] = c;
  } else {
    cellNext_[pred] = c;
  }
  cellPred_[c] = pred;
  if (next == -1) {
    rowLastCell_[row] = c;
  } else {
    cellPred_[next] = c;
  }
  cellNext_[c] = next;
  cellX_[c] = x;
  cellY_[c] = rows_[row].minY;
}

void DetailedPlacement::unplace(int c) {
  int row = cellRow(c);
  int pred = cellPred(c);
  int next = cellNext(c);
  cellRow_[c] = -1;
  if (pred == -1) {
    rowFirstCell_[row] = next;
  } else {
    cellNext_[pred] = next;
  }
  cellPred_[c] = -1;
  if (next == -1) {
    rowLastCell_[row] = pred;
  } else {
    cellPred_[next] = pred;
  }
  cellNext_[c] = -1;
}

void DetailedPlacement::insert(int c, int row, int pred) {
  if (!canInsert(c, row, pred)) {
    throw std::runtime_error("Cannot insert this cell here");
  }
  Point pos = positionOnInsert(c, row, pred);
  unplace(c);
  place(c, row, pred, pos.x);
}

void DetailedPlacement::swap(int c1, int c2) {
  if (!canSwap(c1, c2)) {
    throw std::runtime_error("Cannot swap these cells");
  }
  auto [pos1, pos2] = positionsOnSwap(c1, c2);
  int r1 = cellRow(c1);
  int r2 = cellRow(c2);
  int p1 = cellPred(c1);
  int p2 = cellPred(c2);
  int x1 = pos1.x;
  int x2 = pos2.x;
  unplace(c1);
  unplace(c2);
  if (p1 == c2) {
    place(c1, r2, p2, x1);
    place(c2, r1, c1, x2);
  } else if (p2 == c1) {
    place(c2, r1, p1, x2);
    place(c1, r2, c2, x1);
  } else {
    place(c1, r2, p2, x1);
    place(c2, r1, p1, x2);
  }
}

std::pair<Point, Point> DetailedPlacement::positionsOnSwap(int c1,
                                                           int c2) const {
  Point p1 = cellPos(c1);
  Point p2 = cellPos(c2);
  int x1, x2;
  if (cellPred(c1) == c2) {
    x1 = p2.x;
    x2 = p2.x + cellWidth(c1);
  } else if (cellPred(c2) == c1) {
    x2 = p1.x;
    x1 = p1.x + cellWidth(c2);
  } else {
    x1 = (boundaryBefore(c2) + boundaryAfter(c2) - cellWidth(c1)) / 2;
    x2 = (boundaryBefore(c1) + boundaryAfter(c1) - cellWidth(c2)) / 2;
  }
  return std::make_pair(Point(x1, p2.y), Point(x2, p1.y));
}

Point DetailedPlacement::positionOnInsert(int c, int row, int pred) const {
  int x = (siteEnd(row, pred) - cellWidth(c) + siteBegin(row, pred)) / 2;
  int y = rowY(row);
  return Point(x, y);
}

void DetailedPlacement::check() const {
  if (rows_.size() != nbRows()) {
    throw std::runtime_error("Row size mismatch");
  }
  if (rowFirstCell_.size() != nbRows()) {
    throw std::runtime_error("Row size mismatch");
  }
  if (rowLastCell_.size() != nbRows()) {
    throw std::runtime_error("Row size mismatch");
  }
  if (cellWidth_.size() != nbCells()) {
    throw std::runtime_error("Cell size mismatch");
  }
  if (cellPred_.size() != nbCells()) {
    throw std::runtime_error("Cell size mismatch");
  }
  if (cellNext_.size() != nbCells()) {
    throw std::runtime_error("Cell size mismatch");
  }
  if (cellRow_.size() != nbCells()) {
    throw std::runtime_error("Cell size mismatch");
  }
  if (cellX_.size() != nbCells()) {
    throw std::runtime_error("Cell size mismatch");
  }
  if (cellY_.size() != nbCells()) {
    throw std::runtime_error("Cell size mismatch");
  }
  if (cellIndex_.size() != nbCells()) {
    throw std::runtime_error("Cell size mismatch");
  }
  for (int i = 0; i < nbRows(); ++i) {
    int fc = rowFirstCell(i);
    int lc = rowLastCell(i);
    if ((lc == -1) != (fc == -1)) {
      throw std::runtime_error("Inconcistency between first and last cell");
    }
    if (fc == -1) {
      continue;
    }
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