
#include "place_detailed/legalizer.hpp"

#include <algorithm>
#include <cassert>
#include <iostream>
#include <limits>

Legalizer Legalizer::fromIspdCircuit(const Circuit &circuit) {
  Legalizer ret = Legalizer(circuit.computeRows(), circuit.cellWidths,
                            circuit.cellX, circuit.cellY);
  // Represent fixed cells with -1 width so they are not considered
  for (int i = 0; i < circuit.nbCells(); ++i) {
    if (circuit.cellFixed[i]) {
      ret.cellWidth_[i] = -1;
    }
  }
  return ret;
}

Legalizer::Legalizer(const std::vector<Rectangle> &rows,
                     const std::vector<int> &width,
                     const std::vector<int> &targetX,
                     const std::vector<int> &targetY)
    : costModel_(LegalizationModel::L1),
      cellWidth_(width),
      cellTargetX_(targetX),
      cellTargetY_(targetY) {
  assert(width.size() == targetX.size());
  assert(width.size() == targetY.size());
  // Sort the rows
  rows_ = rows;
  std::stable_sort(
      rows_.begin(), rows_.end(), [](Rectangle a, Rectangle b) -> bool {
        return a.minY < b.minY || (a.minY == b.minY && a.minX < b.minX);
      });
  rowToCells_.resize(rows.size());
  cellToX_ = cellTargetX_;
  cellToY_ = cellTargetY_;
  cellToRow_.assign(width.size(), -1);
}

void Legalizer::run() {
  std::vector<int> cellOrder = computeCellOrder(1.0, 0.5, 0.01);

  for (int c : cellOrder) {
    placeCellOptimally(c);
  }
}

void Legalizer::check() const {
  if (cellWidth_.size() != nbCells()) {
    throw std::runtime_error("Number of cell widths does not match");
  }
  if (cellTargetX_.size() != nbCells()) {
    throw std::runtime_error("Number of cell x targets does not match");
  }
  if (cellTargetY_.size() != nbCells()) {
    throw std::runtime_error("Number of cell y targets does not match");
  }
  if (cellToX_.size() != nbCells()) {
    throw std::runtime_error("Number of cell x positions does not match");
  }
  if (cellToY_.size() != nbCells()) {
    throw std::runtime_error("Number of cell y positions does not match");
  }
  if (cellToRow_.size() != nbCells()) {
    throw std::runtime_error("Number of cell row positions does not match");
  }
  if (rowToCells_.size() != nbRows()) {
    throw std::runtime_error("Number of row cells does not match");
  }
  for (int i = 0; i < nbRows(); ++i) {
    for (int c : rowToCells_[i]) {
      if (cellToRow_[c] != i) {
        throw std::runtime_error(
            "Cell allocation does not match row allocation");
      }
      if (cellToX_[c] < rows_[i].minX) {
        throw std::runtime_error("Cell placed before the row");
      }
      if (cellToX_[c] + cellWidth_[c] > rows_[i].maxX) {
        throw std::runtime_error("Cell placed after the row");
      }
    }
  }
  for (int i = 0; i < nbRows(); ++i) {
    for (int j = 0; j + 1 < rowToCells_[i].size(); ++j) {
      int c1 = rowToCells_[i][j];
      int c2 = rowToCells_[i][j + 1];
      if (cellToX_[c1] + cellWidth_[c1] > cellToX_[c2]) {
        throw std::runtime_error("Cell overlap detected");
      }
    }
  }
}

bool Legalizer::placeCellOptimally(int cell) {
  /**
   * Very naive algorithm that tries every possible position
   */
  if (isIgnored(cell)) {
    return true;
  }
  int targetX = cellTargetX_[cell];
  int targetY = cellTargetY_[cell];
  int bestX = 0;
  int bestRow = -1;
  long long bestDist = std::numeric_limits<long long>::max();

  auto tryPlace = [&](int row) {
    int y = rows_[row].minY;
    if (bestRow != -1 && norm(0, y - targetY, costModel_) > bestDist) {
      // Not possible to do better since the rows are sorted
      return true;
    }
    // Find the best position fo the cell
    auto [ok, x] = placeCellOptimally(cell, row);
    if (!ok) {
      // Not possible to place in this row, but cannot stop yet
      return false;
    }
    auto dist = norm(x - targetX, y - targetY, costModel_);
    if (bestRow == -1 || dist < bestDist) {
      bestX = x;
      bestRow = row;
      bestDist = dist;
    }
    // Cannot stop yet
    return false;
  };

  // Try promising candidates first
  int initialRow = closestRow(targetY);
  for (int row = initialRow; row < nbRows(); ++row) {
    bool canStop = tryPlace(row);
    if (canStop) {
      break;
    }
  }
  for (int row = initialRow - 1; row >= 0; --row) {
    bool canStop = tryPlace(row);
    if (canStop) {
      break;
    }
  }

  if (bestRow == -1) {
    report();
    check();
    throw std::runtime_error(
        "Unable to place a cell with the naive legalization algorithm");
  }
  doPlacement(cell, bestRow, bestX);
  return true;
}

std::pair<bool, int> Legalizer::placeCellOptimally(int cell, int row) const {
  int minCoord = firstFreeX(row);
  int maxCoord = rows_[row].maxX - cellWidth_[cell];

  if (minCoord > maxCoord) {
    return std::make_pair(false, 0);
  }

  // Get the best possible placement
  // OPTIMIZE: with LInf cost we could pack further left without degradation
  int coord = cellTargetX_[cell];
  coord = std::min(coord, maxCoord);
  coord = std::max(coord, minCoord);
  return std::make_pair(true, coord);
}

std::vector<int> Legalizer::computeCellOrder(float weightX, float weightWidth,
                                             float weightY) const {
  // Sort the cells by target X coordinate
  std::vector<std::pair<float, int> > sortedCells;
  for (int i = 0; i < nbCells(); ++i) {
    float val = weightX * cellTargetX_[i] + weightWidth * cellWidth_[i] +
                weightY * cellTargetY_[i];
    sortedCells.emplace_back(val, i);
  }
  std::stable_sort(sortedCells.begin(), sortedCells.end());
  std::vector<int> cells;
  for (auto p : sortedCells) {
    cells.push_back(p.second);
  }
  return cells;
}

void Legalizer::doPlacement(int cell, int row, int x) {
  assert(row >= 0 && row < nbRows());
  rowToCells_[row].push_back(cell);
  cellToRow_[cell] = row;
  cellToX_[cell] = x;
  cellToY_[cell] = rows_[row].minY;
}

void Legalizer::undoPlacement(int cell) {
  int row = cellToRow_[cell];
  if (row == -1) return;  // Not placed

  // Remove from the row
  auto pos = std::find(rowToCells_[row].begin(), rowToCells_[row].end(), cell);
  assert(pos != rowToCells_[row].end());
  rowToCells_[row].erase(pos);
  cellToRow_[cell] = -1;

  // Back to default
  cellToX_[cell] = cellTargetX_[cell];
  cellToY_[cell] = cellTargetY_[cell];
}

int Legalizer::firstFreeX(int row) const {
  assert(row >= 0 && row < nbRows());
  if (rowToCells_[row].empty()) {
    return rows_[row].minX;
  }
  int cell = rowToCells_[row].back();
  return cellToX_[cell] + cellWidth_[cell];
}

int Legalizer::closestRow(int y) const {
  auto it = std::lower_bound(rows_.begin(), rows_.end(), y,
                             [](Rectangle r, int v) { return r.minY < v; });
  if (it == rows_.end()) return nbRows() - 1;
  if (it == rows_.begin()) return 0;
  int row = it - rows_.begin();
  assert(row >= 1);
  if (rows_[row].minX - y > y - rows_[row - 1].minX) {
    return row - 1;
  } else {
    return row;
  }
}

std::vector<int> Legalizer::cellLegalX() const {
  std::vector<int> ret(nbCells());
  for (int r = 0; r < nbRows(); ++r) {
    for (int i = 0; i < rowToCells_[r].size(); ++i) {
      int c = rowToCells_[r][i];
      ret[c] = cellToX_[c];
    }
  }
  return ret;
}

std::vector<int> Legalizer::cellLegalY() const {
  std::vector<int> ret(nbCells());
  for (int r = 0; r < nbRows(); ++r) {
    for (int i = 0; i < rowToCells_[r].size(); ++i) {
      int c = rowToCells_[r][i];
      ret[c] = rows_[r].minY;
    }
  }
  return ret;
}

void Legalizer::exportPlacement(Circuit &circuit) {
  std::vector<int> cellX = cellLegalX();
  std::vector<int> cellY = cellLegalY();
  for (int i = 0; i < circuit.nbCells(); ++i) {
    circuit.cellX[i] = cellX[i];
    circuit.cellY[i] = cellY[i];
  }
}

void Legalizer::report(bool verbose) const {
  std::cout << "Legalizer with " << nbCells() << " cells on " << nbRows()
            << " rows" << std::endl;
  if (!verbose) return;
  for (int row = 0; row < nbRows(); ++row) {
    std::cout << "Row " << rows_[row].minY << ", " << rows_[row].minX << " to "
              << rows_[row].maxX << ": ";
    for (int c : rowToCells_[row]) {
      int x = cellToX_[c];
      int w = cellWidth_[c];
      std::cout << c << " (" << x << " - " << x + w << ") ";
    }
    std::cout << std::endl;
  }
}