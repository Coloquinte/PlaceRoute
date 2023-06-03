

#include "place_detailed/abacus_legalizer.hpp"

#include "utils/norm.hpp"

namespace coloquinte {

AbacusLegalizer::AbacusLegalizer(
    const std::vector<Row> &rows, const std::vector<int> &width,
    const std::vector<int> &height,
    const std::vector<CellRowPolarity> &polarities,
    const std::vector<int> &targetX, const std::vector<int> &targetY,
    const std::vector<CellOrientation> &targetOrientation)
    : LegalizerBase(rows, width, height, polarities, targetX, targetY,
                    targetOrientation) {
  rowToCells_.resize(rows_.size());
  for (const Row &row : rows_) {
    rowLegalizers_.emplace_back(row.minX, row.maxX);
  }
}

void AbacusLegalizer::run() {
  for (int c = 0; c < nbCells(); ++c) {
    placeCell(c);
  }
  for (int i = 0; i < nbRows(); ++i) {
    std::vector<int> pl = rowLegalizers_[i].getPlacement();
    assert(pl.size() == rowToCells_[i].size());
    for (int j = 0; j < pl.size(); ++j) {
      int cell = rowToCells_[i][j];
      cellToX_[cell] = pl[j];
      cellToY_[cell] = rows_[i].minY;
      cellToOrientation_[cell] = getOrientation(cell, i);
    }
  }
  check();
}

std::pair<bool, long long> AbacusLegalizer::evaluatePlacement(int cell,
                                                              int row) {
  if (rowLegalizers_[row].remainingSpace() < cellWidth_[cell]) {
    return std::make_pair(false, 0);
  }
  if (getOrientation(cell, row) == CellOrientation::INVALID) {
    return std::make_pair(false, 0);
  }
  int dist = rowLegalizers_[row].getCost(cellWidth_[cell], cellTargetX_[cell]);
  return std::make_pair(true, dist);
}

void AbacusLegalizer::placeCell(int cell) {
  /**
   * Simple algorithm that tries close row first and stops early if no
   * improvement can be found
   */
  int targetX = cellTargetX_[cell];
  int targetY = cellTargetY_[cell];
  int bestX = 0;
  int bestRow = -1;
  long long bestDist = std::numeric_limits<long long>::max();

  auto tryPlace = [&](int row) {
    if (rows_[row].height() != cellHeight_[cell]) {
      // Forbidden to place a row in a different-height cell
      return false;
    }
    long long yDist = cellWidth_[cell] *
                      norm(0, rows_[row].minY - targetY, LegalizationModel::L1);
    if (bestRow != -1 && yDist > bestDist) {
      // Not possible to do better since the rows are sorted
      return true;
    }
    // Find the best position for the cell
    auto [ok, xDist] = evaluatePlacement(cell, row);
    // TODO: extend this to non-L1 cases
    long long dist = xDist + yDist;
    if (!ok) {
      // Not possible to place in this row, but cannot stop yet
      return false;
    }
    if (bestRow == -1 || dist < bestDist) {
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
    return;
  }
  rowLegalizers_[bestRow].push(cellWidth_[cell], targetX);
  rowToCells_[bestRow].push_back(cell);
  cellIsPlaced_[cell] = true;
}

void AbacusLegalizer::check() const {
  LegalizerBase::check();
  for (int i = 0; i < nbRows(); ++i) {
    for (int c : rowToCells_[i]) {
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

}  // namespace coloquinte
