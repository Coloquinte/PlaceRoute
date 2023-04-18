
#include "place_detailed/tetris_legalizer.hpp"

#include <limits>

namespace coloquinte {
TetrisLegalizer::TetrisLegalizer(const std::vector<Row> &rows,
                                 const std::vector<int> &width,
                                 const std::vector<int> &height,
                                 const std::vector<CellRowPolarity> &polarities,
                                 const std::vector<int> &targetX,
                                 const std::vector<int> &targetY)
    : LegalizerBase(rows, width, height, polarities, targetX, targetY) {
  for (Row r : rows_) {
    rowFreePos_.push_back(r.minX);
  }
}

void TetrisLegalizer::run() {
  for (int c = 0; c < nbCells(); ++c) {
    placeCell(c);
  }
}

std::pair<bool, int> TetrisLegalizer::attemptPlacement(int cell, int y) const {
  auto p = getPossibleIntervals(cellWidth_[cell], cellHeight_[cell], y);
  if (p.empty()) {
    return std::make_pair(false, 0);
  }
  int x = cellTargetX_[cell];
  int dest = 0;
  bool found = false;
  for (auto [b, e] : p) {
    int pos = std::clamp(x, b, e);
    if (!found || std::abs(pos - x) < std::abs(dest - x)) {
      dest = pos;
    }
  }
  return std::make_pair(true, dest);
}

void TetrisLegalizer::placeCell(int cell) {
  int targetX = cellTargetX_[cell];
  int targetY = cellTargetY_[cell];
  int bestX = 0;
  int bestY = 0;
  int bestDist = std::numeric_limits<int>::max();
  bool found = false;

  auto tryPlace = [&](int y) -> bool {
    if (found && std::abs(targetY - y) >= bestDist) {
      // No need to continue searching as we cannot find better solution
      return true;
    }
    auto [ok, x] = attemptPlacement(cell, y);
    if (!ok) {
      // Not possible to place in this row
      return false;
    }
    int dist = std::abs(targetX - x) + std::abs(targetY - y);
    if (!found || dist < bestDist) {
      found = true;
      bestDist = dist;
      bestX = x;
      bestY = y;
    }
    return false;
  };

  // Try promising candidates first
  int initialRow = closestRow(targetY);
  for (int row = initialRow; row < nbRows(); ++row) {
    bool canStop = tryPlace(rows_[row].minY);
    if (canStop) {
      break;
    }
  }
  for (int row = initialRow - 1; row >= 0; --row) {
    bool canStop = tryPlace(rows_[row].minY);
    if (canStop) {
      break;
    }
  }

  if (!found) {
    return;
  }
  cellToX_[cell] = bestX;
  cellToY_[cell] = bestY;
  cellIsPlaced_[cell] = true;
  instanciateCell(bestX, bestY, cellWidth_[cell], cellHeight_[cell]);
}

void TetrisLegalizer::instanciateCell(int x, int y, int w, int h) {
  if (h <= 0 || w <= 0) {
    return;
  }
  int rowInd = closestRow(y);
  for (int r = rowInd; r < nbRows(); ++r) {
    if (rows_[r].minY != y) {
      break;
    }
    if (x < rows_[r].maxX && x + w > rows_[r].minX) {
      rowFreePos_[r] = x + w;
    }
  }
  if (h <= rowHeight()) {
    return;
  }
  instanciateCell(x, y + rowHeight(), w, h - rowHeight());
}

std::vector<std::pair<int, int> > TetrisLegalizer::getPossibleIntervals(
    int w, int h, int y) const {
  int rowInd = closestRow(y);
  std::vector<std::pair<int, int> > intervals;
  for (int r = rowInd; r < nbRows(); ++r) {
    if (rows_[r].minY != y) {
      break;
    }
    int b = rowFreePos_[r];
    int e = rows_[r].maxX - w;
    if (e >= b) {
      intervals.emplace_back(b, e);
    }
  }
  if (h <= rowHeight() || intervals.empty()) {
    return intervals;
  }
  std::vector<std::pair<int, int> > other =
      getPossibleIntervals(w, h - rowHeight(), y + rowHeight());

  std::vector<std::pair<int, int> > ret;
  for (auto [b1, e1] : intervals) {
    for (auto [b2, e2] : other) {
      if (b1 <= e2 && b2 <= e1) {
        ret.emplace_back(std::max(b1, b2), std::min(e1, e2));
      }
    }
  }
  return ret;
}
}  // namespace coloquinte
