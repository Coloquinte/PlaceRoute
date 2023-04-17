
#include "place_detailed/tetris_legalizer.hpp"

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

void TetrisLegalizer::placeCell(int cell) {
  int targetX = cellTargetX_[cell];
  int targetY = cellTargetY_[cell];
  int bestX = 0;
  int bestRow = -1;

  auto tryPlace = [](int row) {

  };
}

std::vector<std::pair<int, int> > TetrisLegalizer::getPossibleIntervals(
    int w, int h, int y) const {
  int rowInd = closestRow(y);
  std::vector<std::pair<int, int> > ret;
  ret.emplace_back(rows_[rowInd].minX, rows_[rowInd].maxX);
  return ret;
  // TODO
}
}  // namespace coloquinte
