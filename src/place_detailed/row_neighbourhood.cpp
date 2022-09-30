#include "place_detailed/row_neighbourhood.hpp"

#include <algorithm>
#include <stdexcept>

namespace coloquinte {
RowNeighbourhood::RowNeighbourhood(const std::vector<Rectangle> &rows,
                                   int nbNeighbourRows) {
  simpleSetup(rows, nbNeighbourRows);
}

void RowNeighbourhood::check() const {
  if (rowsBelow_.size() != nbRows()) {
    throw std::runtime_error("Inconsistent number of rows in neighbourhood");
  }
  if (rowsAbove_.size() != nbRows()) {
    throw std::runtime_error("Inconsistent number of rows in neighbourhood");
  }
  if (rowsLeft_.size() != nbRows()) {
    throw std::runtime_error("Inconsistent number of rows in neighbourhood");
  }
  if (rowsRight_.size() != nbRows()) {
    throw std::runtime_error("Inconsistent number of rows in neighbourhood");
  }
}

namespace {
std::vector<int> keepFirstK(const std::vector<int> &inds, int nb) {
  if (inds.size() <= nb) return inds;
  return std::vector<int>(inds.begin(), inds.begin() + nb);
}

bool orderBelow(std::pair<int, Rectangle> a, std::pair<int, Rectangle> b) {
  int ay = a.second.minY;
  int by = b.second.minY;
  return ay > by || (ay == by && a.first < b.first);
}

bool orderAbove(std::pair<int, Rectangle> a, std::pair<int, Rectangle> b) {
  int ay = a.second.minY;
  int by = b.second.minY;
  return ay < by || (ay == by && a.first < b.first);
}
}  // namespace

void RowNeighbourhood::simpleSetup(const std::vector<Rectangle> &rows,
                                   int nbNeighbourRows) {
  rowsBelow_ = rowsBelow(rows, nbNeighbourRows);
  rowsAbove_ = rowsAbove(rows, nbNeighbourRows);
  rowsLeft_.resize(rows.size());
  rowsRight_.resize(rows.size());

  for (int row = 0; row < rows.size(); ++row) {
    rowsLeft_[row] = keepFirstK(rowsLeft(rows[row], rows), nbNeighbourRows);
    rowsRight_[row] = keepFirstK(rowsRight(rows[row], rows), nbNeighbourRows);
  }
}

bool RowNeighbourhood::isBelow(Rectangle r1, Rectangle r2) {
  // Only rows strictly below
  if (r2.minY <= r1.minY) return false;
  // Only rows that share part of their x range
  if (r2.minX >= r1.maxX) return false;
  if (r2.maxX <= r1.minX) return false;
  return true;
}

bool RowNeighbourhood::isAbove(Rectangle r1, Rectangle r2) {
  // Only rows strictly above
  if (r2.minY >= r1.minY) return false;
  // Only rows that share part of their x range
  if (r2.minX >= r1.maxX) return false;
  if (r2.maxX <= r1.minX) return false;
  return true;
}

bool RowNeighbourhood::isLeft(Rectangle r1, Rectangle r2) {
  // Only rows strictly on the left
  return r2.minX >= r1.maxX;
}

bool RowNeighbourhood::isRight(Rectangle r1, Rectangle r2) {
  // Only rows strictly on the right
  return r2.maxX <= r1.minX;
}

std::vector<std::vector<int> > RowNeighbourhood::rowsBelow(
    const std::vector<Rectangle> &rows, int nbNeighbourRows) {
  std::vector<std::vector<int> > ret(rows.size());
  std::vector<std::pair<int, Rectangle> > sortedRows;
  for (int i = 0; i < rows.size(); ++i) {
    sortedRows.emplace_back(i, rows[i]);
  }
  std::sort(sortedRows.begin(), sortedRows.end(), orderBelow);
  for (int i = 0; i < sortedRows.size(); ++i) {
    auto [ind1, row1] = sortedRows[i];
    int nbFound = 0;
    for (int j = i + 1; j < sortedRows.size(); ++j) {
      auto [ind2, row2] = sortedRows[j];
      if (isBelow(row2, row1)) {
        ret[ind1].push_back(ind2);
        ++nbFound;
      }
      if (nbFound >= nbNeighbourRows) {
        break;
      }
    }
  }
  return ret;
}

std::vector<std::vector<int> > RowNeighbourhood::rowsAbove(
    const std::vector<Rectangle> &rows, int nbNeighbourRows) {
  std::vector<std::vector<int> > ret(rows.size());
  std::vector<std::pair<int, Rectangle> > sortedRows;
  for (int i = 0; i < rows.size(); ++i) {
    sortedRows.emplace_back(i, rows[i]);
  }
  std::sort(sortedRows.begin(), sortedRows.end(), orderAbove);
  for (int i = 0; i < sortedRows.size(); ++i) {
    auto [ind1, row1] = sortedRows[i];
    int nbFound = 0;
    for (int j = i + 1; j < sortedRows.size(); ++j) {
      auto [ind2, row2] = sortedRows[j];
      if (isAbove(row2, row1)) {
        ret[ind1].push_back(ind2);
        ++nbFound;
      }
      if (nbFound >= nbNeighbourRows) {
        break;
      }
    }
  }
  return ret;
}

std::vector<int> RowNeighbourhood::rowsLeft(
    Rectangle row, const std::vector<Rectangle> &rows) {
  std::vector<std::pair<int, Rectangle> > sortedRows;
  for (int other = 0; other < rows.size(); ++other) {
    if (isLeft(rows[other], row)) sortedRows.emplace_back(other, rows[other]);
  }
  // Sort by distance to the row's left side
  std::sort(sortedRows.begin(), sortedRows.end(),
            [row](std::pair<int, Rectangle> a, std::pair<int, Rectangle> b) {
              int ax = a.second.maxX;
              int bx = b.second.maxX;
              int ay = a.second.minY;
              int by = b.second.minY;
              int rx = row.minX;
              int ry = row.minY;
              int ad = std::abs(ax - rx) + std::abs(ay - ry);
              int bd = std::abs(bx - rx) + std::abs(by - ry);
              return ad < bd || (ad == bd && a.first < b.first);
            });
  std::vector<int> ret;
  for (auto p : sortedRows) {
    ret.push_back(p.first);
  }
  return ret;
}

std::vector<int> RowNeighbourhood::rowsRight(
    Rectangle row, const std::vector<Rectangle> &rows) {
  std::vector<std::pair<int, Rectangle> > sortedRows;
  for (int other = 0; other < rows.size(); ++other) {
    if (isRight(rows[other], row)) sortedRows.emplace_back(other, rows[other]);
  }
  // Sort by distance to the row's right side
  std::sort(sortedRows.begin(), sortedRows.end(),
            [row](std::pair<int, Rectangle> a, std::pair<int, Rectangle> b) {
              int ax = a.second.minX;
              int bx = b.second.minX;
              int ay = a.second.minY;
              int by = b.second.minY;
              int rx = row.maxX;
              int ry = row.minY;
              int ad = std::abs(ax - rx) + std::abs(ay - ry);
              int bd = std::abs(bx - rx) + std::abs(by - ry);
              return ad < bd || (ad == bd && a.first < b.first);
            });
  std::vector<int> ret;
  for (auto p : sortedRows) {
    ret.push_back(p.first);
  }
  return ret;
}
}  // Namespace coloquinte
