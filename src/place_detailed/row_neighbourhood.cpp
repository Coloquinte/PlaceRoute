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
}  // namespace

void RowNeighbourhood::simpleSetup(const std::vector<Rectangle> &rows,
                                   int nbNeighbourRows) {
  rowsBelow_.resize(rows.size());
  rowsAbove_.resize(rows.size());
  rowsLeft_.resize(rows.size());
  rowsRight_.resize(rows.size());
  for (int row = 0; row < rows.size(); ++row) {
    rowsBelow_[row] = keepFirstK(rowsBelow(rows[row], rows), nbNeighbourRows);
    rowsAbove_[row] = keepFirstK(rowsAbove(rows[row], rows), nbNeighbourRows);
    rowsLeft_[row] = keepFirstK(rowsLeft(rows[row], rows), nbNeighbourRows);
    rowsRight_[row] = keepFirstK(rowsRight(rows[row], rows), nbNeighbourRows);
  }
}

std::vector<int> RowNeighbourhood::rowsBelow(
    Rectangle row, const std::vector<Rectangle> &rows) {
  std::vector<std::pair<int, Rectangle> > sortedRows;
  for (int other = 0; other < rows.size(); ++other) {
    // Only rows strictly below
    if (row.minY <= rows[other].minY) continue;
    // Only rows that share part of their x range
    if (row.minX >= rows[other].maxX) continue;
    if (row.maxX <= rows[other].minX) continue;
    sortedRows.emplace_back(other, rows[other]);
  }
  std::sort(sortedRows.begin(), sortedRows.end(),
            [](std::pair<int, Rectangle> a, std::pair<int, Rectangle> b) {
              int ay = a.second.minY;
              int by = b.second.minY;
              return ay > by || (ay == by && a.first < b.first);
            });
  std::vector<int> ret;
  for (auto p : sortedRows) {
    ret.push_back(p.first);
  }
  return ret;
}

std::vector<int> RowNeighbourhood::rowsAbove(
    Rectangle row, const std::vector<Rectangle> &rows) {
  std::vector<std::pair<int, Rectangle> > sortedRows;
  for (int other = 0; other < rows.size(); ++other) {
    // Only rows strictly above
    if (row.minY >= rows[other].minY) continue;
    // Only rows that share part of their x range
    if (row.minX >= rows[other].maxX) continue;
    if (row.maxX <= rows[other].minX) continue;
    sortedRows.emplace_back(other, rows[other]);
  }
  std::sort(sortedRows.begin(), sortedRows.end(),
            [](std::pair<int, Rectangle> a, std::pair<int, Rectangle> b) {
              int ay = a.second.minY;
              int by = b.second.minY;
              return ay < by || (ay == by && a.first < b.first);
            });
  std::vector<int> ret;
  for (auto p : sortedRows) {
    ret.push_back(p.first);
  }
  return ret;
}

std::vector<int> RowNeighbourhood::rowsLeft(
    Rectangle row, const std::vector<Rectangle> &rows) {
  std::vector<std::pair<int, Rectangle> > sortedRows;
  for (int other = 0; other < rows.size(); ++other) {
    // Only rows strictly on the left
    if (row.minX < rows[other].maxX) continue;
    sortedRows.emplace_back(other, rows[other]);
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
    // Only rows strictly on the right
    if (row.maxX > rows[other].minX) continue;
    sortedRows.emplace_back(other, rows[other]);
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
