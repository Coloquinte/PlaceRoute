#include "place_detailed/row_neighbourhood.hpp"

#include <algorithm>
#include <stdexcept>

namespace coloquinte {
RowNeighbourhood::RowNeighbourhood(const std::vector<Rectangle> &rows,
                                   int nbNeighbourRows) {
  simpleSetup(rows, nbNeighbourRows);
}

RowNeighbourhood::RowNeighbourhood(const std::vector<Row> &rows,
                                   int nbNeighbourRows) {
  std::vector<Rectangle> r;
  r.reserve(rows.size());
  for (Rectangle row : rows) {
    r.push_back(row);
  }
  simpleSetup(r, nbNeighbourRows);
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
  if (inds.size() <= nb) {
    return inds;
  }
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

bool orderSide(std::pair<int, Rectangle> a, std::pair<int, Rectangle> b) {
  int ax = a.second.minX;
  int bx = b.second.minX;
  int ay = a.second.minY;
  int by = b.second.minY;
  return ay < by || (ay == by && ax < bx);
}
}  // namespace

void RowNeighbourhood::simpleSetup(const std::vector<Rectangle> &rows,
                                   int nbNeighbourRows) {
  rowsBelow_ = rowsBelow(rows, nbNeighbourRows);
  rowsAbove_ = rowsAbove(rows, nbNeighbourRows);
  buildRowsSides(rows, nbNeighbourRows);
}

bool RowNeighbourhood::isBelow(Rectangle r1, Rectangle r2) {
  // Only rows strictly below
  if (r2.minY <= r1.minY) {
    return false;
  }
  // Only rows that share part of their x range
  if (r2.minX >= r1.maxX) {
    return false;
  }
  if (r2.maxX <= r1.minX) {
    return false;
  }
  return true;
}

bool RowNeighbourhood::isAbove(Rectangle r1, Rectangle r2) {
  // Only rows strictly above
  if (r2.minY >= r1.minY) {
    return false;
  }
  // Only rows that share part of their x range
  if (r2.minX >= r1.maxX) {
    return false;
  }
  if (r2.maxX <= r1.minX) {
    return false;
  }
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
  sortedRows.reserve(rows.size());

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
  sortedRows.reserve(rows.size());

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

std::vector<int> RowNeighbourhood::buildLeftFrom(
    Rectangle row, const std::vector<Rectangle> &rows, int ind) const {
  std::vector<int> candidates;
  candidates.push_back(ind);
  for (int c : rowsAbove(ind)) {
    candidates.push_back(c);
  }
  for (int c : rowsBelow(ind)) {
    candidates.push_back(c);
  }
  std::vector<std::pair<int, Rectangle> > sortedRows;
  for (int c : candidates) {
    if (isLeft(rows[c], row)) {
      sortedRows.emplace_back(c, rows[c]);
    }
  }
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
  ret.reserve(sortedRows.size());

  for (auto p : sortedRows) {
    ret.push_back(p.first);
  }
  return ret;
}

std::vector<int> RowNeighbourhood::buildRightFrom(
    Rectangle row, const std::vector<Rectangle> &rows, int ind) const {
  std::vector<int> candidates;
  candidates.push_back(ind);
  for (int c : rowsAbove(ind)) {
    candidates.push_back(c);
  }
  for (int c : rowsBelow(ind)) {
    candidates.push_back(c);
  }
  std::vector<std::pair<int, Rectangle> > sortedRows;
  for (int c : candidates) {
    if (isRight(rows[c], row)) {
      sortedRows.emplace_back(c, rows[c]);
    }
  }
  std::sort(sortedRows.begin(), sortedRows.end(),
            [row](std::pair<int, Rectangle> a, std::pair<int, Rectangle> b) {
              int ax = a.second.maxX;
              int bx = b.second.maxX;
              int ay = a.second.minY;
              int by = b.second.minY;
              int rx = row.maxX;
              int ry = row.minY;
              int ad = std::abs(ax - rx) + std::abs(ay - ry);
              int bd = std::abs(bx - rx) + std::abs(by - ry);
              return ad < bd || (ad == bd && a.first < b.first);
            });
  std::vector<int> ret;
  ret.reserve(sortedRows.size());

  for (auto p : sortedRows) {
    ret.push_back(p.first);
  }
  return ret;
}

void RowNeighbourhood::buildRowsSides(const std::vector<Rectangle> &rows,
                                      int nbNeighbourRows) {
  rowsLeft_.clear();
  rowsLeft_.resize(rows.size());
  rowsRight_.clear();
  rowsRight_.resize(rows.size());
  std::vector<std::pair<int, Rectangle> > sortedRows;
  sortedRows.reserve(rows.size());

  for (int i = 0; i < rows.size(); ++i) {
    sortedRows.emplace_back(i, rows[i]);
  }
  std::sort(sortedRows.begin(), sortedRows.end(), orderSide);
  for (int i = 0; i + 1 < sortedRows.size(); ++i) {
    auto [ind1, row1] = sortedRows[i];
    auto [ind2, row2] = sortedRows[i + 1];
    if (isLeft(row1, row2)) {
      rowsLeft_[ind2] =
          keepFirstK(buildLeftFrom(row2, rows, ind1), nbNeighbourRows);
    }
    if (isRight(row2, row1)) {
      rowsRight_[ind1] =
          keepFirstK(buildRightFrom(row1, rows, ind2), nbNeighbourRows);
    }
  }
}

}  // Namespace coloquinte
