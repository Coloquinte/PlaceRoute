#include "place_detailed/row_legalizer.hpp"

#include <cassert>
#include <cmath>
#include <numeric>
#include <utility>

namespace coloquinte {
inline int RowLegalizer::getDisplacement(int width, int targetPos,
                                         bool update) {
  int targetAbsPos = targetPos - usedSpace();
  int slope = -width;

  int cur_pos = end_;
  int cur_cost = 0;

  std::vector<Bound> passed_bounds;

  // While the slope is negative or the position is not legal yet
  while (not bounds.empty() and
         ((slope < 0 and bounds.top().absolutePos > targetAbsPos) or
          bounds.top().absolutePos > end_ - usedSpace() - width)) {
    int old_pos = cur_pos;
    cur_pos = bounds.top().absolutePos;
    cur_cost += (old_pos - cur_pos) * (slope + width);
    slope += bounds.top().weight;

    // Remember which bounds we encountered in order to reset the object to its
    // initial state
    if (not update) passed_bounds.push_back(bounds.top());
    bounds.pop();
  }
  // Always before the end, after the beginning, and not pushed beyond the
  // target position
  int finalAbsPos =
      std::min(end_ - usedSpace() - width,
               std::max(begin_, slope >= 0 ? cur_pos : targetAbsPos));

  cur_cost += (cur_pos - finalAbsPos) * (slope + width);

  assert(finalAbsPos >= begin_);
  assert(finalAbsPos <= end_ - usedSpace() - width);

  if (update) {
    cumWidth_.push_back(width + usedSpace());
    constrainingPos_.push_back(finalAbsPos);
    if (slope > 0) {  // Remaining capacity of an encountered bound
      bounds.push(Bound(slope, cur_pos));
    }
    // The new bound, minus what it absorbs of the remaining slope
    if (targetAbsPos > begin_) {
      bounds.push(Bound(2 * width + std::min(slope, 0), targetAbsPos));
    }
  } else {
    for (Bound b : passed_bounds) {
      bounds.push(b);
    }
  }

  return cur_cost +
         width * std::abs(finalAbsPos -
                          targetAbsPos);  // Add the cost of the new cell
}

int RowLegalizer::getCost(int width, int targetPos) {
  return getDisplacement(width, targetPos, false);
}

int RowLegalizer::push(int width, int targetPos) {
  return getDisplacement(width, targetPos, true);
}

void RowLegalizer::clear() {
  cumWidth_.assign(1, 0);
  bounds = std::priority_queue<Bound>();
  constrainingPos_.clear();
}

std::vector<int> RowLegalizer::getPlacement() const {
  auto finalAbsPos = constrainingPos_;
  std::partial_sum(finalAbsPos.rbegin(), finalAbsPos.rend(),
                   finalAbsPos.rbegin(),
                   [](int a, int b) -> int { return std::min(a, b); });

  std::vector<int> ret(finalAbsPos.size());
  for (int i = 0; i < finalAbsPos.size(); ++i) {
    ret[i] = finalAbsPos[i] + cumWidth_[i];

    assert(finalAbsPos[i] >= begin_);
    assert(finalAbsPos[i] + cumWidth_[i + 1] <= end_);
  }
  return ret;
}
}