#include "place_detailed/row_legalizer.hpp"

inline int OSRP_leg::get_displacement(legalizable_task const newly_pushed,
                                      bool update) {
  int target_abs_pos = newly_pushed.target_pos - current_width();
  int width = newly_pushed.width;
  int slope = -width;

  int cur_pos = end;
  int cur_cost = 0;

  std::vector<OSRP_bound> passed_bounds;

  while (not bounds.empty() and
         ((slope < 0 and bounds.top().absolute_pos >
                             target_abs_pos)  // Not reached equilibrium
          or bounds.top().absolute_pos >
                 end - current_width() - width)  // Still not a legal position
  ) {
    int old_pos = cur_pos;
    cur_pos = bounds.top().absolute_pos;
    cur_cost +=
        (old_pos - cur_pos) *
        (slope + width);  // The additional cost for the other cells encountered
    slope += bounds.top().weight;

    // Remember which bounds we encountered in order to reset the object to its
    // initial state
    if (not update) passed_bounds.push_back(bounds.top());
    bounds.pop();
  }

  int final_abs_pos = std::min(
      end - current_width() -
          width,  // Always before the end and after the beginning
      std::max(begin, slope >= 0
                          ? cur_pos
                          : target_abs_pos)  // but did we stop before reaching
                                             // the target position?
  );

  cur_cost +=
      (cur_pos - final_abs_pos) *
      (slope + width);  // The additional cost for the other cells encountered

  assert(final_abs_pos >= begin);
  assert(final_abs_pos <= end - current_width() - width);

  if (update) {
    prev_width.push_back(width + current_width());
    cells.push_back(newly_pushed.ind);
    constraining_pos.push_back(final_abs_pos);
    if (slope > 0) {  // Remaining capacity of an encountered bound
      bounds.push(OSRP_bound(slope, cur_pos));
    }
    // The new bound, minus what it absorbs of the remaining slope
    if (target_abs_pos > begin) {
      bounds.push(OSRP_bound(2 * width + std::min(slope, 0), target_abs_pos));
    }
  } else {
    for (OSRP_bound b : passed_bounds) {
      bounds.push(b);
    }
  }

  return cur_cost +
         width * std::abs(final_abs_pos -
                          target_abs_pos);  // Add the cost of the new cell
}

inline std::vector<std::pair<int, int> > OSRP_leg::get_placement() const {
  auto final_abs_pos = constraining_pos;
  std::partial_sum(final_abs_pos.rbegin(), final_abs_pos.rend(),
                   final_abs_pos.rbegin(),
                   [](int a, int b) -> T { return std::min(a, b); });

  std::vector<std::pair<int, int> > ret(cells.size());
  for (int i = 0; i < cells.size(); ++i) {
    ret[i] = std::make_pair(cells[i], final_abs_pos[i] + prev_width[i]);

    assert(final_abs_pos[i] >= begin);
    assert(final_abs_pos[i] + prev_width[i + 1] <= end);
  }
  return ret;
}