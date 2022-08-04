#pragma once

#include <vector>
#include <queue>


struct legalizable_task {
  int width;
  int target_pos;
  int ind;
  legalizable_task(int w, int p, int i) : width(w), target_pos(p), ind(i) {}
  bool operator<(legalizable_task const o) const {
    return target_pos < o.target_pos;
  }
};

/**
 * @brief Obtain the positions minimizing total weighted displacement along a
 *row.
 *
 * It is an ordered single row problem/fixed order single machine scheduling
 *problem, solved by the clumping/specialized cascading descent algorithm.
 *
 * The cost model is linear in the distance to the target position, weighted by
 *the width of the cells
 **/
class OSRP_leg {
  struct OSRP_bound {
    int absolute_pos;  // Will be the target absolute position of the cell
    int weight;        // Will be the width of the cell

    bool operator<(OSRP_bound const o) const {
      return absolute_pos < o.absolute_pos;
    }
    OSRP_bound(int w, int abs_pos) : absolute_pos(abs_pos), weight(w) {}
  };

  int begin, end;

  std::vector<int> cells;  // The indexes in the circuit
  std::vector<int>
      constraining_pos;  // Where the cells have been pushed and constrain the
                         // positions of preceding cells
  std::vector<int> prev_width;  // Cumulative width of the cells: calculates the
                                // absolute position of new cells

  std::priority_queue<OSRP_bound> bounds;

  // Get the cost of pushing a cell on the row
  int get_displacement(legalizable_task const newly_pushed, bool update);

 public:
  int current_width() const { return prev_width.back(); }
  int remaining_space() const { return end - begin - current_width(); }
  int last_available_pos() const {
    return constraining_pos.back() + current_width();
  }

  int get_cost(legalizable_task const task) {
    return get_displacement(task, false);
  }
  void push(legalizable_task const task) { get_displacement(task, true); }

  // Initialize
  OSRP_leg(int b, int e) : begin(b), end(e), prev_width(1, 0) {}
  OSRP_leg() {}

  // Get the resulting placement
  std::vector<std::pair<int, int> > get_placement() const;
};
