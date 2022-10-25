
#include "place_global/transportation.hpp"

#include <cassert>
#include <limits>
#include <numeric>
#include <queue>
#include <stdexcept>
#include <utility>
#include <vector>

namespace coloquinte {

class CurrentAllocation {
  static constexpr int NULL_IND = -1;

  // Internal data structures

  // Priority queue element to determine the source to be used between regions
  struct MovableSource {
    int source;
    float cost;
    bool operator<(MovableSource const o) const {
      return cost > o.cost  // Sorted by cost
             || (cost == o.cost &&
                 source < o.source);  // And by index to limit the number of
                                      // fractional elements between two regions
    }
    MovableSource(int s, float c) : source(s), cost(c) {}
  };

  // Member data

  // The current state
  // For each region, for each source, the capacity allocated by the region
  std::vector<std::vector<long long> > sr_allocations_;
  // The costs from a region to a source
  std::vector<std::vector<float> > sr_costs_;
  // The demands of the sources
  std::vector<long long> s_demands_;
  // The remaining capacities of the regions
  std::vector<long long> r_capacities_;

  // Shortest path data

  // The costs of allocating to a region
  std::vector<float> r_costs_;
  // The parents of the regions i.e. the regions where we push sources first (or
  // null_ind)
  std::vector<int> r_parents_;
  // The source involved in these edges
  std::vector<int> r_sources_;
  // The capacities of the edges to the parents, or of the region if no parent
  std::vector<long long> arc_capacities_;

  // What is the best source to move to go from region k1 to region k2?
  std::vector<std::vector<std::priority_queue<MovableSource> > >
      best_interregions_costs_;
  int dijkstra_cnt_;

  // Helper functions

  // Number of regions
  int region_cnt() const {
    assert(sr_costs_.size() == sr_allocations_.size());
    return sr_costs_.size();
  }

  // Update the edge between two regions
  void update_edge(int r1, int r2);
  // Add a source to all heaps of a region; returns if we need to update a path
  bool add_source_to_heaps(int r, int source);
  // Initialize the heaps of a region
  void create_heaps(int reg);

  // Run the shortest path algorithm to update the cost of each region
  void dijkstra_update();

  // Update the edge and returns if we need to rerun Dijkstra
  bool push_edge(int reg, long long flow);
  // Updates a full path when pushing an element; returns if we need to rerun
  // Dijkstra
  bool push_path(int pushed_reg, long long demanded, long long& flow);

 public:
  // Add a new source to the transportation problem; should be done in
  // decreasing order of demand to keep low complexity
  void add_source(int elt_ind);

  CurrentAllocation(const std::vector<long long>& caps,
                    const std::vector<long long>& demands,
                    const std::vector<std::vector<float> >& costs)
      : sr_allocations_(caps.size()),
        sr_costs_(costs),
        s_demands_(demands),
        r_capacities_(caps),
        r_costs_(caps.size(), 0.0),
        r_parents_(caps.size(), NULL_IND),
        r_sources_(caps.size(), NULL_IND),
        arc_capacities_(caps),
        best_interregions_costs_(
            caps.size(),
            std::vector<std::priority_queue<MovableSource> >(caps.size())),
        dijkstra_cnt_(0) {
    assert(caps.size() > 0);
    assert(costs.size() == caps.size());
    dijkstra_update();
  }

  std::vector<std::vector<long long> > get_allocations() const {
    return sr_allocations_;
  }
  int get_iterations_cnt() const { return dijkstra_cnt_; }
};

void CurrentAllocation::update_edge(int r1, int r2) {
  while (not best_interregions_costs_[r1][r2].empty() and
         sr_allocations_[r1][best_interregions_costs_[r1][r2].top().source] ==
             0) {
    best_interregions_costs_[r1][r2].pop();
  }

  if (not best_interregions_costs_[r1][r2].empty()) {
    // There is an edge
    MovableSource cur = best_interregions_costs_[r1][r2].top();
    float new_cost = r_costs_[r2] + cur.cost;
    if (new_cost < r_costs_[r1]) {
      r_costs_[r1] = cur.cost;
      r_sources_[r1] = cur.source;
      r_parents_[r1] = r2;
      arc_capacities_[r1] = sr_allocations_[r1][cur.source];
    }
  }
}

bool CurrentAllocation::add_source_to_heaps(int r, int source) {
  bool need_rerun = false;
  for (int i = 0; i < region_cnt(); ++i) {
    if (i == r) continue;
    best_interregions_costs_[r][i].push(
        MovableSource(source, sr_costs_[i][source] - sr_costs_[r][source]));
    while (sr_allocations_[r][best_interregions_costs_[r][i].top().source] ==
           0) {
      best_interregions_costs_[r][i].pop();
    }
    need_rerun =
        (best_interregions_costs_[r][i].top().source == source) or need_rerun;
  }
  return need_rerun;
}

void CurrentAllocation::create_heaps(int reg) {
  // Get all relevant elements
  std::vector<std::vector<MovableSource> > interregion_costs(region_cnt());
  for (int i = 0; i < sr_allocations_[reg].size(); ++i) {
    if (sr_allocations_[reg][i] > 0) {
      for (int oreg = 0; oreg < region_cnt(); ++oreg) {
        if (oreg == reg) continue;
        interregion_costs[oreg].push_back(
            MovableSource(i, sr_costs_[oreg][i] - sr_costs_[reg][i]));
      }
    }
  }
  // Create the heaps
  for (int oreg = 0; oreg < region_cnt(); ++oreg) {
    best_interregions_costs_[reg][oreg] = std::priority_queue<MovableSource>(
        interregion_costs[oreg].begin(), interregion_costs[oreg].end());
  }
}

// Returns if the path has been modified so that we would need to rerun Dijkstra
bool CurrentAllocation::push_edge(int reg, long long flow) {
  int cur_source = r_sources_[reg];

  // Does this edge allocates a new source in the destination region? If yes,
  // update the corresponding heaps
  bool already_present = sr_allocations_[r_parents_[reg]][cur_source] > 0;

  // Deallocating from the first region is handled by the get_edge function:
  // just substract the flow
  sr_allocations_[reg][cur_source] -= flow;
  sr_allocations_[r_parents_[reg]][cur_source] += flow;

  // The source to be pushed was indeed present in the region
  assert(sr_allocations_[reg][cur_source] >= 0);
  // The region is full, which explains why we need to push
  assert(r_capacities_[reg] == 0);
  // The flow is not bigger than what can be sent
  assert(flow <= arc_capacities_[reg]);

  // Just update the capacity if it turns out that we don't need to run Dijkstra
  arc_capacities_[reg] = sr_allocations_[reg][cur_source];

  if (arc_capacities_[reg] == 0) {
    // The source may have been deleted from a region: rerun Dijkstra at the end
    return true;
  } else if (not already_present and r_capacities_[r_parents_[reg]] == 0) {
    // A new source is allocated to a full region: rerun Dijkstra at the end if
    // it changed the heap's top
    return add_source_to_heaps(r_parents_[reg], cur_source);
  } else {
    // The edge is still present with the same cost and non-zero updated
    // capacity The path still exists: no need to rerun Dijkstra yet
    return false;
  }
}

void CurrentAllocation::dijkstra_update() {
  // Simple case of the regions with remaining capacity
  std::vector<int> visited(region_cnt(), 0);
  int visited_cnt = 0;
  for (int i = 0; i < region_cnt(); ++i) {
    r_sources_[i] = NULL_IND;
    r_parents_[i] = NULL_IND;
    if (r_capacities_[i] > 0) {
      r_costs_[i] = 0.0;
      arc_capacities_[i] = r_capacities_[i];

      visited[i] = 1;
      ++visited_cnt;
    } else {
      r_costs_[i] = std::numeric_limits<float>::infinity();
      arc_capacities_[i] = 0;
    }
  }
  if (visited_cnt == region_cnt()) {
    return;
  }
  // Get the costs for every non-visited region
  for (int i = 0; i < region_cnt(); ++i)
    // For every region that is not visited yet
    if (visited[i] == 0) {
      for (int j = 0; j < region_cnt(); ++j)
        // For every already visited region
        if (visited[j] == 1) {
          // Get the best interregion cost
          update_edge(i, j);
        }
    }
  while (visited_cnt < region_cnt()) {
    // Find the region with the lowest cost to visit; mark it visited
    int best_reg = NULL_IND;
    float best_cost = std::numeric_limits<float>::infinity();
    for (int i = 0; i < region_cnt(); ++i) {
      // For every region that is not visited yet
      if (visited[i] == 0) {
        if (r_costs_[i] < best_cost) {
          best_cost = r_costs_[i];
          best_reg = i;
        }
      }
    }
    if (best_reg == NULL_IND)
      break;  // Some regions are unreachable, typically because they have zero
              // capacity at the beginning
    visited[best_reg] = 1;
    ++visited_cnt;
    // Update the cost for every unvisited region
    for (int i = 0; i < region_cnt(); ++i) {
      if (visited[i] == 0) {
        // For every region that is not visited yet
        update_edge(i, best_reg);
      }
    }
  }
}

bool CurrentAllocation::push_path(int pushed_reg, long long demanded,
                                  long long& flow) {
  // Get the final flow sent, which is smaller than the capacities on the path
  flow = demanded;
  for (int reg = pushed_reg; reg != NULL_IND; reg = r_parents_[reg]) {
    flow = std::min(flow, arc_capacities_[reg]);
  }

  bool rerun_dijkstra = false;
  // Update the path between the regions
  int reg = pushed_reg;
  for (; r_parents_[reg] != NULL_IND; reg = r_parents_[reg]) {
    assert(r_capacities_[reg] == 0);
    rerun_dijkstra = push_edge(reg, flow) or rerun_dijkstra;
  }

  assert(r_capacities_[reg] > 0);
  assert(arc_capacities_[reg] == r_capacities_[reg]);
  assert(r_capacities_[reg] >= flow);

  // Update the capacities at the end
  r_capacities_[reg] -= flow;
  arc_capacities_[reg] -= flow;

  // The last region on the path is the one that satisfies the demand
  if (r_capacities_[reg] ==
      0) {  // If we just consumed the available capacity, it becomes useful to
            // move sources off this region: build the heap
    create_heaps(reg);
    rerun_dijkstra = true;
  }

  assert(flow > 0);

  // If an edge changes cost or a region is full,
  // we need to update the costs, parents, sources and arc_capacities using a
  // Dijkstra but later
  return rerun_dijkstra;
}

void CurrentAllocation::add_source(
    int elt_ind) {  // long long demand, std::vector<float> const & costs){
  for (int i = 0; i < region_cnt(); ++i) {
    sr_allocations_[i].push_back(0);
  }

  bool need_rerun = false;
  long long demand = s_demands_[elt_ind];

  while (demand > 0) {
    // In case we modified the structures earlier
    if (need_rerun) {
      dijkstra_update();
      need_rerun = false;
    }

    ++dijkstra_cnt_;
    int best_reg = NULL_IND;
    float best_cost = std::numeric_limits<float>::infinity();
    for (int reg = 0; reg < region_cnt(); ++reg) {
      // Find the region which gets the source
      if (r_costs_[reg] + sr_costs_[reg][elt_ind] < best_cost) {
        best_reg = reg;
        best_cost = r_costs_[reg] + sr_costs_[reg][elt_ind];
      }
    }
    if (best_reg == NULL_IND) {
      throw std::runtime_error("No reachable region found\n");
    }

    long long flow = 0;
    // Tells us whether we need to update the data structures
    need_rerun = push_path(best_reg, demand, flow);
    demand -= flow;

    // Lazily store the change
    sr_allocations_[best_reg][elt_ind] += flow;
  }

  // Set the source's demand
  for (int i = 0; i < region_cnt(); ++i) {
    if (r_capacities_[i] == 0 and sr_allocations_[i][elt_ind] > 0) {
      need_rerun = add_source_to_heaps(i, elt_ind) or need_rerun;
    }
  }
  // We leave a clean set with correct paths for the next iteration
  if (need_rerun) dijkstra_update();
}

std::vector<std::vector<long long> > solveTransportation(
    const std::vector<long long>& capacities,
    const std::vector<long long>& demands,
    const std::vector<std::vector<float> >& costs) {
  CurrentAllocation transporter(capacities, demands, costs);

  for (int i = 0; i < demands.size(); ++i) {
    transporter.add_source(i);
  }

  return transporter.get_allocations();
}

std::vector<long long> solveTransportation1D(std::vector<t1D_elt> sources,
                                             std::vector<t1D_elt> sinks) {
  /* Description of the algorithm:
   *
   *    For each cell, put it in its optimal region or the last region where a
   * cell is if there is no space in it.
   *
   *    Push all changes in the derivative of the cost function to a priority
   * queue; those changes occur
   *         * when evicting the preceding cell from a region (most such changes
   * are 0 and not considered, hence the complexity)
   *         * when moving to a non-full region
   *
   *    While the new cell overlaps with a new region, get the new slope
   * (derivative) at this point and push all preceding cell until this region is
   * freed or the slope becomes 0 (in which case the new region is now occupied)
   */

  struct Bound {
    long long pos;
    int slope_diff;
    bool operator<(const Bound o) const { return pos < o.pos; }
  };

  std::priority_queue<Bound> bounds;
  std::vector<long long> constraining_pos;
  std::vector<long long> prev_cap(1, 0), prev_dem(1, 0);
  for (auto const s : sinks) {
    prev_cap.push_back(s.second + prev_cap.back());
  }
  for (auto const s : sources) {
    prev_dem.push_back(s.second + prev_dem.back());
  }
  // The sinks have enough capacity to hold the whole demand
  assert(prev_cap.back() >= prev_dem.back());

  const long long min_abs_pos = 0,
                  max_abs_pos = prev_cap.back() - prev_dem.back();
  assert(min_abs_pos <= max_abs_pos);

  auto push_bound = [&](long long p, int s) {
    assert(s >= 0);
    if (p > min_abs_pos) {
      Bound b;
      b.pos = p;
      b.slope_diff = s;
      bounds.push(b);
    }
  };

  // Distance to the right - distance to the left
  auto get_slope = [&](int src, int boundary) {
    assert(boundary + 1 < sinks.size());
    assert(src < sources.size());
    return std::abs(sources[src].first - sinks[boundary + 1].first) -
           std::abs(sources[src].first - sinks[boundary].first);
  };

  long long cur_abs_pos = min_abs_pos;
  int opt_r = 0, next_r = 0, first_free_r = 0;

  for (int i = 0; i < sources.size(); ++i) {
    // Update the optimal region
    while (opt_r + 1 < sinks.size() and
           (sinks[opt_r].first + sinks[opt_r + 1].first) / 2 <
               sources[i].first) {
      ++opt_r;
    }
    // Update the next region
    int prev_next_r = next_r;
    while (next_r < sinks.size() and sinks[next_r].first <= sources[i].first) {
      ++next_r;
    }

    int dest_reg = std::max(first_free_r, opt_r);
    assert(dest_reg < sinks.size());

    if (i > 0) {
      // Push bounds due to changing the source crossing the boundary j/j+1
      // Linear amortized complexity accross all sources (next_r grows)
      // get_slope(i-1, j) - get_slope(i, j) == 0 if j >= next_r
      // get_slope(i-1, j) - get_slope(i, j) == 0 if j < prev_next_r-1

      for (int j = std::max(prev_next_r, 1) - 1;
           j < std::min(first_free_r, next_r + 1); ++j) {
        assert(get_slope(i, j) <= get_slope(i - 1, j));
        push_bound(prev_cap[j + 1] - prev_dem[i],
                   get_slope(i - 1, j) - get_slope(i, j));
      }
    }
    // Add the bounds due to crossing the boundaries alone
    for (int j = first_free_r; j < opt_r; ++j) {
      assert(get_slope(i, j) <= 0);
      push_bound(prev_cap[j + 1] - prev_dem[i], -get_slope(i, j));
    }

    first_free_r = std::max(first_free_r, opt_r);
    // Just after the previous cell or at the beginning of the destination
    // region
    long long this_abs_pos =
        std::max(cur_abs_pos, prev_cap[first_free_r] - prev_dem[i]);

    while (first_free_r + 1 < sinks.size() and
           this_abs_pos > std::max(prev_cap[first_free_r + 1] - prev_dem[i + 1],
                                   min_abs_pos)) {
      // Absolute position that wouldn't make the cell fit in the region, and we
      // are not in the last region yet
      long long end_pos =
          std::max(prev_cap[first_free_r + 1] - prev_dem[i + 1], min_abs_pos);

      int add_slope = get_slope(i, first_free_r);
      int slope = add_slope;

      while (not bounds.empty() and slope >= 0 and bounds.top().pos > end_pos) {
        this_abs_pos = bounds.top().pos;
        slope -= bounds.top().slope_diff;
        bounds.pop();
      }
      if (slope >=
          0) {  // We still push: the cell completely escapes the region
        this_abs_pos = end_pos;
        push_bound(end_pos, add_slope - slope);
      } else {  // Ok, absorbed the whole slope: push what remains and we still
                // occupy the next region
        push_bound(this_abs_pos, -slope);
        ++first_free_r;
      }
    }
    cur_abs_pos = this_abs_pos;
    constraining_pos.push_back(this_abs_pos);
  }

  assert(constraining_pos.size() == sources.size());
  if (not constraining_pos.empty()) {
    // Calculate the final constraining_pos
    constraining_pos.back() = std::min(max_abs_pos, constraining_pos.back());
  }

  std::partial_sum(
      constraining_pos.rbegin(), constraining_pos.rend(),
      constraining_pos.rbegin(),
      [](long long a, long long b) -> long long { return std::min(a, b); });

  for (int i = 0; i < constraining_pos.size(); ++i) {
    constraining_pos[i] += prev_dem[i];
  }

  return constraining_pos;
}

}  // namespace coloquinte