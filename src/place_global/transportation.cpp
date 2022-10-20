
#include <cassert>
#include <limits>
#include <queue>
#include <stdexcept>
#include <vector>

class current_allocation {
  static constexpr int null_ind = -1;

  // Internal data structures

  // Priority queue element to determine the source to be used between regions
  struct movable_source {
    int source;
    float cost;
    bool operator<(movable_source const o) const {
      return cost > o.cost  // Sorted by cost
             || (cost == o.cost &&
                 source < o.source);  // And by index to limit the number of
                                      // fractional elements between two regions
    }
    movable_source(int s, float c) : source(s), cost(c) {}
  };

  // Member data

  // The current state
  std::vector<std::vector<long long> >
      sr_allocations;  // For each region, for each source, the capacity
                       // allocated by the region
  std::vector<std::vector<float> >
      sr_costs;                      // The costs from a region to a source
  std::vector<long long> s_demands;  // The demands of the sources
  std::vector<long long>
      r_capacities;  // The remaining capacities of the regions

  // Shortest path data
  // The costs of allocating to a region
  std::vector<float> r_costs;
  // The parents of the regions i.e. the regions where we push sources first (or
  // null_ind)
  std::vector<int> r_parents;
  // The source involved in these edges
  std::vector<int> r_sources;
  // The capacities of the edges to the parents, or of the region if no parent
  std::vector<long long> arc_capacities;

  // What is the best source to move to go from region k1 to region k2?
  std::vector<std::vector<std::priority_queue<movable_source> > >
      best_interregions_costs;
  int dijkstra_cnt;

  // Helper functions

  // Number of regions
  int region_cnt() const {
    assert(sr_costs.size() == sr_allocations.size());
    return sr_costs.size();
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

  current_allocation(std::vector<long long> caps,
                     std::vector<long long> demands,
                     std::vector<std::vector<float> > costs)
      : sr_allocations(caps.size()),
        sr_costs(costs),
        s_demands(demands),
        r_capacities(caps),
        r_costs(caps.size(), 0.0),
        r_parents(caps.size(), null_ind),
        r_sources(caps.size(), null_ind),
        arc_capacities(caps),
        best_interregions_costs(
            caps.size(),
            std::vector<std::priority_queue<movable_source> >(caps.size())),
        dijkstra_cnt(0) {
    assert(caps.size() > 0);
    assert(costs.size() == caps.size());
    dijkstra_update();
  }

  std::vector<std::vector<long long> > get_allocations() const {
    return sr_allocations;
  }
  int get_iterations_cnt() const { return dijkstra_cnt; }
};

void current_allocation::update_edge(int r1, int r2) {
  while (not best_interregions_costs[r1][r2].empty() and
         sr_allocations[r1][best_interregions_costs[r1][r2].top().source] ==
             0) {
    best_interregions_costs[r1][r2].pop();
  }

  if (not best_interregions_costs[r1][r2].empty()) {
    // There is an edge
    movable_source cur = best_interregions_costs[r1][r2].top();
    float new_cost = r_costs[r2] + cur.cost;
    if (new_cost < r_costs[r1]) {
      r_costs[r1] = cur.cost;
      r_sources[r1] = cur.source;
      r_parents[r1] = r2;
      arc_capacities[r1] = sr_allocations[r1][cur.source];
    }
  }
}

bool current_allocation::add_source_to_heaps(int r, int source) {
  bool need_rerun = false;
  for (int i = 0; i < region_cnt(); ++i) {
    if (i == r) continue;
    best_interregions_costs[r][i].push(
        movable_source(source, sr_costs[i][source] - sr_costs[r][source]));
    while (sr_allocations[r][best_interregions_costs[r][i].top().source] == 0) {
      best_interregions_costs[r][i].pop();
    }
    need_rerun =
        (best_interregions_costs[r][i].top().source == source) or need_rerun;
  }
  return need_rerun;
}

void current_allocation::create_heaps(int reg) {
  // Get all relevant elements
  std::vector<std::vector<movable_source> > interregion_costs(region_cnt());
  for (int i = 0; i < sr_allocations[reg].size(); ++i) {
    if (sr_allocations[reg][i] > 0) {
      for (int oreg = 0; oreg < region_cnt(); ++oreg) {
        if (oreg == reg) continue;
        interregion_costs[oreg].push_back(
            movable_source(i, sr_costs[oreg][i] - sr_costs[reg][i]));
      }
    }
  }
  // Create the heaps
  for (int oreg = 0; oreg < region_cnt(); ++oreg) {
    best_interregions_costs[reg][oreg] = std::priority_queue<movable_source>(
        interregion_costs[oreg].begin(), interregion_costs[oreg].end());
  }
}

// Returns if the path has been modified so that we would need to rerun Dijkstra
bool current_allocation::push_edge(int reg, long long flow) {
  int cur_source = r_sources[reg];

  // Does this edge allocates a new source in the destination region? If yes,
  // update the corresponding heaps
  bool already_present = sr_allocations[r_parents[reg]][cur_source] > 0;

  // Deallocating from the first region is handled by the get_edge function:
  // just substract the flow
  sr_allocations[reg][cur_source] -= flow;
  sr_allocations[r_parents[reg]][cur_source] += flow;

  assert(sr_allocations[reg][cur_source] >=
         0);  // The source to be pushed was indeed present in the region
  assert(r_capacities[reg] ==
         0);  // The region is full, which explains why we need to push
  assert(flow <=
         arc_capacities[reg]);  // The flow is not bigger than what can be sent

  arc_capacities[reg] =
      sr_allocations[reg]
                    [cur_source];  // Just update the capacity if it turns out
                                   // that we don't need to run Dijkstra

  if (arc_capacities[reg] == 0) {
    // The source may have been deleted from a region: rerun Dijkstra at the end
    return true;
  } else if (not already_present and r_capacities[r_parents[reg]] == 0) {
    // A new source is allocated to a full region: rerun Dijkstra at the end if
    // it changed the heap's top
    return add_source_to_heaps(r_parents[reg], cur_source);
  } else {
    // The edge is still present with the same cost and non-zero updated
    // capacity The path still exists: no need to rerun Dijkstra yet
    return false;
  }
}

void current_allocation::dijkstra_update() {
  // Simple case of the regions with remaining capacity
  std::vector<int> visited(region_cnt(), 0);
  int visited_cnt = 0;
  for (int i = 0; i < region_cnt(); ++i) {
    r_sources[i] = null_ind;
    r_parents[i] = null_ind;
    if (r_capacities[i] > 0) {
      r_costs[i] = 0.0;
      arc_capacities[i] = r_capacities[i];

      visited[i] = 1;
      ++visited_cnt;
    } else {
      r_costs[i] = std::numeric_limits<float>::infinity();
      arc_capacities[i] = 0;
    }
  }
  // if(visited_cnt <= 0) throw std::runtime_error("Capacity problem: no region
  // has been marked as reachable\n");
  if (visited_cnt == region_cnt()) {
    return;
  }
  // Get the costs for every non-visited region
  for (int i = 0; i < region_cnt(); ++i)
    if (visited[i] == 0) {  // For every region that is not visited yet
      for (int j = 0; j < region_cnt(); ++j)
        if (visited[j] == 1) {  // For every already visited region
          // Get the best interregion cost
          update_edge(i, j);
        }
    }
  while (visited_cnt < region_cnt()) {
    // Find the region with the lowest cost to visit; mark it visited
    int best_reg = null_ind;
    float best_cost = std::numeric_limits<float>::infinity();
    for (int i = 0; i < region_cnt(); ++i)
      if (visited[i] == 0) {  // For every region that is not visited yet
        if (r_costs[i] < best_cost) {
          best_cost = r_costs[i];
          best_reg = i;
        }
      }
    if (best_reg == null_ind)
      break;  // Some regions are unreachable, typically because they have zero
              // capacity at the beginning
    visited[best_reg] = 1;
    ++visited_cnt;
    // Update the cost for every unvisited region
    for (int i = 0; i < region_cnt(); ++i)
      if (visited[i] == 0) {  // For every region that is not visited yet
        update_edge(i, best_reg);
      }
  }
}

bool current_allocation::push_path(int pushed_reg, long long demanded,
                                   long long& flow) {
  // Get the final flow sent, which is smaller than the capacities on the path
  flow = demanded;
  for (int reg = pushed_reg; reg != null_ind; reg = r_parents[reg]) {
    flow = std::min(flow, arc_capacities[reg]);
  }

  bool rerun_dijkstra = false;
  // Update the path between the regions
  int reg = pushed_reg;
  for (; r_parents[reg] != null_ind; reg = r_parents[reg]) {
    assert(r_capacities[reg] == 0);
    rerun_dijkstra = push_edge(reg, flow) or rerun_dijkstra;
  }

  assert(r_capacities[reg] > 0);
  assert(arc_capacities[reg] == r_capacities[reg]);
  assert(r_capacities[reg] >= flow);

  // Update the capacities at the end
  r_capacities[reg] -= flow;
  arc_capacities[reg] -= flow;

  // The last region on the path is the one that satisfies the demand
  if (r_capacities[reg] ==
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

void current_allocation::add_source(
    int elt_ind) {  // long long demand, std::vector<float> const & costs){
  for (int i = 0; i < region_cnt(); ++i) {
    sr_allocations[i].push_back(0);
  }

  bool need_rerun = false;
  long long demand = s_demands[elt_ind];

  while (demand > 0) {
    // In case we modified the structures earlier
    if (need_rerun) {
      dijkstra_update();
      need_rerun = false;
    }

    ++dijkstra_cnt;
    int best_reg = null_ind;
    float best_cost = std::numeric_limits<float>::infinity();
    for (int reg = 0; reg < region_cnt(); ++reg) {
      // Find the region which gets the source
      if (r_costs[reg] + sr_costs[reg][elt_ind] < best_cost) {
        best_reg = reg;
        best_cost = r_costs[reg] + sr_costs[reg][elt_ind];
      }
    }
    if (best_reg == null_ind) {
      throw std::runtime_error("No reachable region found\n");
    }

    long long flow = 0;
    // Tells us whether we need to update the data structures
    need_rerun = push_path(best_reg, demand, flow);
    demand -= flow;

    // Lazily store the change
    sr_allocations[best_reg][elt_ind] += flow;
  }

  // Set the source's demand
  for (int i = 0; i < region_cnt(); ++i) {
    if (r_capacities[i] == 0 and sr_allocations[i][elt_ind] > 0) {
      need_rerun = add_source_to_heaps(i, elt_ind) or need_rerun;
    }
  }
  // We leave a clean set with correct paths for the next iteration
  if (need_rerun) dijkstra_update();
}

std::vector<std::vector<long long> > solveTransportation(std::vector<long long> const & capacities, std::vector<long long> const & demands, std::vector<std::vector<float> > const & costs){
    current_allocation transporter(capacities, demands, costs);

    for(int i=0; i<demands.size(); ++i){
        transporter.add_source(i);
    }

    return transporter.get_allocations();
}
