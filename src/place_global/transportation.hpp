#pragma once

#include <queue>
#include <vector>

namespace coloquinte {

/**
 * @brief A transportation solver for problems with few sinks and many sources
 */
class TransportationProblem {
 public:
  /**
   * @brief Initialize the datastructure
   */
  TransportationProblem(const std::vector<long long>& capacities,
                       const std::vector<long long>& demands,
                       const std::vector<std::vector<long long> >& costs);

  TransportationProblem(const std::vector<long long>& capacities,
                       const std::vector<long long>& demands,
                       const std::vector<std::vector<float> >& costs);

  /**
   * @brief Number of sources in the problem
   */
  int nbSources() const { return demands_.size(); }

  /**
   * @brief Number of sinks in the problem
   */
  int nbSinks() const { return capacities_.size(); }

  /**
   * @brief Compute the total demand from the sources
   */
  long long totalDemand() const;

  /**
   * @brief Compute the total capacity at the sinks
   */
  long long totalCapacity() const;

  /**
   * @brief Demand of all sources
   */
  const std::vector<long long>& demands() const { return demands_; }

  /**
   * @brief Capacities for all sinks
   */
  const std::vector<long long>& capacities() const { return capacities_; }

  /**
   * @brief Costs between sinks and sources
   */
  const std::vector<std::vector<long long> >& costs() const { return costs_; }

  /**
   * @brief Allocations between sinks and sources
   */
  const std::vector<std::vector<long long> >& allocations() const {
    return allocations_;
  }

  /**
   * @brief Demand for this source
   */
  long long demand(int src) const { return demands_[src]; }

  /**
   * @brief Capacity for this sink
   */
  long long capacity(int snk) const { return capacities_[snk]; }

  /**
   * @brief Allocation between one sink and source
   */
  long long allocation(int snk, int src) const {
    return allocations_[snk][src];
  }

  /**
   * @brief Cost between one sink and source
   */
  long long cost(int snk, int src) const { return costs_[snk][src]; }

  /**
   * @brief Cost moving a source between sinks
   */
  long long movingCost(int src, int snk1, int snk2) const;

  /**
   * @brief Cost between one sink and source, with scaling removed
   */
  double originalCost(int snk, int src) const {
    return costs_[snk][src] / conversionFactor_;
  }

  /**
   * @brief Demand allocated in the current solution for this source
   */
  long long allocatedDemand(int src) const;

  /**
   * @brief Capacity allocated in the current solution for this sink
   */
  long long allocatedCapacity(int snk) const;

  /**
   * @brief Cost of the current allocation
   */
  double allocationCost() const;

  /**
   * @brief Add a dummy source with 0 cost
   */
  void addSource(long long demand);

  /**
   * @brief Add a dummy sink with 0 cost
   */
  void addSink(long long capa);

  /**
   * @brief Return true if the problem is balanced
   */
  bool isBalanced() const { return totalDemand() == totalCapacity(); }

  /**
   * @brief Return true if the allocation is feasible
   */
  bool isFeasible() const;

  /**
   * @brief Ensure that the capacity is greater or equal than the demand
   */
  void increaseCapacity();

  /**
   * @brief Ensure that the demand is greater or equal than the capacity
   */
  void increaseDemand();

  /**
   * @brief Ensure that the capacity is greater or equal than the demand with an additional sink
   */
  void addDummyCapacity();

  /**
   * @brief Ensure that the demand is greater or equal than the capacity with an additional source
   */
  void addDummyDemand();

  /**
   * @brief Ensure a feasible solution (demands + capacities met)
   */
  void makeFeasible();

  /**
   * @brief Initialize with an initial solution
   */
  void setAllocations(const std::vector<std::vector<long long> >& allocations);

  /**
   * @brief Initialize with an initial solution from a target assignment
   */
  void setAssignment(const std::vector<int>& assignment);

  /**
   * @brief Obtain an assignment corresponding to the allocation (one bin index
   * for each source)
   */
  std::vector<int> toAssignment() const;

  /**
   * @brief Reset the allocation to all-zeros
   */
  void resetAllocations();

  /**
   * @brief Solve the transportation problem
   */
  void solve();

  /**
   * @brief Check the datastructure
   */
  void check() const;

 private:
  /**
   * @brief Convert floating-point costs to a fixed-point equivalent
   */
  void costsFromIntegers(const std::vector<std::vector<float> >& costs);

 private:
  // Problem data
  std::vector<long long> demands_;
  std::vector<long long> capacities_;
  std::vector<std::vector<long long> > costs_;

  // Current solution
  std::vector<std::vector<long long> > allocations_;

  // Optional conversion factor between internal fixed-point costs and external
  double conversionFactor_;

  friend class TransportationCycleCanceling;
  friend class TransportationSuccessiveShortestPath;
};

}  // namespace coloquinte