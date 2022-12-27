
#include "place_global/transportation.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <limits>
#include <numeric>
#include <queue>
#include <stdexcept>
#include <utility>
#include <vector>

namespace coloquinte {
struct CostElt {
  CostElt(CostType cost, int elt) : elt(elt), cost(cost) {}

  bool operator<(const CostElt o) const { return cost > o.cost; }

  int elt;
  CostType cost;
};
using PrioQueue = std::priority_queue<CostElt>;

class TransportationSuccessiveShortestPath {
 public:
  /**
   * @brief Solve the problem using the cycle canceling method
   */
  static void solve(TransportationProblem& pb);

  /**
   * @brief Init from a transportation problem
   */
  TransportationSuccessiveShortestPath(TransportationProblem& pb);

  /**
   * @brief Run the whole algorithm
   */
  void run();

  /**
   * @brief Find the best source to be sent between two sinks
   */
  inline int sentSource(int snk1, int snk2) const {
    assert(snk1 != snk2);
    return queues_[snk1][snk2].top().elt;
  }

  /**
   * @brief Find the cost of sending the best source between two sinks
   */
  inline CostType movingCost(int snk1, int snk2) const {
    if (snk1 == snk2) return 0LL;
    return queues_[snk1][snk2].top().cost;
  }

  /**
   * @brief Available quantity to send between two sinks
   */
  inline DemandType sentQuantity(int snk1, int snk2) const {
    int src = sentSource(snk1, snk2);
    return pb_.allocation(snk1, src);
  }

  /**
   * @brief Find a good ordering for the sources (high demand first)
   */
  std::vector<int> sortedSourcesByDemand() const;

  /**
   * @brief Find a good ordering for the sources (high cost discrepancy first)
   */
  std::vector<int> sortedSourcesByCost() const;

  /**
   * @brief Find the best sink to send the source to
   */
  int bestSink(int src) const;

  /**
   * @brief Send as much as possible of the source to the sink; return how much
   * was sent
   */
  void sendSource(int src);

  /**
   * @brief Send as much as possible of the source to the sink; return how much
   * was sent
   */
  DemandType sendSource(int src, int sink, DemandType quantity);

 private:
  /**
   * @brief Setup the priority queues for a full sink
   */
  void initQueues(int sink);

  /**
   * @brief Update a sink after a potential source removal
   */
  void updateSinkQueues(int sink, int src);

  /**
   * @brief Update a sink before a potential source addition
   */
  void updateDestQueues(int sink, int src);

  /**
   * @brief Compute all costs
   */
  void updateTree();

  /**
   * @brief Update the costs after a change in parent
   */
  void updateCosts(int sink);

 private:
  TransportationProblem& pb_;

  // Queues between sinks
  std::vector<std::vector<PrioQueue> > queues_;

  std::vector<DemandType> remainingCapa_;
  std::vector<CostType> sendingCost_;
  std::vector<int> sinkParent_;
};

TransportationProblem::TransportationProblem(
    const std::vector<DemandType>& capacities,
    const std::vector<DemandType>& demands,
    const std::vector<std::vector<CostType> >& costs)
    : capacities_(capacities), demands_(demands), costs_(costs) {
  resetAllocations();
  check();
}

TransportationProblem::TransportationProblem(
    const std::vector<DemandType>& capacities,
    const std::vector<DemandType>& demands,
    const std::vector<std::vector<float> >& costs)
    : capacities_(capacities), demands_(demands) {
  costsFromIntegers(costs);
  resetAllocations();
  check();
}

void TransportationProblem::resetAllocations() {
  allocations_.assign(nbSinks(), std::vector<DemandType>(nbSources(), 0LL));
}

void TransportationProblem::costsFromIntegers(
    const std::vector<std::vector<float> >& costs) {
  float maxVal = 1.0e-8f;
  for (auto& c : costs) {
    for (float d : c) {
      maxVal = std::max(d, maxVal);
    }
  }
  double maxLong = static_cast<double>(std::numeric_limits<CostType>::max());
  conversionFactor_ = maxLong / maxVal;
  conversionFactor_ /= 4.0;
  conversionFactor_ /= costs.size();
  costs_.clear();
  costs_.resize(costs.size());
  for (int i = 0; i < costs.size(); ++i) {
    for (int j = 0; j < costs[i].size(); ++j) {
      costs_[i].push_back(
          static_cast<CostType>(std::round(costs[i][j] * conversionFactor_)));
    }
  }
}

void TransportationProblem::check() const {
  if (demands_.size() != nbSources()) {
    throw std::runtime_error("Inconsistant demands");
  }
  if (capacities_.size() != nbSinks()) {
    throw std::runtime_error("Inconsistant capacities");
  }
  for (DemandType c : demands_) {
    if (c <= 0) {
      throw std::runtime_error("Demands must be non-negative");
    }
  }
  for (DemandType c : capacities_) {
    if (c <= 0) {
      throw std::runtime_error("Capacities must be non-negative");
    }
  }

  if (costs_.size() != nbSinks()) {
    throw std::runtime_error(
        "Inconsistant cost size (first dimension is sinks)");
  }
  for (const auto& c : costs_) {
    if (c.size() != nbSources()) {
      throw std::runtime_error(
          "Inconsistant cost size (second dimension is sources)");
    }
  }

  if (allocations_.size() != nbSinks()) {
    throw std::runtime_error(
        "Inconsistant allocation size (first dimension is sinks)");
  }
  for (const auto& c : allocations_) {
    if (c.size() != nbSources()) {
      throw std::runtime_error(
          "Inconsistant allocation size (second dimension is sources)");
    }
  }
}

DemandType TransportationProblem::totalDemand() const {
  return std::accumulate(demands_.begin(), demands_.end(), 0LL);
}

DemandType TransportationProblem::totalCapacity() const {
  return std::accumulate(capacities_.begin(), capacities_.end(), 0LL);
}

void TransportationProblem::addSource(DemandType demand) {
  demands_.push_back(demand);
  for (auto& c : costs_) {
    c.push_back(0.0f);
  }
  for (auto& c : allocations_) {
    c.push_back(0LL);
  }
}

void TransportationProblem::addSink(DemandType capa) {
  capacities_.push_back(capa);
  costs_.emplace_back(nbSources(), 0.0f);
  allocations_.emplace_back(nbSources(), 0LL);
}

DemandType TransportationProblem::allocatedCapacity(int snk) const {
  DemandType tot = 0LL;
  for (int src = 0; src < nbSources(); ++src) {
    tot += allocations_[snk][src];
  }
  return tot;
}

DemandType TransportationProblem::allocatedDemand(int src) const {
  DemandType tot = 0LL;
  for (int snk = 0; snk < nbSinks(); ++snk) {
    tot += allocations_[snk][src];
  }
  return tot;
}

double TransportationProblem::allocationCost() const {
  double ret = 0.0;
  for (int snk = 0; snk < nbSinks(); ++snk) {
    for (int src = 0; src < nbSources(); ++src) {
      ret += allocations_[snk][src] * static_cast<double>(costs_[snk][src]);
    }
  }
  return ret / totalDemand() / conversionFactor_;
}

inline CostType TransportationProblem::movingCost(int src, int snk1,
                                                  int snk2) const {
  return costs_[snk2][src] - costs_[snk1][src];
}

bool TransportationProblem::isFeasible() const {
  if (totalDemand() > totalCapacity()) {
    throw std::runtime_error(
        "Cannot check feasibility of an unbalanced problem");
  }
  for (int src = 0; src < nbSources(); ++src) {
    if (allocatedDemand(src) != demand(src)) {
      return false;
    }
  }
  for (int snk = 0; snk < nbSinks(); ++snk) {
    if (allocatedCapacity(snk) > capacity(snk)) {
      return false;
    }
  }
  return true;
}

void TransportationProblem::addDummyCapacity() {
  DemandType dem = totalDemand();
  DemandType cap = totalCapacity();
  if (dem > cap) {
    addSink(dem - cap);
  }
}

void TransportationProblem::addDummyDemand() {
  DemandType dem = totalDemand();
  DemandType cap = totalCapacity();
  if (dem < cap) {
    addSource(cap - dem);
  }
}

void TransportationProblem::increaseCapacity() {
  DemandType missing = totalDemand() - totalCapacity();
  if (missing <= 0LL) {
    return;
  }
  DemandType added = missing / nbSinks();
  for (int i = 0; i < nbSinks(); ++i) {
    capacities_[i] += added;
  }
  missing = missing - added * nbSinks();
  assert(missing >= 0LL && missing < nbSinks());
  for (int i = 0; i < missing; ++i) {
    capacities_[i] += 1;
  }
}

void TransportationProblem::increaseDemand() {
  DemandType missing = totalCapacity() - totalDemand();
  if (missing <= 0LL) {
    return;
  }
  DemandType added = missing / nbSources();
  for (int i = 0; i < nbSources(); ++i) {
    demands_[i] += added;
  }
  missing = missing - added * nbSources();
  assert(missing >= 0LL && missing < nbSources());
  for (int i = 0; i < missing; ++i) {
    demands_[i] += 1;
  }
}

void TransportationProblem::makeFeasible() {
  if (totalDemand() > totalCapacity()) {
    throw std::runtime_error(
        "Cannot make a problem feasible if demand is greater than capacity");
  }
  std::vector<DemandType> remainingCapa = capacities_;
  for (int src = 0; src < nbSources(); ++src) {
    DemandType demand = demands_[src];
    DemandType remaining = demand;
    // Try to follow the initial solution given
    for (int snk = 0; snk < nbSinks(); ++snk) {
      DemandType alloc = allocations_[snk][src];
      alloc = std::min(alloc, remaining);
      alloc = std::min(alloc, remainingCapa[snk]);
      allocations_[snk][src] = alloc;
      remaining -= alloc;
      remainingCapa[snk] -= alloc;
    }
    // Allocate greedily what's left
    for (int snk = 0; snk < nbSinks(); ++snk) {
      DemandType alloc = std::min(remaining, remainingCapa[snk]);
      allocations_[snk][src] += alloc;
      remaining -= alloc;
      remainingCapa[snk] -= alloc;
    }
    if (remaining < 0LL) {
      // Should never happen if preconditions are met
      throw std::runtime_error("Could not satisfy demand for one of the cells");
    }
  }
}

void TransportationProblem::setAllocations(
    const std::vector<std::vector<DemandType> >& allocations) {
  allocations_ = allocations;
  check();
}

void TransportationProblem::setAssignment(const std::vector<int>& assignment) {
  allocations_.assign(nbSinks(), std::vector<DemandType>(nbSources(), 0LL));
  if (assignment.size() > nbSources()) {
    throw std::runtime_error(
        "Assignment should be no larger than the number of sources");
  }
  for (int src = 0; src < assignment.size(); ++src) {
    int sink = assignment[src];
    if (sink < 0 || sink >= nbSinks()) {
      throw std::runtime_error("Assignment should be a valid sink index");
    }
    allocations_[sink][src] = demands_[src];
  }
}

std::vector<int> TransportationProblem::toAssignment() const {
  std::vector<int> ret(nbSources());
  for (int src = 0; src < nbSources(); ++src) {
    int bestSink = 0;
    DemandType bestAlloc = -1LL;
    for (int sink = 0; sink < nbSinks(); ++sink) {
      if (allocations_[sink][src] > bestAlloc) {
        bestSink = sink;
        bestAlloc = allocations_[sink][src];
      }
    }
    ret[src] = bestSink;
  }
  return ret;
}

void TransportationProblem::solve() {
  TransportationSuccessiveShortestPath::solve(*this);
}

void TransportationSuccessiveShortestPath::solve(TransportationProblem& pb) {
  TransportationSuccessiveShortestPath solver(pb);
  solver.run();
}

TransportationSuccessiveShortestPath::TransportationSuccessiveShortestPath(
    TransportationProblem& pb)
    : pb_(pb) {
  remainingCapa_ = pb_.capacities_;
  sendingCost_.assign(pb_.nbSinks(), 0LL);
  sinkParent_.assign(pb_.nbSinks(), -1);
  queues_.resize(pb_.nbSinks());
}

void TransportationSuccessiveShortestPath::run() {
  std::vector<int> sources = sortedSourcesByDemand();
  pb_.resetAllocations();
  for (int src : sources) {
    sendSource(src);
  }
}

std::vector<int> TransportationSuccessiveShortestPath::sortedSourcesByDemand()
    const {
  std::vector<std::pair<DemandType, int> > sources;
  for (int i = 0; i < pb_.nbSources(); ++i) {
    sources.emplace_back(-pb_.demand(i), i);
  }
  std::sort(sources.begin(), sources.end());
  std::vector<int> ret;
  for (auto [d, i] : sources) {
    ret.push_back(i);
  }
  return ret;
}

int TransportationSuccessiveShortestPath::bestSink(int src) const {
  int ret = 0;
  CostType bestCost = std::numeric_limits<CostType>::max();
  for (int i = 0; i < pb_.nbSinks(); ++i) {
    CostType cost = sendingCost_[i] + pb_.cost(i, src);
    if (cost < bestCost) {
      bestCost = cost;
      ret = i;
    }
  }
  return ret;
}

void TransportationSuccessiveShortestPath::sendSource(int src) {
  DemandType remaining = pb_.demand(src);
  while (remaining > 0LL) {
    int sink = bestSink(src);
    DemandType sent = sendSource(src, sink, remaining);
    assert(sent > 0LL);
    remaining -= sent;
  }
}

void TransportationSuccessiveShortestPath::updateTree() {
  std::vector<char> toVisit(pb_.nbSinks(), false);
  sendingCost_.assign(pb_.nbSinks(), std::numeric_limits<CostType>::max());
  for (int i = 0; i < pb_.nbSinks(); ++i) {
    if (remainingCapa_[i] > 0LL) {
      sendingCost_[i] = 0LL;
      toVisit[i] = true;
    }
    sinkParent_[i] = -1;
  }

  while (true) {
    int bestVisit = -1;
    CostType bestCost = std::numeric_limits<CostType>::max();
    for (int i = 0; i < pb_.nbSinks(); ++i) {
      if (!toVisit[i]) {
        continue;
      }
      if (sendingCost_[i] < bestCost) {
        bestCost = sendingCost_[i];
        bestVisit = i;
      }
    }
    if (bestVisit == -1) {
      break;
    }
    for (int i = 0; i < pb_.nbSinks(); ++i) {
      if (remainingCapa_[i] > 0LL) {
        continue;
      }
      CostType newCost = movingCost(i, bestVisit) + sendingCost_[bestVisit];
      if (newCost < sendingCost_[i]) {
        sinkParent_[i] = bestVisit;
        sendingCost_[i] = newCost;
        toVisit[i] = true;
      }
    }
    toVisit[bestVisit] = false;
  }
}

DemandType TransportationSuccessiveShortestPath::sendSource(
    int src, int sink, DemandType quantity) {
  // How much we can send
  assert(quantity > 0LL);
  DemandType maxSent = quantity;
  int snk1 = sink;
  while (sinkParent_[snk1] != -1) {
    int snk2 = sinkParent_[snk1];
    maxSent = std::min(maxSent, sentQuantity(snk1, snk2));
    assert(maxSent > 0LL);
    snk1 = snk2;
  }
  maxSent = std::min(maxSent, remainingCapa_[snk1]);
  assert(maxSent > 0LL);

  // Actually send
  snk1 = sink;
  int sentSrc = src;
  bool needUpdate = false;
  while (sinkParent_[snk1] != -1) {
    assert(remainingCapa_[snk1] == 0LL);
    int snk2 = sinkParent_[snk1];
    CostType oldCost = movingCost(snk1, snk2);
    updateDestQueues(snk1, sentSrc);
    pb_.allocations_[snk1][sentSrc] += maxSent;
    sentSrc = sentSource(snk1, snk2);
    pb_.allocations_[snk1][sentSrc] -= maxSent;
    updateSinkQueues(snk1, sentSrc);
    CostType newCost = movingCost(snk1, snk2);
    needUpdate |= newCost > oldCost;
    snk1 = snk2;
  }
  pb_.allocations_[snk1][sentSrc] += maxSent;

  remainingCapa_[snk1] -= maxSent;
  if (remainingCapa_[snk1] == 0LL) {
    initQueues(snk1);
    needUpdate = true;
  }
  if (needUpdate) {
    updateTree();
  }

  return maxSent;
}

void TransportationSuccessiveShortestPath::initQueues(int sink) {
  assert(remainingCapa_[sink] == 0LL);
  queues_[sink].clear();
  queues_[sink].resize(pb_.nbSinks());
  std::vector<int> sources;
  for (int src = 0; src < pb_.nbSources(); ++src) {
    if (pb_.allocation(sink, src) == 0LL) {
      continue;
    }
    sources.push_back(src);
  }
  for (int dest = 0; dest < pb_.nbSinks(); ++dest) {
    if (sink == dest) continue;
    std::vector<CostElt> elements;
    for (int src : sources) {
      elements.emplace_back(pb_.movingCost(src, sink, dest), src);
    }
    queues_[sink][dest] = PrioQueue(elements.begin(), elements.end());
  }
}

void TransportationSuccessiveShortestPath::updateSinkQueues(int sink, int src) {
  if (pb_.allocations_[sink][src] != 0LL) {
    return;
  }
  for (int dst = 0; dst < pb_.nbSinks(); ++dst) {
    if (dst == sink) continue;
    while (true) {
      if (queues_[sink][dst].empty()) {
        break;
      }
      int bestSource = sentSource(sink, dst);
      if (pb_.allocations_[sink][bestSource] != 0LL) {
        break;
      } else {
        assert(!queues_[sink][dst].empty());
        queues_[sink][dst].pop();
      }
    }
  }
}

void TransportationSuccessiveShortestPath::updateDestQueues(int sink, int src) {
  if (pb_.allocations_[sink][src] != 0LL) {
    return;
  }
  for (int dst = 0; dst < pb_.nbSinks(); ++dst) {
    if (dst == sink) continue;
    assert(!queues_[sink][dst].empty());
    CostType cost = pb_.movingCost(src, sink, dst);
    queues_[sink][dst].emplace(cost, src);
  }
}
}  // namespace coloquinte