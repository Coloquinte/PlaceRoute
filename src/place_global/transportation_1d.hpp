
#include <iostream>
#include <queue>
#include <tuple>
#include <vector>

class Transportation1d;
class Transportation1dSorter;
class Transportation1dSolver;

/**
 * @brief A transportation solver for unidimensional optimal transport
 */
class Transportation1d {
 public:
  using Solution = std::vector<std::tuple<int, int, long long>>;

  /**
   * @brief Initialize the datastructure
   *
   * @param u Source positions
   * @param v Sink positions
   * @param s Source supplies
   * @param d Sink demands
   */
  Transportation1d(const std::vector<long long> &u,
                   const std::vector<long long> &v,
                   const std::vector<long long> &s,
                   const std::vector<long long> &d);
  /**
   * @brief Initialize the datastructure
   *
   * @param u Source positions
   * @param v Sink positions
   * @param s Source supplies
   * @param d Sink demands
   */
  Transportation1d(std::vector<long long> &&u, std::vector<long long> &&v,
                   std::vector<long long> &&s, std::vector<long long> &&d);

  /**
   * @brief Solve a 1D transportation problem
   */
  Solution solve();

  /**
   * @brief Give a rounded solution for a 1D transportation problem; exact only
   * for assignment problems
   */
  std::vector<int> assign();

  /**
   * @brief Number of sources in the problem
   */
  int nbSources() const { return u.size(); }

  /**
   * @brief Number of sinks in the problem
   */
  int nbSinks() const { return v.size(); }

  /**
   * @brief Total source supply for the problem
   */
  long long totalSupply() const;

  /**
   * @brief Total sink demand for the problem
   */
  long long totalDemand() const;

  /**
   * @brief Add more demand if there is not enough
   */
  void balanceDemand();

  /**
   * @brief Compute the cost of allocating a source to a sink
   */
  long long cost(int i, int j) const { return std::abs(u[i] - v[j]); }

  /**
   * @brief Position of the sources
   */
  const std::vector<long long> &sourcePosition() const { return u; }

  /**
   * @brief Position of the sinks
   */
  const std::vector<long long> &sinkPosition() const { return v; }

  /**
   * @brief Supply of the sources
   */
  const std::vector<long long> &sourceSupply() const { return s; }

  /**
   * @brief Demand of the sinks
   */
  const std::vector<long long> &sinkDemand() const { return d; }

  /**
   * @brief Read a serialized problem
   */
  static Transportation1d read(std::istream &f);

  /**
   * @brief Read a serialized solution
   */
  static Solution readSolution(std::istream &f);

  /**
   * @brief Write a serialized problem
   */
  void write(std::ostream &f) const;

  /**
   * @brief Write a serialized solution
   */
  static void writeSolution(const Solution &sol, std::ostream &f);

  /**
   * @brief Write a serialized assignment
   */
  static void writeAssignment(const std::vector<int> &sol, std::ostream &f);

  /**
   * @brief Check that a solution is valid
   */
  void checkSolutionValid(const Solution &sol) const;

  /**
   * @brief Compute the cost for a solution
   */
  long long cost(const Solution &sol) const;

  /**
   * @brief Basic checking
   */
  void check() const;

  /**
   * @brief Check that the problem is sorted
   */
  void checkSorted() const;

  /**
   * @brief Check that the problem is strictly sorted
   */
  void checkStrictlySorted() const;

  /**
   * @brief Check that all capacities are non-zero
   */
  void checkNonZeroCapacities() const;

 protected:
  std::vector<long long> u;
  std::vector<long long> v;
  std::vector<long long> s;
  std::vector<long long> d;

  friend class Transportation1dSorter;
};

class Transportation1dSorter {
 public:
  using Solution = Transportation1d::Solution;

  /**
   * @brief Initialize the datastructure
   */
  Transportation1dSorter(const std::vector<long long> &u,
                         const std::vector<long long> &v,
                         const std::vector<long long> &s,
                         const std::vector<long long> &d);

  Transportation1dSolver convert(const Transportation1d &pb) const;

  /**
   * Convert back a solution given by the solver through preprocessing
   */
  Solution convertSolutionBack(const Solution &sol) const;

  /**
   * Convert back an assignment given by the solver through preprocessing
   */
  std::vector<int> convertAssignmentBack(const std::vector<int> &a) const;

 private:
  std::vector<int> srcOrder;
  std::vector<int> snkOrder;
};

/**
 * @brief Implementation of the transportation solver for unidimensional optimal
 * transport
 */
class Transportation1dSolver : public Transportation1d {
 public:
  using Event = std::pair<long long, long long>;
  using PrioQueue = std::priority_queue<Event>;
  /**
   * @brief Initialize the datastructure
   *
   * @param u Source positions
   * @param v Sink positions
   * @param s Source supplies
   * @param d Sink demands
   */
  Transportation1dSolver(std::vector<long long> &&u, std::vector<long long> &&v,
                         std::vector<long long> &&s,
                         std::vector<long long> &&d);

  /**
   * @brief Compute the total supply at the sources
   */
  long long totalSupply() const { return S.back(); }

  /**
   * @brief Compute the total demand at the sinks
   */
  long long totalDemand() const { return D.back(); };

  /**
   * @brief Compute the change in reduced cost at a boundary
   */
  long long delta(int i, int j) const {
    return cost(i, j + 1) + cost(i + 1, j) - cost(i + 1, j + 1) - cost(i, j);
  }

  /**
   * @brief Run the whole optimization with efficient algorithm
   */
  void run();

  /**
   * @brief Obtain the current solution
   */
  Solution computeSolution() const;

  /**
   * @brief Obtain an assignment close to the fractional solution
   */
  std::vector<int> computeAssignment() const;

  /**
   * @brief Check that the given solution is optimal
   */
  void checkSolutionOptimal(const Solution &sol) const;

  /**
   * @brief Check the datastructure
   */
  void check() const;

 private:
  /**
   * @brief Initialize the additional datastructures
   */
  void setupData();

  /**
   * @brief Flush the positions according to the non-overlap constraints
   */
  void flushPositions();

  /**
   * @brief Push a single source
   */
  void push(int i);

  /**
   * @brief Push a single source
   */
  void pushOnce(int i);

  /**
   * @brief Push to the last sink to be used
   */
  void pushToLastSink(int i);

  /**
   * @brief Push to the next sink to be used
   */
  void pushToNewSink(int i);

  /**
   * @brief Push the events corresponding to a given source
   */
  void pushNewSourceEvents(int i);

  /**
   * @brief Push the events corresponding to entering a new sink
   */
  void pushNewSinkEvents(int i, int j);

  /**
   * @brief Get the reduced cost at the current position
   */
  long long getSlope(bool pop = false);

  /**
   * @brief Update the optimal sink
   */
  void updateOptimalSink(int i);

 private:
  // Additional data
  std::vector<long long> S;
  std::vector<long long> D;

  // Solution and state
  std::vector<long long> p;
  PrioQueue events;
  long long lastPosition;
  int lastOccupiedSink;
  int optimalSink;
};
