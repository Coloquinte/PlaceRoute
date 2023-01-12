
#include <eigen3/Eigen/Cholesky>
#include <eigen3/Eigen/Dense>
#include <random>

template <class T>
class InteriorPointTransportation {
 public:
  using Matrix = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;
  using Vector = Eigen::Matrix<T, Eigen::Dynamic, 1>;
  template <int N>
  using Position = Eigen::Matrix<T, Eigen::Dynamic, N>;

  InteriorPointTransportation(const Vector &dem, const Vector &cap,
                              const Matrix &c);

  template <int N>
  static InteriorPointTransportation<T> createEuclidean(
      const Vector &dem, const Vector &cap, const Position<N> &srcPos,
      const Position<N> &snkPos, T quadratic = 0.0);
  template <int N>
  static InteriorPointTransportation<T> createManhattan(
      const Vector &dem, const Vector &cap, const Position<N> &srcPos,
      const Position<N> &snkPos, T quadratic = 0.0);

  static InteriorPointTransportation<T> makeRandom(int nbSources, int nbSinks,
                                                   int seed = 1);

  int nbSources() const { return demands.rows(); }

  int nbSinks() const { return capacities.rows(); }

  int nbVariables() const { return nbSources() * nbSinks(); }

  int nbConstraints(bool ignoreRedundant = true) const {
    return nbSources() + nbSinks() - ignoreRedundant;
  }

  T totalDemand() const { return demands.sum(); }
  T totalCapacity() const { return capacities.sum(); }

  void normalizeCosts();
  void normalizeDemands();

  Vector initialPrimalSolution() const;
  Vector initialDualSolution(T margin = 1.0e-2) const;
  void checkPrimalSolution(const Vector &x, T margin = 1.0e-6) const;
  void checkDualSolution(const Vector &s, T margin = 1.0e-6) const;

  Vector costVector() const;
  Vector constraintVector() const;

  Vector applyA(const Vector &x) const;
  Vector applyAt(const Vector &x) const;
  Vector applyAGAt(const Vector &g, const Vector &x) const;
  Vector applyAGAtInv(const Vector &g, const Vector &x) const;

  Matrix makeDenseA() const;
  Matrix makeDenseAGAt(Vector g) const;

  void check() const;

 private:
  Vector demands;
  Vector capacities;
  Matrix costs;
};

template <class T>
InteriorPointTransportation<T>::InteriorPointTransportation(const Vector &dem,
                                                            const Vector &cap,
                                                            const Matrix &c) {
  this->demands = dem;
  this->capacities = cap;
  this->costs = c;
  normalizeDemands();
  normalizeCosts();
}

template <class T>
void InteriorPointTransportation<T>::normalizeCosts() {
  for (int i = 0; i < costs.rows(); ++i) {
    costs.row(i) = costs.row(i).array() - costs.row(i).minCoeff();
  }
  costs /= costs.mean();
}

template <class T>
void InteriorPointTransportation<T>::normalizeDemands() {
  T factor = static_cast<T>(1.0) / demands.mean();
  demands *= factor;
  capacities *= factor;
}

template <class T>
void InteriorPointTransportation<T>::check() const {
  if (costs.rows() != nbSources()) {
    throw std::runtime_error(
        "Transportation problem: Cost matrix rows should match the number of "
        "sources");
  }
  if (costs.cols() != nbSinks()) {
    throw std::runtime_error(
        "Transportation problem: Cost matrix columns should match the number "
        "of sinks");
  }
  T dem = totalDemand();
  T cap = totalCapacity();
  if (std::abs(dem - cap) >
      std::max(1.0e-8, 1.0e-6 * (std::abs(dem) + std::abs(cap)))) {
    throw std::runtime_error("Total demand and capacity should be close");
  }
  if ((costs.array() < 0).any()) {
    throw std::runtime_error(
        "Transportation problem: there should be no negative cost");
  }
  if (std::abs(costs.mean() - 1.0) > 1.0e-6) {
    throw std::runtime_error(
        "Transportation problem: cost mean should be close to 1");
  }
  if ((demands.array() <= 0).any()) {
    throw std::runtime_error(
        "Transportation problem: there should be only positive demands");
  }
  if ((capacities.array() <= 0).any()) {
    throw std::runtime_error(
        "Transportation problem: there should be only positive capacities");
  }
}

template <class T>
Eigen::Matrix<T, Eigen::Dynamic, 1>
InteriorPointTransportation<T>::initialPrimalSolution() const {
  Vector proportions = capacities / capacities.sum();
  return (proportions * demands.transpose()).reshaped();
}

template <class T>
Eigen::Matrix<T, Eigen::Dynamic, 1>
InteriorPointTransportation<T>::initialDualSolution(T margin) const {
  return Vector::Constant(nbConstraints(), -margin);
}

template <class T>
void InteriorPointTransportation<T>::checkPrimalSolution(const Vector &x,
                                                         T margin) const {
  if ((x.array() < -margin).any()) {
    throw std::runtime_error("Invalid primal solution: negative variables");
  }
  Vector constraints = applyA(x);
  if (((constraints.array() - constraintVector().array()).abs() > margin)
          .any()) {
    throw std::runtime_error("Invalid primal solution: constraints not met");
  }
}

template <class T>
void InteriorPointTransportation<T>::checkDualSolution(const Vector &y,
                                                       T margin) const {
  Vector s = costs.reshaped() - applyAt(y);
  if ((s.array() < -margin).any()) {
    throw std::runtime_error("Invalid dual solution: negative reduced costs");
  }
}

template <class T>
Eigen::Matrix<T, Eigen::Dynamic, 1> InteriorPointTransportation<T>::costVector()
    const {
  return costs.reshaped();
}

template <class T>
Eigen::Matrix<T, Eigen::Dynamic, 1>
InteriorPointTransportation<T>::constraintVector() const {
  Vector ret(nbConstraints());
  ret.segment(0, nbSources()) = demands;
  ret.segment(nbSources(), nbSinks() - 1) =
      capacities.segment(0, nbSinks() - 1);
  return ret;
}

template <class T>
Eigen::Matrix<T, Eigen::Dynamic, 1> InteriorPointTransportation<T>::applyA(
    const Vector &x) const {
  assert(x.rows() == nbVariables());
  auto view = x.reshaped(nbSinks(), nbSources());
  Vector ret(nbConstraints());
  ret.segment(0, nbSources()) = view.colwise().sum();
  ret.segment(nbSources(), nbSinks() - 1) =
      view.rowwise().sum().segment(0, nbSinks() - 1);
  return ret;
}

template <class T>
Eigen::Matrix<T, Eigen::Dynamic, 1> InteriorPointTransportation<T>::applyAt(
    const Vector &x) const {
  assert(x.rows() == nbConstraints());
  Vector ret = Vector::Zero(nbVariables());
  auto view = ret.reshaped(nbSinks(), nbSources());
  auto xs = x.segment(0, nbSources());
  auto xt = x.segment(nbSources(), nbSinks() - 1);
  view.topRows(nbSinks() - 1).colwise() += xt;
  view.rowwise() += xs.transpose();
  return ret;
}

template <class T>
Eigen::Matrix<T, Eigen::Dynamic, 1> InteriorPointTransportation<T>::applyAGAt(
    const Vector &g, const Vector &x) const {
  assert(g.rows() == nbVariables());
  assert(x.rows() == nbConstraints());
  auto view = g.reshaped(nbSinks(), nbSources());
  Vector gt = view.topRows(nbSinks() - 1).rowwise().sum();
  Vector gs = view.colwise().sum();
  auto G = view.topRows(nbSinks() - 1);
  auto xs = x.segment(0, nbSources());
  auto xt = x.segment(nbSources(), nbSinks() - 1);
  Vector ret(nbConstraints());
  ret.segment(0, nbSources()) = xs.cwiseProduct(gs) + G.transpose() * xt;
  ret.segment(nbSources(), nbSinks() - 1) = xt.cwiseProduct(gt) + G * xs;
  return ret;
}

template <class T>
Eigen::Matrix<T, Eigen::Dynamic, 1>
InteriorPointTransportation<T>::applyAGAtInv(const Vector &g,
                                             const Vector &x) const {
  assert(g.rows() == nbVariables());
  assert(x.rows() == nbConstraints());
  auto view = g.reshaped(nbSinks(), nbSources());
  Vector gt = view.topRows(nbSinks() - 1).rowwise().sum();
  Vector gs = view.colwise().sum();
  auto G = view.topRows(nbSinks() - 1);

  Matrix L = G * gs.cwiseInverse().asDiagonal();
  auto Ds = gs.asDiagonal();
  Matrix Dt = gt.asDiagonal().toDenseMatrix() - (L * Ds * L.transpose()).eval();

  Vector xs = x.segment(0, nbSources());
  Vector xt = x.segment(nbSources(), nbSinks() - 1);

  // First triangular solve
  xt = xt - L * xs;

  // Diagonal solve
  xs = xs.cwiseProduct(gs.cwiseInverse());
  xt = Dt.llt().solve(xt);

  // Second triangular solve
  xs = xs - L.transpose() * xt;

  Vector ret(nbConstraints());
  ret.segment(0, nbSources()) = xs;
  ret.segment(nbSources(), nbSinks() - 1) = xt;
  return ret;
}

template <class T>
Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>
InteriorPointTransportation<T>::makeDenseA() const {
  Matrix a = Matrix::Zero(nbConstraints(), nbVariables());
  for (int i = 0; i < nbSources(); ++i) {
    for (int j = 0; j < nbSinks(); ++j) {
      a(i, i * nbSinks() + j) = 1.0;
    }
  }
  for (int j = 0; j + 1 < nbSinks(); ++j) {
    for (int i = 0; i < nbSources(); ++i) {
      a(nbSources() + j, i * nbSinks() + j) = 1.0;
    }
  }
  return a;
}

template <class T>
Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>
InteriorPointTransportation<T>::makeDenseAGAt(Vector g) const {
  Matrix a = makeDenseA();
  return a * g.asDiagonal() * a.transpose();
}

template <class T>
InteriorPointTransportation<T> InteriorPointTransportation<T>::makeRandom(
    int nbSources, int nbSinks, int seed) {
  std::mt19937 rgen(seed);
  std::uniform_real_distribution<T> dist;
  Vector demands(nbSources);
  Vector capacities(nbSinks);
  Matrix costs(nbSources, nbSinks);
  for (int i = 0; i < nbSources; ++i) {
    demands(i) = dist(rgen);
  }
  for (int j = 0; j < nbSinks; ++j) {
    capacities(j) = dist(rgen);
  }
  for (int i = 0; i < nbSources; ++i) {
    for (int j = 0; j < nbSinks; ++j) {
      costs(i, j) = dist(rgen);
    }
  }
  capacities *= (demands.sum() / capacities.sum());
  return InteriorPointTransportation<T>(demands, capacities, costs);
}
