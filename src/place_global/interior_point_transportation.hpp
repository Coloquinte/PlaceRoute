
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

  template <int N, int NRM>
  static InteriorPointTransportation<T> createFromDistance(
      const Vector &dem, const Vector &cap, const Position<N> &srcPos,
      const Position<N> &snkPos, T quadratic = 0.0);
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
  template <int N, int NRM>
  static InteriorPointTransportation<T> makeRandomFromDistance(int nbSources,
                                                               int nbSinks,
                                                               int seed = 1);

  int nbSources() const { return demands.rows(); }

  int nbSinks() const { return capacities.rows(); }

  int nbVariables() const { return nbSources() * nbSinks(); }

  int nbConstraints(bool ignoreRedundant = true) const {
    return nbSources() + nbSinks() - ignoreRedundant;
  }

  T totalDemand() const { return demands.sum(); }
  T totalCapacity() const { return capacities.sum(); }
  Vector demandUsage(const Vector &x) const;
  Vector capacityUsage(const Vector &x) const;

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

  Vector solve(int maxIter = 100, T eps = 1.0e-9, T theta = 0.9995) const;

  void check() const;

 private:
  void newtonDirection(const Vector &r_b, const Vector &r_c,
                       const Vector &r_x_s, const Vector &x, const Vector &s,
                       Vector &dx, Vector &dy, Vector &ds) const;
  void stepSize(const Vector &x, const Vector &s, const Vector &d_x,
                const Vector &d_s, T eta, T &alpha_x, T &alpha_s) const;

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
  return Vector::Constant(nbConstraints(), -margin / 2.0);
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
  Vector s = costVector() - applyAt(y);
  if ((s.array() < -margin).any()) {
    throw std::runtime_error("Invalid dual solution: negative reduced costs");
  }
}

template <class T>
Eigen::Matrix<T, Eigen::Dynamic, 1> InteriorPointTransportation<T>::costVector()
    const {
  return costs.transpose().reshaped();
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
Eigen::Matrix<T, Eigen::Dynamic, 1> InteriorPointTransportation<T>::demandUsage(
    const Vector &x) const {
  auto view = x.reshaped(nbSinks(), nbSources());
  return view.colwise().sum();
}

template <class T>
Eigen::Matrix<T, Eigen::Dynamic, 1>
InteriorPointTransportation<T>::capacityUsage(const Vector &x) const {
  auto view = x.reshaped(nbSinks(), nbSources());
  return view.rowwise().sum();
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
  xt = Dt.lu().solve(xt);

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
Eigen::Matrix<T, Eigen::Dynamic, 1> InteriorPointTransportation<T>::solve(
    int maxIter, T eps, T theta) const {
  Vector c = costVector();
  Vector b = constraintVector();

  Vector x = initialPrimalSolution();
  Vector y = initialDualSolution(1.0);
  Vector s = c - applyAt(y);

  T bc = 1 + std::max(b.norm(), c.norm());

  for (int niter = 0; niter < maxIter; ++niter) {
    // Compute residuals and update mu
    Vector r_b = applyA(x) - b;
    Vector r_c = applyAt(y) + s - c;
    Vector r_x_s = x.cwiseProduct(s);
    T mu = r_x_s.mean();

    // Check relative decrease in residual, for purposes of convergence test
    T sqResidual = r_b.squaredNorm() + r_c.squaredNorm() + r_x_s.squaredNorm();
    T residual = std::sqrt(sqResidual) / bc;
    if (residual < eps) {
      break;
    }

    // ----- Predictor step -----

    // Get affine-scaling direction
    Vector dx_aff, dy_aff, ds_aff;
    newtonDirection(r_b, r_c, r_x_s, x, s, dx_aff, dy_aff, ds_aff);

    // Get affine-scaling step length
    T alpha_x_aff, alpha_s_aff;
    stepSize(x, s, dx_aff, ds_aff, 1, alpha_x_aff, alpha_s_aff);
    T mu_aff = (x + alpha_x_aff * dx_aff).dot(s + alpha_s_aff * ds_aff) /
               nbVariables();

    // Set central parameter
    T sigma = std::pow(mu_aff / mu, 3);

    // ----- Corrector step -----

    // Set up right hand sides
    r_x_s = r_x_s + dx_aff.cwiseProduct(ds_aff) -
            Vector::Constant(nbVariables(), sigma * mu);

    // Get corrector's direction
    Vector dx_cc, dy_cc, ds_cc;
    newtonDirection(r_b, r_c, r_x_s, x, s, dx_cc, dy_cc, ds_cc);

    // Compute search direction and step
    Vector dx = dx_aff + dx_cc;
    Vector dy = dy_aff + dy_cc;
    Vector ds = ds_aff + ds_cc;

    T alpha_x, alpha_s;
    stepSize(x, s, dx, ds, theta, alpha_x, alpha_s);

    // Update iterates
    x = x + alpha_x * dx;
    y = y + alpha_s * dy;
    s = s + alpha_s * ds;
  }

  return x;
}

template <class T>
void InteriorPointTransportation<T>::newtonDirection(
    const Vector &r_b, const Vector &r_c, const Vector &r_x_s, const Vector &x,
    const Vector &s, Vector &dx, Vector &dy, Vector &ds) const {
  // Block cholesky solve for D dx + A dy = u, At dx = v
  Vector u = -r_c + r_x_s.cwiseProduct(x.cwiseInverse());
  Vector v = -r_b;
  Vector D = -s.cwiseProduct(x.cwiseInverse());

  // First triangular solve
  v = v - applyA(u.cwiseProduct(D.cwiseInverse()));

  // Diagonal solve
  u = u.cwiseProduct(D.cwiseInverse());
  dy = -applyAGAtInv(D.cwiseInverse(), v);

  // Second triangular solve
  dx = u - applyAt(dy).cwiseProduct(D.cwiseInverse());

  // Final step
  ds = -(r_x_s + s.cwiseProduct(dx)).cwiseProduct(x.cwiseInverse());
}

template <class T>
void InteriorPointTransportation<T>::stepSize(const Vector &x, const Vector &s,
                                              const Vector &d_x,
                                              const Vector &d_s, T eta,
                                              T &alpha_x, T &alpha_s) const {
  alpha_x =
      -1.0 / std::min((d_x.cwiseProduct(x.cwiseInverse())).minCoeff(), (T)-1.0);
  alpha_x = std::min((T)1.0, eta * alpha_x);
  alpha_s =
      -1.0 / std::min((d_s.cwiseProduct(s.cwiseInverse())).minCoeff(), (T)-1.0);
  alpha_s = std::min((T)1.0, eta * alpha_s);
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

template <class T>
template <int N, int NRM>
InteriorPointTransportation<T>
InteriorPointTransportation<T>::makeRandomFromDistance(int nbSources,
                                                       int nbSinks, int seed) {
  std::mt19937 rgen(seed);
  std::uniform_real_distribution<T> dist;
  Vector demands(nbSources);
  Vector capacities(nbSinks);
  Position<N> srcPos(nbSources, N);
  Position<N> snkPos(nbSinks, N);
  for (int i = 0; i < nbSources; ++i) {
    demands(i) = dist(rgen);
  }
  for (int j = 0; j < nbSinks; ++j) {
    capacities(j) = dist(rgen);
  }
  for (int i = 0; i < nbSources; ++i) {
    for (int k = 0; k < N; ++k) {
      srcPos(i, k) = dist(rgen);
    }
  }
  for (int j = 0; j < nbSinks; ++j) {
    for (int k = 0; k < N; ++k) {
      snkPos(j, k) = dist(rgen);
    }
  }
  capacities *= (demands.sum() / capacities.sum());
  return InteriorPointTransportation<T>::createFromDistance<N, NRM>(
      demands, capacities, srcPos, snkPos);
}

template <class T>
template <int N, int NRM>
InteriorPointTransportation<T>
InteriorPointTransportation<T>::createFromDistance(const Vector &dem,
                                                   const Vector &cap,
                                                   const Position<N> &srcPos,
                                                   const Position<N> &snkPos,
                                                   T quadratic) {
  if (srcPos.rows() != dem.rows()) {
    throw std::runtime_error(
        "Transportation problem should have as many source positions as "
        "demands.");
  }
  if (snkPos.rows() != cap.rows()) {
    throw std::runtime_error(
        "Transportation problem should have as many sink positions as "
        "capacities.");
  }
  Matrix costs(dem.rows(), cap.rows());
  for (int i = 0; i < dem.rows(); ++i) {
    for (int j = 0; j < cap.rows(); ++j) {
      costs(i, j) = (srcPos.row(i) - snkPos.row(j)).template lpNorm<NRM>();
    }
  }
  // TODO: add quadratic case
  return InteriorPointTransportation<T>(dem, cap, costs);
}

template <class T>
template <int N>
InteriorPointTransportation<T> InteriorPointTransportation<T>::createEuclidean(
    const Vector &dem, const Vector &cap, const Position<N> &srcPos,
    const Position<N> &snkPos, T quadratic) {
  return createFromDistance<T, N, 2>(dem, cap, srcPos, snkPos, quadratic);
}
template <class T>
template <int N>
InteriorPointTransportation<T> InteriorPointTransportation<T>::createManhattan(
    const Vector &dem, const Vector &cap, const Position<N> &srcPos,
    const Position<N> &snkPos, T quadratic) {
  return createFromDistance<T, N, 1>(dem, cap, srcPos, snkPos, quadratic);
}
