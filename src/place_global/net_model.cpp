#include "place_global/net_model.hpp"

#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/IterativeLinearSolvers>
#include <eigen3/Eigen/Sparse>
#include <limits>

NetModel NetModel::xTopology(const Circuit &circuit) {
  return fromData(circuit.cellWidths, circuit.cellX, circuit.cellFixed,
                  circuit.netLimits, circuit.pinCells, circuit.pinXOffsets);
}

NetModel NetModel::yTopology(const Circuit &circuit) {
  return fromData(circuit.cellHeights, circuit.cellY, circuit.cellFixed,
                  circuit.netLimits, circuit.pinCells, circuit.pinYOffsets);
}

NetModel NetModel::fromData(const std::vector<int> &cellSizes,
                            const std::vector<int> &pl,
                            const std::vector<char> &cellFixed,
                            const std::vector<int> &netLimits,
                            const std::vector<int> &pinCells,
                            const std::vector<int> &pinOffsets) {
  NetModel ret(cellSizes.size());
  for (int i = 0; i + 1 < netLimits.size(); ++i) {
    float minPos = std::numeric_limits<float>::infinity();
    float maxPos = -std::numeric_limits<float>::infinity();
    int b = netLimits[i];
    int e = netLimits[i + 1];
    std::vector<int> cells;
    std::vector<float> offsets;
    for (int j = b; j < e; ++j) {
      int cell = pinCells[j];
      int offset = pinOffsets[j];
      if (cellFixed[cell]) {
        int pos = pl[cell] + offset;
        minPos = std::min(minPos, (float)pos);
        maxPos = std::max(maxPos, (float)pos);
      } else {
        cells.push_back(cell);
        // Offset to center of cell instead of lower-left
        offsets.push_back(offset - 0.5f * cellSizes[cell]);
      }
    }
    ret.addNet(cells, offsets, minPos, maxPos);
  }
  ret.check();
  return ret;
}

NetModel::NetModel(int nbCells) : nbCells_(nbCells) { netLimits_.push_back(0); }

void NetModel::addNet(const std::vector<int> &cells) {
  std::vector<float> offsets(cells.size(), 0.0f);
  addNet(cells, offsets);
}

void NetModel::addNet(const std::vector<int> &cells,
                      const std::vector<float> &pinOffsets, float weight) {
  assert(cells.size() == pinOffsets.size());
  if (cells.size() <= 1) {
    return;
  }
  netLimits_.push_back(netLimits_.back() + cells.size());
  netWeight_.push_back(weight);
  netCells_.insert(netCells_.end(), cells.begin(), cells.end());
  netPinOffsets_.insert(netPinOffsets_.end(), pinOffsets.begin(),
                        pinOffsets.end());
}

void NetModel::addNet(const std::vector<int> &cells,
                      const std::vector<float> &pinOffsets, float minPin,
                      float maxPin, float weight) {
  if (cells.empty()) {
    return;
  }
  if (std::isfinite(minPin)) {
    std::vector<int> fCells = cells;
    std::vector<float> fPinOffsets = pinOffsets;
    fCells.push_back(-1);
    fPinOffsets.push_back(minPin);
    if (maxPin != minPin) {
      fCells.push_back(-1);
      fPinOffsets.push_back(maxPin);
    }
    addNet(fCells, fPinOffsets, weight);
  } else {
    addNet(cells, pinOffsets, weight);
  }
}

void NetModel::check() const {
  assert(netCells_.size() == netPinOffsets_.size());
  assert(netLimits_.size() == nbNets() + 1);
  assert(netWeight_.size() == nbNets());
  for (int c : netCells_) {
    assert(c >= -1);
    assert(c < nbCells());
  }
  for (int i = 0; i < nbNets(); ++i) {
    // At least two cells per net
    assert(netLimits_[i] + 2 <= netLimits_[i + 1]);
  }
  assert(netLimits_.front() == 0);
  assert(netLimits_.back() == netCells_.size());
}

std::vector<float> NetModel::pinPositions(const std::vector<float> &pl) const {
  std::vector<float> ret(netCells_.size());
  for (int i = 0; i < netCells_.size(); ++i) {
    ret[i] = pl[netCells_[i]] + netPinOffsets_[i];
  }
  return ret;
}

float NetModel::value(const std::vector<float> &pl) const {
  std::vector<float> pos = pinPositions(pl);
  float ret = 0.0f;
  for (int i = 0; i < nbNets(); ++i) {
    int b = netLimits_[i];
    int e = netLimits_[i + 1];
    float mn = std::numeric_limits<float>::infinity();
    float mx = -std::numeric_limits<float>::infinity();
    for (int j = b; j < e; ++j) {
      mn = std::min(mn, pos[j]);
      mx = std::max(mx, pos[j]);
    }
    ret += (mx - mn);
  }
  return ret;
}

std::tuple<int, int, float, float> NetModel::minPin(
    int net, const std::vector<float> &pl) const {
  int bestI = 0;
  int bestC = -1;
  float bestO = std::numeric_limits<float>::infinity();
  float bestPos = bestO;
  for (int i = 0; i < nbPins(net); ++i) {
    float pos = pinPosition(net, i, pl);
    if (pos < bestPos) {
      bestI = i;
      bestC = pinCell(net, i);
      bestO = pinOffset(net, i);
      bestPos = pos;
    }
  }
  return std::make_tuple(bestI, bestC, bestO, bestPos);
}

std::tuple<int, int, float, float> NetModel::maxPin(
    int net, const std::vector<float> &pl) const {
  int bestI = 0;
  int bestC = -1;
  float bestO = -std::numeric_limits<float>::infinity();
  float bestPos = bestO;
  for (int i = 0; i < nbPins(net); ++i) {
    float pos = pinPosition(net, i, pl);
    if (pos > bestPos) {
      bestI = i;
      bestC = pinCell(net, i);
      bestO = pinOffset(net, i);
      bestPos = pos;
    }
  }
  return std::make_tuple(bestI, bestC, bestO, bestPos);
}

/**
 * Matrix creation for quadratic approximations
 */
class MatrixCreator {
 public:
  explicit MatrixCreator(const NetModel &topo)
      : topo_(topo),
        nbCells_(topo.nbCells()),
        nbSupps_(0),
        rhs_(topo.nbCells()),
        initial_(topo.nbCells()) {}

  int nbCells() const { return nbCells_; }
  int matSize() const { return nbCells_ + nbSupps_; }

  const std::vector<Eigen::Triplet<float> > &mat() const { return mat_; }
  const std::vector<float> &rhs() const { return rhs_; }
  const std::vector<float> &initial() const { return initial_; }

  int addCell(float initialPos);

  std::vector<float> solve();

  static MatrixCreator createStar(const NetModel &topo);
  static MatrixCreator createStar(const NetModel &topo,
                                  const std::vector<float> &pl, float epsilon);
  static MatrixCreator createB2B(const NetModel &topo,
                                 const std::vector<float> &pl, float epsilon);

  void addPenalty(const std::vector<float> &netPlacement,
                  const std::vector<float> &placementTarget,
                  const std::vector<float> &penaltyStrength,
                  float cutoffDistance);

  void addBipoint(int net);
  void addClique(int net);
  void addStar(int net);

  void addBipoint(int net, const std::vector<float> &pl, float epsilon);
  void addClique(int net, const std::vector<float> &pl, float epsilon);
  void addStar(int net, const std::vector<float> &pl, float epsilon);
  void addB2B(int net, const std::vector<float> &pl, float epsilon);

  void check() const;

 private:
  void addPin(int c1, int c2, float offs1, float offs2, float weight);
  void addMovingPin(int c1, int c2, float offs1, float offs2, float weight);
  void addFixedPin(int c1, float offs1, float pos, float weight);

 private:
  const NetModel &topo_;
  int nbCells_;
  int nbSupps_;
  std::vector<Eigen::Triplet<float> > mat_;
  std::vector<float> rhs_;
  std::vector<float> initial_;
};

void MatrixCreator::addMovingPin(int c1, int c2, float offs1, float offs2,
                                 float weight) {
  if (c1 == c2) return;
  assert(c1 >= 0);
  assert(c2 >= 0);
  mat_.emplace_back(c1, c2, -weight);
  mat_.emplace_back(c2, c1, -weight);
  mat_.emplace_back(c1, c1, weight);
  mat_.emplace_back(c2, c2, weight);
  rhs_[c1] += weight * (offs2 - offs1);
  rhs_[c2] += weight * (offs1 - offs2);
}

void MatrixCreator::addFixedPin(int c1, float offs1, float pos, float weight) {
  assert(c1 >= 0);
  mat_.emplace_back(c1, c1, weight);
  rhs_[c1] += weight * (pos - offs1);
}

void MatrixCreator::addPin(int c1, int c2, float offs1, float offs2,
                           float weight) {
  if (c1 == c2) return;
  if (c1 == -1) {
    addFixedPin(c2, offs2, offs1, weight);
    return;
  }
  if (c2 == -1) {
    addFixedPin(c1, offs1, offs2, weight);
    return;
  }
  addMovingPin(c1, c2, offs1, offs2, weight);
}

void MatrixCreator::addPenalty(const std::vector<float> &netPlacement,
                               const std::vector<float> &placementTarget,
                               const std::vector<float> &penaltyStrength,
                               float cutoffDistance) {
  assert(netPlacement.size() == nbCells());
  assert(placementTarget.size() == nbCells());
  assert(penaltyStrength.size() == nbCells());
  for (int i = 0; i < nbCells_; ++i) {
    float dist = std::abs(netPlacement[i] - placementTarget[i]);
    float strength = penaltyStrength[i] / std::max(dist, cutoffDistance);
    addFixedPin(i, 0.0f, placementTarget[i], strength);
  }
}

int MatrixCreator::addCell(float initialPos) {
  initial_.push_back(initialPos);
  rhs_.push_back(0.0f);
  int ret = nbCells_ + nbSupps_;
  nbSupps_++;
  return ret;
}

MatrixCreator MatrixCreator::createStar(const NetModel &topo) {
  MatrixCreator ret(topo);
  for (int i = 0; i < topo.nbNets(); ++i) {
    ret.addStar(i);
  }
  return ret;
}

MatrixCreator MatrixCreator::createStar(const NetModel &topo,
                                        const std::vector<float> &pl,
                                        float epsilon) {
  MatrixCreator ret(topo);
  for (int i = 0; i < topo.nbNets(); ++i) {
    ret.addStar(i, pl, epsilon);
  }
  return ret;
}

MatrixCreator MatrixCreator::createB2B(const NetModel &topo,
                                       const std::vector<float> &pl,
                                       float epsilon) {
  MatrixCreator ret(topo);
  for (int i = 0; i < topo.nbNets(); ++i) {
    ret.addB2B(i, pl, epsilon);
  }
  return ret;
}

void MatrixCreator::addBipoint(int net) {
  addPin(topo_.pinCell(net, 0), topo_.pinCell(net, 1), topo_.pinOffset(net, 0),
         topo_.pinOffset(net, 1), topo_.netWeight(net));
}

void MatrixCreator::addClique(int net) {
  int nb = topo_.nbPins(net);
  float w = 2.0f * topo_.netWeight(net) / (nb * (nb - 1));
  for (int i = 0; i + 1 < nb; ++i) {
    for (int j = i + 1; j < nb; ++j) {
      addPin(topo_.pinCell(net, i), topo_.pinCell(net, j),
             topo_.pinOffset(net, i), topo_.pinOffset(net, j), w);
    }
  }
}

void MatrixCreator::addStar(int net) {
  int nb = topo_.nbPins(net);
  if (nb <= 2) {
    addBipoint(net);
  } else {
    float w = topo_.netWeight(net) / nb;
    int c = addCell(0.0f);
    for (int i = 0; i < nb; ++i) {
      addPin(topo_.pinCell(net, i), c, topo_.pinOffset(net, i), 0.0f, w);
    }
  }
}

void MatrixCreator::addBipoint(int net, const std::vector<float> &pl,
                               float epsilon) {
  float w = topo_.netWeight(net) /
            std::max(epsilon, std::abs(topo_.pinPosition(net, 0, pl) -
                                       topo_.pinPosition(net, 1, pl)));
  addPin(topo_.pinCell(net, 0), topo_.pinCell(net, 1), topo_.pinOffset(net, 0),
         topo_.pinOffset(net, 1), w);
}

void MatrixCreator::addClique(int net, const std::vector<float> &pl,
                              float epsilon) {
  int nb = topo_.nbPins(net);
  float w = 2.0f * topo_.netWeight(net) / (nb * (nb - 1));
  for (int i = 0; i + 1 < nb; ++i) {
    for (int j = i + 1; j < nb; ++j) {
      float distW =
          w / std::max(epsilon, std::abs(topo_.pinPosition(net, i, pl) -
                                         topo_.pinPosition(net, j, pl)));
      addPin(topo_.pinCell(net, i), topo_.pinCell(net, j),
             topo_.pinOffset(net, i), topo_.pinOffset(net, j), distW);
    }
  }
}

void MatrixCreator::addStar(int net, const std::vector<float> &pl,
                            float epsilon) {
  if (topo_.nbPins(net) <= 2) {
    addBipoint(net, pl, epsilon);
  } else {
    auto [minI, minCell, minOffset, minPos] = topo_.minPin(net, pl);
    auto [maxI, maxCell, maxOffset, maxPos] = topo_.maxPin(net, pl);
    float starPos = 0.5f * (minPos + maxPos);
    int starC = addCell(starPos);
    for (int i = 0; i < topo_.nbPins(net); ++i) {
      int c = topo_.pinCell(net, i);
      float pos = topo_.pinPosition(net, i, pl);
      if (i == minI || i == maxI) {
        // Strength of 1 when at the current position
        float w =
            topo_.netWeight(net) / std::max(epsilon, std::abs(pos - starPos));
        addPin(c, starC, topo_.pinOffset(net, i), 0.0f, w);
      } else {
        // Strength of 0 when at the current position, 1 when pushed to a
        // boundary of the net
        float dist = std::min(maxPos - pos, pos - minPos);
        float w = topo_.netWeight(net) / std::max(epsilon, dist);
        addPin(c, starC, topo_.pinOffset(net, i), pos - starPos, w);
      }
    }
  }
}

void MatrixCreator::addB2B(int net, const std::vector<float> &pl,
                           float epsilon) {
  auto [minI, minCell, minOffset, minPos] = topo_.minPin(net, pl);
  auto [maxI, maxCell, maxOffset, maxPos] = topo_.maxPin(net, pl);
  float w = topo_.netWeight(net) / (topo_.nbPins(net) - 1);
  for (int i = 0; i < topo_.nbPins(net); ++i) {
    float pos = topo_.pinPosition(net, i, pl);
    int c = topo_.pinCell(net, i);
    if (i == minI) continue;
    float distMin = w / std::max(epsilon, std::abs(pos - minPos));
    addPin(c, minCell, topo_.pinOffset(net, i), minOffset, distMin);
    if (i == maxI) continue;
    float distMax = w / std::max(epsilon, std::abs(pos - maxPos));
    addPin(c, maxCell, topo_.pinOffset(net, i), maxOffset, distMax);
  }
}

std::vector<float> MatrixCreator::solve() {
  check();
  Eigen::SparseMatrix<float> mat(matSize(), matSize());
  mat.setFromTriplets(mat_.begin(), mat_.end());
  Eigen::Map<Eigen::Matrix<float, -1, 1> > rhs(rhs_.data(), rhs_.size());
  Eigen::Map<Eigen::Matrix<float, -1, 1> > initial(initial_.data(),
                                                   initial_.size());
  Eigen::ConjugateGradient<Eigen::SparseMatrix<float> > solver;
  solver.compute(mat);
  solver.setTolerance(1.0e-6f);
  solver.setMaxIterations(1000);
  Eigen::Matrix<float, -1, 1> res = solver.solveWithGuess(rhs, initial);
  // Copy to a std::vector and remove the fake cells
  std::vector<float> ret;
  ret.resize(matSize());
  Eigen::Matrix<float, -1, 1>::Map(ret.data(), ret.size()) = res;
  ret.resize(nbCells_);
  return ret;
}

void MatrixCreator::check() const {
  assert(matSize() == initial_.size());
  assert(matSize() == rhs_.size());
}

std::vector<float> NetModel::solveStar() const {
  MatrixCreator builder = MatrixCreator::createStar(*this);
  return builder.solve();
}

std::vector<float> NetModel::solveStar(const std::vector<float> &placement,
                                       float epsilon) const {
  MatrixCreator builder = MatrixCreator::createStar(*this, placement, epsilon);
  return builder.solve();
}

std::vector<float> NetModel::solveStar(
    const std::vector<float> &netPlacement, float epsilon,
    const std::vector<float> &placementTarget,
    const std::vector<float> &penaltyStrength, float cutoffDistance) const {
  MatrixCreator builder =
      MatrixCreator::createStar(*this, netPlacement, epsilon);
  builder.addPenalty(netPlacement, placementTarget, penaltyStrength,
                     cutoffDistance);
  return builder.solve();
}

std::vector<float> NetModel::solveB2B(const std::vector<float> &placement,
                                      float epsilon) const {
  MatrixCreator builder = MatrixCreator::createB2B(*this, placement, epsilon);
  return builder.solve();
}

std::vector<float> NetModel::solveB2B(const std::vector<float> &netPlacement,
                                      float epsilon,
                                      const std::vector<float> &placementTarget,
                                      const std::vector<float> &penaltyStrength,
                                      float cutoffDistance) const {
  MatrixCreator builder =
      MatrixCreator::createB2B(*this, netPlacement, epsilon);
  builder.addPenalty(netPlacement, placementTarget, penaltyStrength,
                     cutoffDistance);
  return builder.solve();
}
