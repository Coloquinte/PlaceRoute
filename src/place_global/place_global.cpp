#include "place_global.hpp"

#include <cassert>
#include <chrono>
#include <cmath>
#include <future>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <utility>

#include "density_legalizer.hpp"
#include "net_model.hpp"

namespace coloquinte {

namespace {
/**
 * Mix two placements with a weight: 0 for the first, 1 for the second, or
 * in-between
 */
std::vector<float> blendPlacement(const std::vector<float> &v1,
                                  const std::vector<float> &v2,
                                  float blending) {
  if (blending == 0.0f) {
    return v1;
  }
  if (blending == 1.0f) {
    return v2;
  }
  std::vector<float> ret;
  assert(v1.size() == v2.size());
  ret.reserve(v1.size());
  for (int i = 0; i < v1.size(); ++i) {
    ret.push_back((1.0f - blending) * v1[i] + blending * v2[i]);
  }
  return ret;
}
}  // namespace

void GlobalPlacer::place(Circuit &circuit, const ColoquinteParameters &params,
                         const std::optional<PlacementCallback> &callback) {
  params.check();
  std::cout << "Global placement starting" << std::endl;
  auto startTime = std::chrono::steady_clock::now();
  GlobalPlacer pl(circuit, params);
  pl.callback_ = callback;
  pl.run();
  auto endTime = std::chrono::steady_clock::now();
  std::chrono::duration<float> duration = endTime - startTime;
  std::cout << std::fixed << std::setprecision(2) << "Global placement done in "
            << duration.count() << "s" << std::endl;
  pl.exportPlacement(circuit);
}

GlobalPlacer::GlobalPlacer(Circuit &circuit, const ColoquinteParameters &params)
    : leg_(DensityLegalizer::fromIspdCircuit(
          circuit, params.global.roughLegalization.binSize,
          params.global.roughLegalization.sideMargin)),
      xtopo_(NetModel::xTopology(circuit)),
      ytopo_(NetModel::yTopology(circuit)),
      params_(params),
      circuit_(circuit) {
  rgen_.seed(params_.seed);
  averageCellLength_ = computeAverageCellSize();
  perCellPenalty_ = computePerCellPenalty();
  auto rlp = params.global.roughLegalization;
  LegalizationModel m = rlp.costModel;
  DensityLegalizer::Parameters legParams;
  legParams.nbSteps = rlp.nbSteps;
  legParams.costModel = rlp.costModel;
  legParams.lineReoptSize = rlp.lineReoptSize;
  legParams.lineReoptOverlap = rlp.lineReoptOverlap;
  legParams.diagReoptSize = rlp.diagReoptSize;
  legParams.diagReoptOverlap = rlp.diagReoptOverlap;
  legParams.squareReoptSize = rlp.squareReoptSize;
  legParams.squareReoptOverlap = rlp.squareReoptOverlap;
  legParams.unidimensionalTransport =
      rlp.unidimensionalTransport && m == LegalizationModel::L1;
  legParams.coarseningLimit = rlp.coarseningLimit;
  if (m == LegalizationModel::L1 || m == LegalizationModel::L2 ||
      m == LegalizationModel::LInf) {
    float dist = leg_.placementArea().width() + leg_.placementArea().height();
    legParams.quadraticPenaltyFactor = rlp.quadraticPenalty / dist;
  }
  leg_.setParams(legParams);
}

void GlobalPlacer::exportPlacement(Circuit &circuit) const {
  assert(xtopo_.nbCells() == circuit.nbCells());
  assert(ytopo_.nbCells() == circuit.nbCells());
  assert(leg_.nbCells() == circuit.nbCells());
  float w = params_.global.exportBlending;
  std::vector<float> xplace = blendPlacement(xPlacementLB_, xPlacementUB_, w);
  std::vector<float> yplace = blendPlacement(yPlacementLB_, yPlacementUB_, w);
  exportPlacement(circuit, xplace, yplace);
}

void GlobalPlacer::exportPlacement(Circuit &circuit,
                                   const std::vector<float> &xplace,
                                   const std::vector<float> &yplace) {
  assert(xplace.size() == circuit.nbCells());
  assert(yplace.size() == circuit.nbCells());

  for (int i = 0; i < circuit.nbCells(); ++i) {
    if (circuit.isFixed(i)) {
      continue;
    }
    circuit.cellX_[i] = std::round(xplace[i] - 0.5 * circuit.placedWidth(i));
    circuit.cellY_[i] = std::round(yplace[i] - 0.5 * circuit.placedHeight(i));
  }
}

float GlobalPlacer::computeAverageCellSize() const {
  float totalDemand = leg_.totalDemand();
  float avgDemand = totalDemand == 0.0 ? 0.0 : totalDemand / leg_.nbCells();
  return std::sqrt(avgDemand);
}

std::vector<float> GlobalPlacer::computePerCellPenalty() const {
  int nbCells = leg_.nbNonEmptyCells();
  float meanArea = leg_.totalDemand() / std::max(1, nbCells);
  std::vector<float> ret;
  ret.reserve(leg_.nbCells());

  for (int i = 0; i < leg_.nbCells(); ++i) {
    ret.push_back(std::pow(leg_.cellDemand(i) / meanArea,
                           params_.global.penalty.areaExponent));
  }
  return ret;
}

std::vector<float> GlobalPlacer::computeIterationPerCellPenalty() {
  std::vector<float> penalty = perCellPenalty_;
  for (float &s : penalty) {
    s *= penalty_;
  }
  if (params_.global.noise > 0.0) {
    std::uniform_real_distribution<float> dist(0.0f, params_.global.noise);
    for (float &s : penalty) {
      s *= (1.0f + dist(rgen_));
    }
  }
  return penalty;
}

void GlobalPlacer::run() {
  runInitialLB();
  penalty_ = params_.global.penalty.initialValue;
  approximationDistance_ = initialApproximationDistance();
  penaltyCutoffDistance_ = initialPenaltyCutoffDistance();

  float lb = valueLB();
  float ub = std::numeric_limits<float>::infinity();
  int firstStep = params_.global.nbInitialSteps + 1;
  for (step_ = firstStep; step_ <= params_.global.maxNbSteps; ++step_) {
    std::cout << "#" << step_ << ":" << std::flush;
    runUB();
    ub = valueUB();
    std::cout << std::defaultfloat << std::setprecision(4) << "\tUB " << ub;

    float dist = leg_.meanDistance();
    std::cout << std::fixed << std::setprecision(1) << "\tDist " << dist;
    std::cout << std::flush;

    float gap = (ub - lb) / ub;
    // Stop if distance or the difference between LB and UB is small enough
    if (gap < params_.global.gapTolerance || dist < distanceTolerance()) {
      std::cout << std::endl;
      break;
    }
    for (int i = 0; i < params_.global.nbStepsBeforeRoughLegalization; ++i) {
      if (i != 0) {
        std::cout << "#" << step_ << ":\t........\t........";
        std::cout << std::flush;
      }
      runLB();
      lb = valueLB();
      std::cout << std::defaultfloat << std::setprecision(4) << "\tLB " << lb
                << std::endl;
    }
    penalty_ *= params_.global.penalty.updateFactor;
    penaltyCutoffDistance_ *= params_.global.penalty.cutoffDistanceUpdateFactor;
    approximationDistance_ *=
        params_.global.continuousModel.approximationDistanceUpdateFactor;
  }
  runUB();
}

float GlobalPlacer::valueLB() const {
  return xtopo_.value(xPlacementLB_) + ytopo_.value(yPlacementLB_);
}

float GlobalPlacer::valueUB() const {
  return xtopo_.value(xPlacementUB_) + ytopo_.value(yPlacementUB_);
}

void GlobalPlacer::runInitialLB() {
  NetModel::Parameters params;
  params.netModel = params_.global.continuousModel.netModel;
  params.tolerance =
      params_.global.continuousModel.conjugateGradientErrorTolerance;
  params.maxNbIterations =
      params_.global.continuousModel.maxNbConjugateGradientSteps;
  xPlacementLB_ = xtopo_.solveStar(params);
  yPlacementLB_ = ytopo_.solveStar(params);
  std::cout << std::defaultfloat << std::setprecision(4) << "#0:\tLB "
            << valueLB() << std::endl;
  callback(PlacementStep::LowerBound, xPlacementLB_, yPlacementLB_);
  for (step_ = 1; step_ <= params_.global.nbInitialSteps; ++step_) {
    xPlacementLB_ = xtopo_.solve(xPlacementLB_, params);
    yPlacementLB_ = ytopo_.solve(yPlacementLB_, params);
    std::cout << std::defaultfloat << std::setprecision(4) << "#" << step_
              << ":\tLB " << valueLB() << std::endl;
    callback(PlacementStep::LowerBound, xPlacementLB_, yPlacementLB_);
  }
  // Simplify blending solutions by having a UB immediately
  xPlacementUB_ = xPlacementLB_;
  yPlacementUB_ = yPlacementLB_;
}

void GlobalPlacer::runLB() {
  // Compute the parameters for the continuous model solver
  NetModel::Parameters params;
  params.netModel = params_.global.continuousModel.netModel;
  params.approximationDistance = approximationDistance_;
  params.penaltyCutoffDistance = penaltyCutoffDistance_;
  params.tolerance =
      params_.global.continuousModel.conjugateGradientErrorTolerance;
  params.maxNbIterations =
      params_.global.continuousModel.maxNbConjugateGradientSteps;

  // Compute the per-cell penalty with randomization
  std::vector<float> penalty = computeIterationPerCellPenalty();

  float w = params_.global.penalty.targetBlending;
  std::vector<float> xTarget = blendPlacement(xPlacementLB_, xPlacementUB_, w);
  std::vector<float> yTarget = blendPlacement(yPlacementLB_, yPlacementUB_, w);

  // Solve the continuous model (x and y independently)
  std::future<std::vector<float> > x =
      std::async(std::launch::async, &NetModel::solveWithPenalty, &xtopo_,
                 xPlacementLB_, xTarget, penalty, params);
  std::future<std::vector<float> > y =
      std::async(std::launch::async, &NetModel::solveWithPenalty, &ytopo_,
                 yPlacementLB_, yTarget, penalty, params);
  xPlacementLB_ = x.get();
  yPlacementLB_ = y.get();
  callback(PlacementStep::LowerBound, xPlacementLB_, yPlacementLB_);
}

void GlobalPlacer::runUB() {
  float w = params_.global.roughLegalization.targetBlending;
  std::vector<float> xTarget = blendPlacement(xPlacementLB_, xPlacementUB_, w);
  std::vector<float> yTarget = blendPlacement(yPlacementLB_, yPlacementUB_, w);
  leg_.updateCellTargetX(xTarget);
  leg_.updateCellTargetY(yTarget);
  leg_.run();
  xPlacementUB_ = leg_.spreadCoordX(xTarget);
  yPlacementUB_ = leg_.spreadCoordY(yTarget);
  callback(PlacementStep::UpperBound, xPlacementUB_, yPlacementUB_);
}

void GlobalPlacer::callback(PlacementStep step,
                            const std::vector<float> &xplace,
                            const std::vector<float> &yplace) {
  if (!callback_.has_value()) return;
  exportPlacement(circuit_, xplace, yplace);
  callback_.value()(step);
}
}  // namespace coloquinte