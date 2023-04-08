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

void GlobalPlacerParameters::check() const {
  if (maxNbSteps < 0) {
    throw std::runtime_error("Invalid number of steps");
  }
  if (nbInitialSteps < 0) {
    throw std::runtime_error("Invalid number of initial steps");
  }
  if (nbInitialSteps >= maxNbSteps) {
    throw std::runtime_error(
        "Number of initial steps should be lower than max number");
  }
  if (nbStepsPerLegalization < 1) {
    throw std::runtime_error(
        "Number of steps per legalization should be positive");
  }
  if (gapTolerance < 0.0f || gapTolerance > 1.0f) {
    throw std::runtime_error(
        "Invalid gap tolerance (should be between 0 and 1)");
  }
  if (distanceTolerance < 0.0f) {
    throw std::runtime_error(
        "Invalid distance tolerance (should be non-negative");
  }
  if (penaltyCutoffDistance < 1.0e-6) {
    throw std::runtime_error("Too small cutoff distance may lead to issues");
  }
  if (penaltyCutoffDistanceUpdateFactor < 0.8 ||
      penaltyCutoffDistanceUpdateFactor > 1.2) {
    throw std::runtime_error(
        "Penalty cutoff update factor should be close to 1");
  }
  if (penaltyAreaExponent < 0.49 || penaltyAreaExponent > 1.01) {
    throw std::runtime_error(
        "Penalty area exponent should be between 0.5 and 1");
  }
  if (initialPenalty <= 0.0f) {
    throw std::runtime_error("Initial penalty should be positive");
  }
  if (penaltyUpdateFactor <= 1.0f || penaltyUpdateFactor >= 2.0f) {
    throw std::runtime_error(
        "Penalty update factor should be between one and two");
  }
  if (penaltyTargetBlending < 0.1f || penaltyTargetBlending > 1.1f) {
    throw std::runtime_error(
        "Penalty target blending should generally be between 0.5 and 1");
  }
  if (approximationDistance < 1.0e-6) {
    throw std::runtime_error(
        "Too small approximation distance may lead to issues");
  }
  if (approximationDistanceUpdateFactor < 0.8 ||
      approximationDistanceUpdateFactor > 1.2) {
    throw std::runtime_error(
        "Approximation distance update factor should be close to 1");
  }
  if (approximationDistance > 1.0e3) {
    throw std::runtime_error(
        "Too large approximation distance is highly imprecise");
  }
  if (maxNbConjugateGradientSteps <= 0) {
    throw std::runtime_error(
        "Must have positive number of steps during conjugate gradients");
  }
  if (conjugateGradientErrorTolerance < 1.0e-8) {
    throw std::runtime_error("Too small error tolerance may lead to issues");
  }
  if (conjugateGradientErrorTolerance > 1.0) {
    throw std::runtime_error("Too large error tolerance is highly imprecise");
  }
  if (roughLegalizationNbSteps < 0) {
    throw std::runtime_error(
        "Must have non-negative number of steps for rough legalization");
  }
  if (roughLegalizationBinSize < 1.0f) {
    throw std::runtime_error(
        "Bin size should generally be larger than 1 (one standard cell)");
  }
  if (roughLegalizationBinSize > 25.0f) {
    throw std::runtime_error(
        "Bin size should not be too large (10 should be enough)");
  }
  if (roughLegalizationLineReoptSize < 1 ||
      roughLegalizationDiagReoptSize < 1 ||
      roughLegalizationSquareReoptSize < 1) {
    throw std::runtime_error(
        "Rough legalization reopt size should be at least 1");
  }
  if (roughLegalizationLineReoptOverlap < 1 ||
      roughLegalizationDiagReoptOverlap < 1 ||
      roughLegalizationSquareReoptOverlap < 1) {
    throw std::runtime_error(
        "Rough legalization reopt overlap should be at least 1");
  }
  if (roughLegalizationLineReoptSize > 64 ||
      roughLegalizationDiagReoptSize > 64 ||
      roughLegalizationSquareReoptSize > 8) {
    throw std::runtime_error("Rough legalization reopt should be small");
  }
  if (roughLegalizationLineReoptSize < 2 &&
      roughLegalizationDiagReoptSize < 2 &&
      roughLegalizationSquareReoptSize < 2 &&
      !roughLegalizationUnidimensionalTransport) {
    throw std::runtime_error(
        "At least one rough legalization reopt value should be 2 or more");
  }
  if (roughLegalizationLineReoptSize > 1 &&
      roughLegalizationLineReoptOverlap >= roughLegalizationLineReoptSize) {
    throw std::runtime_error(
        "Rough legalization reopt overlap should be smaller than reopt size");
  }
  if (roughLegalizationDiagReoptSize > 1 &&
      roughLegalizationDiagReoptOverlap >= roughLegalizationDiagReoptSize) {
    throw std::runtime_error(
        "Rough legalization reopt overlap should be smaller than reopt size");
  }
  if (roughLegalizationSquareReoptSize > 1 &&
      roughLegalizationSquareReoptOverlap >= roughLegalizationSquareReoptSize) {
    throw std::runtime_error(
        "Rough legalization reopt overlap should be smaller than reopt size");
  }
  if (roughLegalizationQuadraticPenalty < 0.0 ||
      roughLegalizationQuadraticPenalty > 1.0) {
    throw std::runtime_error(
        "Rough legalization quadratic penalty should be non-negative and small "
        "(< 1.0)");
  }
  if (roughLegalizationTargetBlending < -0.1 ||
      roughLegalizationTargetBlending > 0.9f) {
    throw std::runtime_error(
        "Rough legalization target blending should generally be between 0 and "
        "0.5");
  }
  if (exportBlending < -0.5f || exportBlending > 1.5f) {
    throw std::runtime_error(
        "Export blending should generally be between 0 and 1");
  }
  if (noise < 0.0 || noise > 2.0) {
    throw std::runtime_error(
        "Noise should be a very small non-negative number");
  }
  if (roughLegalizationCostModel != LegalizationModel::L1 &&
      roughLegalizationUnidimensionalTransport) {
    throw std::runtime_error(
        "Unidimensional transport can only be used with L1 cost model");
  }
}

void GlobalPlacer::place(Circuit &circuit, const GlobalPlacerParameters &params,
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

GlobalPlacer::GlobalPlacer(Circuit &circuit,
                           const GlobalPlacerParameters &params)
    : leg_(DensityLegalizer::fromIspdCircuit(
          circuit, params.roughLegalizationBinSize,
          params.roughLegalizationSideMargin)),
      xtopo_(NetModel::xTopology(circuit)),
      ytopo_(NetModel::yTopology(circuit)),
      params_(params),
      circuit_(circuit) {
  rgen_.seed(params_.seed);
  averageCellLength_ = computeAverageCellSize();
  perCellPenalty_ = computePerCellPenalty();
  DensityLegalizer::Parameters legParams;
  legParams.nbSteps = params.roughLegalizationNbSteps;
  legParams.costModel = params.roughLegalizationCostModel;
  legParams.lineReoptSize = params.roughLegalizationLineReoptSize;
  legParams.lineReoptOverlap = params.roughLegalizationLineReoptOverlap;
  legParams.diagReoptSize = params.roughLegalizationDiagReoptSize;
  legParams.diagReoptOverlap = params.roughLegalizationDiagReoptOverlap;
  legParams.squareReoptSize = params.roughLegalizationSquareReoptSize;
  legParams.squareReoptOverlap = params.roughLegalizationSquareReoptOverlap;
  legParams.unidimensionalTransport =
      params.roughLegalizationUnidimensionalTransport;
  legParams.coarseningLimit = params.roughLegalizationCoarseningLimit;
  LegalizationModel m = params.roughLegalizationCostModel;
  if (m == LegalizationModel::L1 || m == LegalizationModel::L2 ||
      m == LegalizationModel::LInf) {
    float dist = leg_.placementArea().width() + leg_.placementArea().height();
    legParams.quadraticPenaltyFactor =
        params.roughLegalizationQuadraticPenalty / dist;
  }
  leg_.setParams(legParams);
}

void GlobalPlacer::exportPlacement(Circuit &circuit) const {
  assert(xtopo_.nbCells() == circuit.nbCells());
  assert(ytopo_.nbCells() == circuit.nbCells());
  assert(leg_.nbCells() == circuit.nbCells());
  float w = params_.exportBlending;
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
    ret.push_back(
        std::pow(leg_.cellDemand(i) / meanArea, params_.penaltyAreaExponent));
  }
  return ret;
}

std::vector<float> GlobalPlacer::computeIterationPerCellPenalty() {
  std::vector<float> penalty = perCellPenalty_;
  for (float &s : penalty) {
    s *= penalty_;
  }
  if (params_.noise > 0.0) {
    std::uniform_real_distribution<float> dist(0.0f, params_.noise);
    for (float &s : penalty) {
      s *= (1.0f + dist(rgen_));
    }
  }
  return penalty;
}

void GlobalPlacer::run() {
  runInitialLB();
  penalty_ = params_.initialPenalty;
  approximationDistance_ = initialApproximationDistance();
  penaltyCutoffDistance_ = initialPenaltyCutoffDistance();

  float lb = valueLB();
  float ub = std::numeric_limits<float>::infinity();
  int firstStep = params_.nbInitialSteps + 1;
  for (step_ = firstStep; step_ <= params_.maxNbSteps; ++step_) {
    std::cout << "#" << step_ << ":" << std::flush;
    runUB();
    ub = valueUB();
    std::cout << std::defaultfloat << std::setprecision(4) << "\tUB " << ub;

    float dist = leg_.meanDistance();
    std::cout << std::fixed << std::setprecision(1) << "\tDist " << dist;
    std::cout << std::flush;

    float gap = (ub - lb) / ub;
    // Stop if distance or the difference between LB and UB is small enough
    if (gap < params_.gapTolerance || dist < distanceTolerance()) {
      std::cout << std::endl;
      break;
    }
    for (int i = 0; i < params_.nbStepsPerLegalization; ++i) {
      if (i != 0) {
        std::cout << "#" << step_ << ":\t........\t........";
        std::cout << std::flush;
      }
      runLB();
      lb = valueLB();
      std::cout << std::defaultfloat << std::setprecision(4) << "\tLB " << lb
                << std::endl;
    }
    penalty_ *= params_.penaltyUpdateFactor;
    penaltyCutoffDistance_ *= params_.penaltyCutoffDistanceUpdateFactor;
    approximationDistance_ *= params_.approximationDistanceUpdateFactor;
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
  params.netModel = params_.netModel;
  params.tolerance = params_.conjugateGradientErrorTolerance;
  params.maxNbIterations = params_.maxNbConjugateGradientSteps;
  xPlacementLB_ = xtopo_.solveStar(params);
  yPlacementLB_ = ytopo_.solveStar(params);
  std::cout << std::defaultfloat << std::setprecision(4) << "#0:\tLB "
            << valueLB() << std::endl;
  callback(PlacementStep::LowerBound, xPlacementLB_, yPlacementLB_);
  for (step_ = 1; step_ <= params_.nbInitialSteps; ++step_) {
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
  params.netModel = params_.netModel;
  params.approximationDistance = approximationDistance_;
  params.penaltyCutoffDistance = penaltyCutoffDistance_;
  params.tolerance = params_.conjugateGradientErrorTolerance;
  params.maxNbIterations = params_.maxNbConjugateGradientSteps;

  // Compute the per-cell penalty with randomization
  std::vector<float> penalty = computeIterationPerCellPenalty();

  float w = params_.penaltyTargetBlending;
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
  float w = params_.roughLegalizationTargetBlending;
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