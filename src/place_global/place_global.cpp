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
  if (initialPenalty <= 0.0f) {
    throw std::runtime_error("Initial penalty should be positive");
  }
  if (penaltyUpdateFactor <= 1.0f || penaltyUpdateFactor >= 2.0f) {
    throw std::runtime_error(
        "Penalty update factor should be between one and two");
  }
  if (approximationDistance < 1.0e-6) {
    throw std::runtime_error(
        "Too small approximation distance may lead to issues");
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
  if (roughLegalizationReoptLength < 1 ||
      roughLegalizationReoptSquareSize < 1) {
    throw std::runtime_error("Rough legalization reopt should be at least 1");
  }
  if (roughLegalizationReoptLength > 64 ||
      roughLegalizationReoptSquareSize > 8) {
    throw std::runtime_error("Rough legalization reopt should be small");
  }
  if (roughLegalizationReoptLength < 2 &&
      roughLegalizationReoptSquareSize < 2) {
    throw std::runtime_error(
        "At least one rough legalization reopt value should be 2 or more");
  }
  if (roughLegalizationQuadraticPenalty < 0.0 ||
      roughLegalizationQuadraticPenalty > 1.0) {
    throw std::runtime_error(
        "Rough legalization quadratic penalty should be non-negative and small "
        "(< 1.0)");
  }
  if (exportWeighting < -0.5f || exportWeighting > 1.5f) {
    throw std::runtime_error(
        "Export weighting should generally be between 0 and 1");
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
  legParams.reoptimizationLength = params.roughLegalizationReoptLength;
  legParams.reoptimizationSquareSize = params.roughLegalizationReoptSquareSize;
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
  float w = params_.exportWeighting;
  std::vector<float> xplace;
  std::vector<float> yplace;
  for (int i = 0; i < circuit.nbCells(); ++i) {
    xplace.push_back(((1.0f - w) * xPlacementLB_[i] + w * xPlacementUB_[i]));
    yplace.push_back(((1.0f - w) * yPlacementLB_[i] + w * yPlacementUB_[i]));
  }
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
    ret.push_back(std::sqrt(leg_.cellDemand(i) / meanArea));
  }
  return ret;
}

void GlobalPlacer::run() {
  runInitialLB();
  penalty_ = params_.initialPenalty;
  float lb = valueLB();
  for (step_ = params_.nbInitialSteps + 1; step_ <= params_.maxNbSteps;
       ++step_) {
    std::cout << "#" << step_ << std::flush;
    runUB();
    float ub = valueUB();
    float dist = leg_.meanDistance();
    std::cout << std::defaultfloat << std::setprecision(4) << ":\tUB " << ub;
    std::cout << std::fixed << std::setprecision(1) << "\tDist " << dist;
    std::cout << std::flush;

    float gap = (ub - lb) / ub;
    // Stop if distance or the difference between LB and UB is small enough
    if (gap < params_.gapTolerance || dist < distanceTolerance()) {
      std::cout << std::endl;
      break;
    }
    runLB();
    lb = valueLB();
    std::cout << std::defaultfloat << std::setprecision(4) << "\tLB " << lb
              << std::endl;
    penalty_ *= params_.penaltyUpdateFactor;
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
}

void GlobalPlacer::runLB() {
  // Compute the parameters for the continuous model solver
  NetModel::Parameters params;
  params.netModel = params_.netModel;
  params.approximationDistance = approximationDistance();
  params.penaltyCutoffDistance = penaltyCutoffDistance();
  params.tolerance = params_.conjugateGradientErrorTolerance;
  params.maxNbIterations = params_.maxNbConjugateGradientSteps;

  // Compute the per-cell penalty with randomization
  std::vector<float> penalty = perCellPenalty_;
  for (float &s : penalty) {
    s *= penalty_;
  }
  float rand = std::uniform_real_distribution<float>(-1.0e-4f, 1.0e-4f)(rgen_);
  for (float &s : penalty) {
    s *= (1.0f + rand);
  }

  // Solve the continuous model (x and y independently)
  std::future<std::vector<float> > x =
      std::async(std::launch::async, &NetModel::solveWithPenalty, &xtopo_,
                 xPlacementLB_, xPlacementUB_, penalty, params);
  std::future<std::vector<float> > y =
      std::async(std::launch::async, &NetModel::solveWithPenalty, &ytopo_,
                 yPlacementLB_, yPlacementUB_, penalty, params);
  xPlacementLB_ = x.get();
  yPlacementLB_ = y.get();
  callback(PlacementStep::LowerBound, xPlacementLB_, yPlacementLB_);
}

void GlobalPlacer::runUB() {
  leg_.updateCellTargetX(xPlacementLB_);
  leg_.updateCellTargetY(yPlacementLB_);
  leg_.run();
  xPlacementUB_ = leg_.spreadCoordX(xPlacementLB_);
  yPlacementUB_ = leg_.spreadCoordY(yPlacementLB_);
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