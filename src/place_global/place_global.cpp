#include "place_global.hpp"

#include <cassert>
#include <chrono>
#include <cmath>
#include <future>
#include <iostream>
#include <numeric>

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
  if (gapTolerance < 0.0f || gapTolerance > 2.0f) {
    throw std::runtime_error("Invalid gap tolerance");
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
  if (exportWeighting < -0.5f || exportWeighting > 1.5f) {
    throw std::runtime_error(
        "Export weighting should generally be between 0 and 1");
  }
}

void GlobalPlacer::place(Circuit &circuit,
                         const GlobalPlacerParameters &params) {
  params.check();
  GlobalPlacer pl(circuit, params);
  pl.run();
  pl.exportPlacement(circuit);
}

GlobalPlacer::GlobalPlacer(Circuit &circuit,
                           const GlobalPlacerParameters &params)
    : circuit_(circuit),
      leg_(DensityLegalizer::fromIspdCircuit(circuit,
                                             params.roughLegalizationBinSize)),
      xtopo_(NetModel::xTopology(circuit)),
      ytopo_(NetModel::yTopology(circuit)),
      params_(params) {
  rgen_.seed(params_.seed);
  averageCellLength_ = computeAverageCellSize();
  perCellPenalty_ = computePerCellPenalty();
  leg_.setCostModel(params.roughLegalizationCostModel);
  leg_.setNbSteps(params.roughLegalizationNbSteps);
}

void GlobalPlacer::exportPlacement(Circuit &circuit) {
  assert(leg_.nbCells() == circuit.nbCells());
  float w = params_.exportWeighting;
  for (int i = 0; i < circuit.nbCells(); ++i) {
    if (circuit.isFixed(i)) {
      continue;
    }
    circuit.cellX_[i] =
        std::round(((1.0f - w) * xPlacementLB_[i] + w * xPlacementUB_[i]));
    circuit.cellY_[i] =
        std::round(((1.0f - w) * yPlacementLB_[i] + w * yPlacementUB_[i]));
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
  std::cout << "Global placement starting" << std::endl;
  auto startTime = std::chrono::steady_clock::now();
  runInitialLB();
  penalty_ = params_.initialPenalty;
  for (step_ = params_.nbInitialSteps + 1; step_ <= params_.maxNbSteps;
       ++step_) {
    runUB();
    float ub = valueUB();
    runLB();
    float lb = valueLB();
    std::cout << "#" << step_ << ":\tLB " << lb << "\tUB " << ub << std::endl;
    float gap = (ub - lb) / ub;
    if (gap < params_.gapTolerance) {
      break;
    }
    penalty_ *= params_.penaltyUpdateFactor;
  }
  runUB();
  auto endTime = std::chrono::steady_clock::now();
  std::chrono::duration<float> duration = endTime - startTime;
  std::cout << "Global placement done in " << duration.count() << "s"
            << std::endl;
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
  std::cout << "#0:\tLB " << valueLB() << std::endl;
  for (step_ = 1; step_ <= params_.nbInitialSteps; ++step_) {
    xPlacementLB_ = xtopo_.solve(xPlacementLB_, params);
    yPlacementLB_ = ytopo_.solve(yPlacementLB_, params);
    std::cout << "#" << step_ << ":\tLB " << valueLB() << std::endl;
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
}

void GlobalPlacer::runUB() {
  leg_.updateCellTargetX(xPlacementLB_);
  leg_.updateCellTargetY(yPlacementLB_);
  leg_.run();
  xPlacementUB_ = leg_.spreadCoordX(xPlacementLB_);
  yPlacementUB_ = leg_.spreadCoordY(yPlacementLB_);
}
}  // namespace coloquinte