#include "place_global.hpp"

#include <cmath>
#include <future>
#include <iostream>
#include <numeric>

#include "density_legalizer.hpp"
#include "net_model.hpp"

namespace coloquinte {

GlobalPlacerParameters::GlobalPlacerParameters(int effort) {
  maxNbSteps = 30;
  gapTolerance = 0.05;
  penaltyCutoffDistance = 10.0;
  initialPenalty = 0.02;
  penaltyUpdateFactor = 1.2;
  netModel = NetModelOption::BoundToBound;
  approximationDistance = 0.5;
  maxNbConjugateGradientSteps = 1000;
  conjugateGradientErrorTolerance = 1.0e-6;
  roughLegalizationCostModel = LegalizationModel::L1;
  nbRoughLegalizationSteps = 1;
  check();
}

void GlobalPlacerParameters::check() const {
  if (maxNbSteps < 0) {
    throw std::runtime_error("Invalid number of steps");
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
  if (nbRoughLegalizationSteps < 0) {
    throw std::runtime_error(
        "Must have non-negative number of steps for rough legalization");
  }
}

void GlobalPlacer::place(Circuit &circuit,
                         const GlobalPlacerParameters &params) {
  GlobalPlacer pl(circuit, params);
  pl.run();
  pl.leg_.exportPlacement(circuit);
}

GlobalPlacer::GlobalPlacer(Circuit &circuit,
                           const GlobalPlacerParameters &params)
    : circuit_(circuit),
      leg_(DensityLegalizer::fromIspdCircuit(circuit)),
      xtopo_(NetModel::xTopology(circuit)),
      ytopo_(NetModel::yTopology(circuit)),
      params_(params) {
  averageCellLength_ = computeAverageCellSize();
  perCellPenalty_ = computePerCellPenalty();
  leg_.setCostModel(params.roughLegalizationCostModel);
  leg_.setNbSteps(params.nbRoughLegalizationSteps);
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
  for (int i = 0; i < leg_.nbCells(); ++i) {
    ret.push_back(std::sqrt(leg_.cellDemand(i) / meanArea));
  }
  return ret;
}

void GlobalPlacer::run() {
  runInitialLB();
  penalty_ = params_.initialPenalty;
  for (step_ = 1; step_ <= params_.maxNbSteps; ++step_) {
    runUB();
    float ub = valueUB();
    runLB();
    float lb = valueLB();
    std::cout << "#" << step_ << ":\tLB " << lb << "\tUB " << ub << std::endl;
    float gap = (ub - lb) / ub;
    if (gap < params_.gapTolerance) break;
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
  params.tolerance = params_.conjugateGradientErrorTolerance;
  params.maxNbIterations = params_.maxNbConjugateGradientSteps;
  xPlacementLB_ = xtopo_.solveStar(params);
  yPlacementLB_ = ytopo_.solveStar(params);
}

void GlobalPlacer::runLB() {
  std::vector<float> penalty = perCellPenalty_;
  for (float &s : penalty) s *= penalty_;
  NetModel::Parameters params;
  params.approximationDistance = approximationDistance();
  params.penaltyCutoffDistance = penaltyCutoffDistance();
  params.tolerance = params_.conjugateGradientErrorTolerance;
  params.maxNbIterations = params_.maxNbConjugateGradientSteps;
  std::future<std::vector<float> > x = std::async(std::launch::async, &NetModel::solve, &xtopo_, xPlacementLB_, xPlacementUB_, penalty, params);
  std::future<std::vector<float> > y = std::async(std::launch::async, &NetModel::solve, &ytopo_, yPlacementLB_, yPlacementUB_, penalty, params);
  xPlacementLB_ = x.get();
  yPlacementLB_ = y.get();
}

void GlobalPlacer::runUB() {
  leg_.updateCellTargetX(xPlacementLB_);
  leg_.updateCellTargetY(yPlacementLB_);
  leg_.run();
  xPlacementUB_ = leg_.simpleCoordX();
  yPlacementUB_ = leg_.simpleCoordY();
}
}  // namespace coloquinte