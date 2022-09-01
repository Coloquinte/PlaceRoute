#include "place_global.hpp"

#include <cmath>
#include <iostream>
#include <numeric>

#include "density_legalizer.hpp"
#include "net_model.hpp"

namespace coloquinte {

GlobalPlacer::Parameters::Parameters(int effort) {
  maxNbSteps = 30;
  gapTolerance = 0.05;
  penaltyCutoffDistance = 10.0;
  initialPenalty = 0.02;
  penaltyUpdateFactor = 1.2;
  approximationDistance = 0.5;
  maxNbConjugateGradientSteps = 100;
  conjugateGradientErrorTolerance = 1.0e-3;
  check();
}

void GlobalPlacer::Parameters::check() const {
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
    throw std::runtime_error("Penalty update factor should be between one and two");
  }
  if (approximationDistance < 1.0e-6) {
    throw std::runtime_error("Too small approximation distance may lead to issues");
  }
  if (approximationDistance > 1.0e3) {
    throw std::runtime_error("Too large approximation distance is highly imprecise");
  }
  if (maxNbConjugateGradientSteps <= 0) {
    throw std::runtime_error("Must have positive number of steps during conjugate gradients");
  }
  if (conjugateGradientErrorTolerance < 1.0e-6) {
    throw std::runtime_error("Too small error tolerance may lead to issues");
  }
  if (conjugateGradientErrorTolerance > 1.0) {
    throw std::runtime_error("Too large error tolerance is highly imprecise");
  }
}

void GlobalPlacer::place(Circuit &circuit,
                         const GlobalPlacer::Parameters &params) {
  GlobalPlacer pl(circuit, params);
  pl.run();
  pl.leg_.exportPlacement(circuit);
}

GlobalPlacer::GlobalPlacer(Circuit &circuit,
                           const GlobalPlacer::Parameters &params)
    : circuit_(circuit),
      leg_(DensityLegalizer::fromIspdCircuit(circuit)),
      xtopo_(NetModel::xTopology(circuit)),
      ytopo_(NetModel::yTopology(circuit)),
      params_(params) {
  averageCellLength_ = computeAverageCellSize();
  perCellPenalty_ = computePerCellPenalty();
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
    std::cout << "#" << step_ << ":\tLB " << lb << "\tUB " << ub
              << std::endl;
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
  xPlacementLB_ = xtopo_.solveStar();
  yPlacementLB_ = ytopo_.solveStar();
}

void GlobalPlacer::runLB() {
  std::vector<float> penalty = perCellPenalty_;
  for (float &s : penalty) s *= penalty_;
  // TODO: pass other parameters to the solver
  xPlacementLB_ =
      xtopo_.solveB2B(xPlacementLB_, approximationDistance(), xPlacementUB_,
                      penalty, penaltyCutoffDistance());
  yPlacementLB_ =
      ytopo_.solveB2B(yPlacementLB_, approximationDistance(), yPlacementUB_,
                      penalty, penaltyCutoffDistance());
}

void GlobalPlacer::runUB() {
  leg_.updateCellTargetX(xPlacementLB_);
  leg_.updateCellTargetY(yPlacementLB_);
  leg_.run();
  xPlacementUB_ = leg_.simpleCoordX();
  yPlacementUB_ = leg_.simpleCoordY();
}
}  // namespace coloquinte