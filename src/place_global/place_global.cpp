#include "place_global.hpp"

#include <cmath>
#include <iostream>
#include <numeric>

#include "density_legalizer.hpp"
#include "net_model.hpp"

void GlobalPlacer::place(Circuit &circuit, int effort) {
  GlobalPlacer pl(circuit);
  pl.initParameters();
  pl.runInitialLB();
  for (pl.step_ = 1; pl.step_ <= pl.maxNbSteps_; ++pl.step_) {
    pl.runUB();
    pl.runLB();
    std::cout << "#" << pl.step_ << ":\tLB " << pl.valueLB() << "\tUB "
              << pl.valueUB() << std::endl;
    pl.updateParameters();
  }
  pl.runUB();
  pl.leg_.exportPlacement(circuit);
}

GlobalPlacer::GlobalPlacer(Circuit &circuit)
    : circuit_(circuit),
      leg_(DensityLegalizer::fromIspdCircuit(circuit)),
      xtopo_(NetModel::xTopology(circuit)),
      ytopo_(NetModel::yTopology(circuit)) {
  leg_ = DensityLegalizer::fromIspdCircuit(circuit);
  initParameters();
}

std::vector<float> GlobalPlacer::computeBaseForces() const {
  int nbCells = leg_.nbNonEmptyCells();
  float meanArea = leg_.totalDemand() / std::max(1, nbCells);
  std::vector<float> ret;
  for (int i = 0; i < leg_.nbCells(); ++i) {
    ret.push_back(std::sqrt(leg_.cellDemand(i) / meanArea));
  }
  return ret;
}

void GlobalPlacer::initParameters() {
  baseForces_ = computeBaseForces();
  float totalDemand = leg_.totalDemand();
  float avgDemand = totalDemand == 0.0 ? 0.0 : totalDemand / leg_.nbCells();
  float avgDist = std::sqrt(avgDemand);
  epsilon_ = avgDist / 2.0;
  cutoffDistance_ = avgDist * 10.0;
  updateFactor_ = 1.2;
  maxNbSteps_ = 30;
  forceFactor_ = 0.02;
  step_ = 0;
}

void GlobalPlacer::updateParameters() { forceFactor_ *= updateFactor_; }

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
  std::vector<float> strength = baseForces_;
  for (float &s : strength) s *= forceFactor_;
  xPlacementLB_ = xtopo_.solveB2B(xPlacementLB_, epsilon_, xPlacementUB_,
                                  strength, cutoffDistance_);
  yPlacementLB_ = ytopo_.solveB2B(yPlacementLB_, epsilon_, yPlacementUB_,
                                  strength, cutoffDistance_);
}

void GlobalPlacer::runUB() {
  leg_.updateCellTargetX(xPlacementLB_);
  leg_.updateCellTargetY(yPlacementLB_);
  leg_.run();
  xPlacementUB_ = leg_.simpleCoordX();
  yPlacementUB_ = leg_.simpleCoordY();
}