#include "place_global.hpp"

#include <cmath>
#include <iostream>
#include <numeric>

#include "density_legalizer.hpp"
#include "net_model.hpp"

std::vector<float> getBaseForces(const Circuit &circuit) {
  std::vector<float> ret;
  for (int i = 0; i < circuit.nbCells(); ++i) {
    ret.push_back(circuit.area(i));
  }
  float totArea = std::accumulate(ret.begin(), ret.end(), 0.0f);
  for (float &a : ret) {
    a /= totArea;
  }
  return ret;
}

void GlobalPlacer::place(Circuit &circuit) {
  float epsilon = 1.0;
  float cutoffDistance = 1000.0;
  int nbSteps = 10;
  auto baseForces = getBaseForces(circuit);

  auto xtopo = NetModel::xTopology(circuit);
  auto ytopo = NetModel::yTopology(circuit);

  auto xplace = xtopo.solveStar();
  auto yplace = ytopo.solveStar();

  DensityLegalizer leg = DensityLegalizer::fromIspdCircuit(circuit);
  for (int i = 0; i < nbSteps; ++i) {
    float wirelengthPlace = xtopo.value(xplace) + ytopo.value(yplace);
    std::cout << "LB wirelength #" << i << ": " << wirelengthPlace << std::endl;
    leg.updateCellTargetX(xplace);
    leg.updateCellTargetY(yplace);
    leg.run();
    float forceFactor = 0.001 * i;
    auto xtarget = leg.simpleCoordX();
    auto ytarget = leg.simpleCoordY();
    float wirelengthTarget = xtopo.value(xtarget) + ytopo.value(ytarget);
    std::cout << "UB wirelength #" << i << ": " << wirelengthTarget
              << std::endl;
    std::vector<float> strength = baseForces;
    for (float &s : strength) s *= forceFactor;
    xplace = xtopo.solveB2B(xplace, epsilon, xtarget, strength, cutoffDistance);
    yplace = ytopo.solveB2B(yplace, epsilon, ytarget, strength, cutoffDistance);
  }

  std::cout << "LB wirelength final: "
            << xtopo.value(xplace) + ytopo.value(yplace) << std::endl;

  xtopo.exportPlacementX(circuit, xplace);
  ytopo.exportPlacementY(circuit, yplace);

  std::cout << "Circuit wirelength : " << circuit.hpwl() << std::endl;
}
