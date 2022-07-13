#include "density_legalizer.hpp"
#include "wirelength_model.hpp"
#include "place_global.hpp"

#include <iostream>

std::vector<float> getBaseForces(const Circuit &circuit) {
    std::vector<float> ret;
    for (int i = 0; i < circuit.nbCells(); ++i) {
        ret.push_back(circuit.getArea(i));
    }
    float totArea = std::accumulate(ret.begin(), ret.end(), 0.0f);
    for (float &a : ret) {
        a /= totArea;
    }
    return ret;
}


std::pair<std::vector<float>, std::vector<float> > GlobalPlacer::place(const Circuit &circuit) {
    float epsilon = 1.0;
    float cutoffDistance = 1000.0;
    auto baseForces = vectorToTensor(getBaseForces(circuit));

    auto xtopo = NetWirelength::xTopology(circuit);
    auto ytopo = NetWirelength::yTopology(circuit);

    auto xplace = xtopo.starSolve();
    auto yplace = ytopo.starSolve();

    DensityLegalizer leg = DensityLegalizer::fromIspdCircuit(circuit);
    for (int i = 0; i < 100; ++i) {
        leg.updateCellTargetX(tensorToVector(xplace));
        leg.updateCellTargetY(tensorToVector(yplace));
        leg.assign();
        float forceFactor = 0.001 * i;
        auto xtarget = vectorToTensor(leg.simpleCoordX());
        auto ytarget = vectorToTensor(leg.simpleCoordY());
        float wirelengthPlace = xtopo.valueHPWL(xplace) + ytopo.valueHPWL(yplace);
        float wirelengthTarget = xtopo.valueHPWL(xtarget) + ytopo.valueHPWL(ytarget);
        std::cout << "LB wirelength #" << i << ": " << wirelengthPlace << std::endl;
        std::cout << "UB wirelength #" << i << ": " << wirelengthTarget << std::endl;
        xplace = xtopo.b2bSolvePenalized(xplace, epsilon, xtarget, baseForces * forceFactor, cutoffDistance);
        yplace = ytopo.b2bSolvePenalized(yplace, epsilon, ytarget, baseForces * forceFactor, cutoffDistance);
    }

    return std::make_pair(tensorToVector(xplace), tensorToVector(yplace));
}
