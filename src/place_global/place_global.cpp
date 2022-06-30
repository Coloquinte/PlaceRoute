#include "density_legalizer.hpp"
#include "wirelength_model.hpp"
#include "place_global.hpp"


std::pair<std::vector<float>, std::vector<float> > GlobalPlacer::place(const Circuit &circuit) {
    float epsilon = 1.0;

    auto xtopo = NetWirelength::xTopology(circuit);
    auto ytopo = NetWirelength::yTopology(circuit);

    auto xplace = tensorToVector(xtopo.starSolve());
    auto yplace = tensorToVector(ytopo.starSolve());

    DensityLegalizer leg(circuit);
    for (int i = 0; i < 100; ++i) {
        xplace = tensorToVector(xtopo.b2bSolve(vectorToTensor(xplace), epsilon));
        yplace = tensorToVector(ytopo.b2bSolve(vectorToTensor(yplace), epsilon));
        // TODO
        leg.updateCellTargetX(xplace);
        leg.updateCellTargetY(yplace);
    }

    return std::make_pair(xplace, yplace);
}
