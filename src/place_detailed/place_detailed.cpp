#include "legalizer.hpp"
#include "place_detailed.hpp"

#include <iostream>

void DetailedPlacer::place(Circuit &circuit) {
    Legalizer leg = Legalizer::fromIspdCircuit(circuit);
    leg.run();
    leg.exportPlacement(circuit);
}
