#include "place_detailed.hpp"

#include <iostream>

#include "legalizer.hpp"

namespace coloquinte {
void DetailedPlacer::place(Circuit &circuit, int effort) {
  Legalizer leg = Legalizer::fromIspdCircuit(circuit);
  leg.run();
  leg.exportPlacement(circuit);
}
}
