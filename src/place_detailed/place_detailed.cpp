#include "place_detailed.hpp"

#include <iostream>

#include "legalizer.hpp"

void DetailedPlacer::place(Circuit &circuit) {
  Legalizer leg = Legalizer::fromIspdCircuit(circuit);
  leg.run();
  leg.exportPlacement(circuit);
}
