
#include <fstream>

#include "coloquinte.hpp"

namespace coloquinte {

void exportIspdAux(const std::string &filename) {
  std::ofstream f(filename + ".aux");
  f << "RowBasedPlacement : " << filename << ".nodes " << filename << ".nets "
    << filename << ".pl " << filename << ".scl" << std::endl;
}

void exportIspdNodes(const Circuit &circuit, const std::string &filename) {
  std::ofstream f(filename + ".nodes");
  f << "UCLA nodes 1.0\n\n";
  f << "NumNodes : " << circuit.nbCells() << "\n";
  int numTerminals = 0;
  for (int i = 0; i < circuit.nbCells(); ++i) {
    if (circuit.isFixed(i)) {
      ++numTerminals;
    }
  }
  f << "NumTerminals : " << numTerminals << "\n";

  for (int i = 0; i < circuit.nbCells(); ++i) {
    f << "\to" << i << "\t" << circuit.cellWidth_[i] << "\t"
      << circuit.cellHeight_[i];
    if (circuit.isFixed(i)) {
      f << "\tterminal";
    }
    f << "\n";
  }
}

void exportIspdPlace(const Circuit &circuit, const std::string &filename) {
  std::ofstream f(filename + ".pl");
  f << "UCLA pl 1.0\n\n";

  for (int i = 0; i < circuit.nbCells(); ++i) {
    f << "o" << i << "\t" << circuit.x(i) << "\t" << circuit.y(i)
      << "\t: " << toString(circuit.orientation(i)) << "\n";
  }
}

void exportIspdNets(const Circuit &circuit, const std::string &filename) {
  std::ofstream f(filename + ".nets");
  f << "UCLA nets 1.0\n\n";

  f << "NumNets : " << circuit.nbNets() << "\n";
  f << "NumPins : " << circuit.nbPins() << "\n\n";

  for (int i = 0; i < circuit.nbNets(); ++i) {
    f << "NetDegree : " << circuit.nbPinsNet(i) << " n" << i << "\n";
    for (int j = 0; j < circuit.nbPinsNet(i); ++j) {
      int c = circuit.pinCell(i, j);
      f << "\to" << c << " I : ";
      double x = circuit.pinXOffset(i, j) - 0.5 * circuit.cellWidth_[c];
      double y = circuit.pinYOffset(i, j) - 0.5 * circuit.cellHeight_[c];
      f << x << " " << y << "\n";
    }
  }
}

void exportIspdRows(const Circuit &circuit, const std::string &filename) {
  std::ofstream f(filename + ".scl");
  f << "UCLA scl 1.0\n\n";

  f << "NumRows : " << circuit.nbRows() << "\n\n";
  for (int i = 0; i < circuit.nbRows(); ++i) {
    f << "CoreRow Horizontal\n";
    f << "  Coordinate    : " << circuit.rows()[i].minY << "\n";
    f << "  Height        : " << circuit.rows()[i].height() << "\n";
    f << "  Sitewidth     : 1\n";
    f << "  Sitespacing   : 1\n";
    f << "  Siteorient    : 1\n";
    f << "  Sitesymmetry  : 1\n";
    f << "  SubrowOrigin  : " << circuit.rows()[i].minX
      << "     NumSites : " << circuit.rows()[i].width() << "\n";
    f << "End\n";
  }
}

void Circuit::exportIspd(const std::string &filename) const {
  exportIspdAux(filename);
  exportIspdNodes(*this, filename);
  exportIspdPlace(*this, filename);
  exportIspdNets(*this, filename);
  exportIspdRows(*this, filename);
}
}  // namespace coloquinte