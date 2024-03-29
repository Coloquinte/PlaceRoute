#define BOOST_TEST_MODULE EXPANSION

#include <boost/test/unit_test.hpp>

#include "coloquinte.hpp"

using namespace coloquinte;

BOOST_AUTO_TEST_CASE(TestDensityExpansion) {
  Circuit circuit(1);
  circuit.setCellWidth({20});
  circuit.setCellHeight({10});
  circuit.setRows({Row(0, 40, 0, 10, CellOrientation::N)});
  circuit.expandCellsToDensity(0.75, 0.0);
  BOOST_CHECK_EQUAL(circuit.cellWidth()[0], 30);
}

BOOST_AUTO_TEST_CASE(TestDensityExpansionWithMargin) {
  Circuit circuit(1);
  circuit.setCellWidth({20});
  circuit.setCellHeight({10});
  circuit.setRows({Row(0, 50, 0, 10, CellOrientation::N)});
  circuit.expandCellsToDensity(0.75, 0.5);
  BOOST_CHECK_EQUAL(circuit.cellWidth()[0], 30);
}

BOOST_AUTO_TEST_CASE(TestFactorExpansion) {
  Circuit circuit(1);
  circuit.setCellWidth({20});
  circuit.setCellHeight({10});
  circuit.setRows({Row(0, 40, 0, 10, CellOrientation::N)});
  circuit.expandCellsByFactor({1.5f});
  BOOST_CHECK_EQUAL(circuit.cellWidth()[0], 30);
}

BOOST_AUTO_TEST_CASE(TestFactorExpansionWithMargin) {
  Circuit circuit(1);
  circuit.setCellWidth({20});
  circuit.setCellHeight({10});
  circuit.setRows({Row(0, 50, 0, 10, CellOrientation::N)});
  circuit.expandCellsByFactor({1.5f}, 1.0, 0.5);
  BOOST_CHECK_EQUAL(circuit.cellWidth()[0], 30);
}

BOOST_AUTO_TEST_CASE(TestFactorExpansionTooLarge) {
  Circuit circuit(1);
  circuit.setCellWidth({20});
  circuit.setCellHeight({10});
  circuit.setRows({Row(0, 40, 0, 10, CellOrientation::N)});
  circuit.expandCellsByFactor({10.0f});
  BOOST_CHECK_EQUAL(circuit.cellWidth()[0], 40);
}

BOOST_AUTO_TEST_CASE(TestFactorExpansionMultiple) {
  Circuit circuit(2);
  circuit.setCellWidth({20, 12});
  circuit.setCellHeight({10, 10});
  circuit.setRows({Row(0, 100, 0, 10, CellOrientation::N)});
  circuit.expandCellsByFactor({1.5f, 1.34f}, 1.0, 0.0);
  BOOST_CHECK_EQUAL(circuit.cellWidth()[0], 30);
  BOOST_CHECK_EQUAL(circuit.cellWidth()[1], 16);
}

BOOST_AUTO_TEST_CASE(TestExpansionComputation) {
  Circuit circuit(4);
  circuit.setCellWidth({20, 20, 20, 20});
  circuit.setCellHeight({10, 10, 10, 10});
  circuit.setCellX({40, 9, 10, -5});
  circuit.setCellY({-100, 19, 20, 19});
  std::vector<std::pair<Rectangle, float> > congestionMap;
  congestionMap.emplace_back(Rectangle({0, 10, 10, 20}), 1.2);
  congestionMap.emplace_back(Rectangle({-5, 0, 10, 20}), 1.5);
  std::vector<float> expansion = circuit.computeCellExpansion(congestionMap);
  BOOST_CHECK_EQUAL(expansion.size(), circuit.nbCells());
  BOOST_CHECK_EQUAL(expansion[0], 1.0f);
  BOOST_CHECK_EQUAL(expansion[1], 1.2f);
  BOOST_CHECK_EQUAL(expansion[2], 1.0f);
  BOOST_CHECK_EQUAL(expansion[3], 1.5f);
}
