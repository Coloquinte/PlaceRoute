#define BOOST_TEST_MODULE DENSITY_LEGALIZER

#include <boost/test/unit_test.hpp>
#include <cmath>
#include <iostream>
#include <limits>
#include <random>
#include <stdexcept>
#include <vector>

#include "place_global/density_legalizer.hpp"

using namespace coloquinte;

BOOST_AUTO_TEST_CASE(EqualSplit) {
  Rectangle area(0, 10, 0, 5);
  std::vector<int> cellDemand = {4, 4, 4, 4};
  std::vector<float> cellTargetX = {0.0, 0.0, 6.0, 6.0};
  std::vector<float> cellTargetY = {0.0, 0.0, 0.0, 0.0};
  DensityGrid grid(5, area);
  DensityLegalizer leg(grid, cellDemand);
  leg.updateCellDemand(cellDemand);
  leg.updateCellTargetX(cellTargetX);
  leg.updateCellTargetY(cellTargetY);
  leg.run();
  BOOST_CHECK_EQUAL(leg.binCells(0, 0).size(), 2);
  BOOST_CHECK_EQUAL(leg.binCells(1, 0).size(), 2);
  BOOST_CHECK_EQUAL(leg.binUsage(0, 0), 8);
  BOOST_CHECK_EQUAL(leg.binUsage(1, 0), 8);
  BOOST_CHECK_EQUAL(leg.cellBinX(0), 0);
  BOOST_CHECK_EQUAL(leg.cellBinX(1), 0);
  BOOST_CHECK_EQUAL(leg.cellBinX(2), 1);
  BOOST_CHECK_EQUAL(leg.cellBinX(3), 1);
}

BOOST_AUTO_TEST_CASE(OneSide) {
  Rectangle area(0, 8, 0, 2);
  std::vector<int> cellDemand = {4, 4, 4, 4};
  std::vector<float> cellTargetX = {0.0, 6.0, 6.0, 6.0};
  std::vector<float> cellTargetY = {0.0, 0.0, 0.0, 0.0};
  DensityGrid grid(4, area);
  DensityLegalizer leg(grid, cellDemand);
  leg.updateCellDemand(cellDemand);
  leg.updateCellTargetX(cellTargetX);
  leg.updateCellTargetY(cellTargetY);
  leg.run();
  BOOST_CHECK_EQUAL(leg.binCells(0, 0).size(), 2);
  BOOST_CHECK_EQUAL(leg.binCells(1, 0).size(), 2);
  BOOST_CHECK_EQUAL(leg.binCapacity(0, 0), 8);
  BOOST_CHECK_EQUAL(leg.binCapacity(1, 0), 8);
  BOOST_CHECK_EQUAL(leg.binUsage(0, 0), 8);
  BOOST_CHECK_EQUAL(leg.binUsage(1, 0), 8);
  BOOST_CHECK_EQUAL(leg.cellBinX(0), 0);
  BOOST_CHECK_EQUAL(leg.cellBinX(1), 0);
  BOOST_CHECK_EQUAL(leg.cellBinX(2), 1);
  BOOST_CHECK_EQUAL(leg.cellBinX(3), 1);
}

BOOST_AUTO_TEST_CASE(Random) {
  Rectangle area(0, 10, 0, 10);
  std::vector<int> cellDemand = {4, 4, 4, 4};
  std::vector<float> cellTargetX = {0.0, 6.0, 0.0, 6.0};
  std::vector<float> cellTargetY = {0.0, 0.0, 0.0, 0.0};
  DensityGrid grid(5, area);
  DensityLegalizer leg(grid, cellDemand);
  leg.updateCellDemand(cellDemand);
  leg.updateCellTargetX(cellTargetX);
  leg.updateCellTargetY(cellTargetY);
  leg.run();
  BOOST_CHECK_EQUAL(leg.cellBinX(0), 0);
  BOOST_CHECK_EQUAL(leg.cellBinX(1), 1);
  BOOST_CHECK_EQUAL(leg.cellBinX(2), 0);
  BOOST_CHECK_EQUAL(leg.cellBinX(3), 1);
}

BOOST_AUTO_TEST_CASE(Coord) {
  Rectangle area(0, 10, 0, 5);
  DensityGrid grid(1, area);
  HierarchicalDensityPlacement pl(grid, 0);
  pl.refineX();
  pl.refineY();
  pl.refineX();
  BOOST_CHECK_EQUAL(pl.nbBinsX(), 4);
  BOOST_CHECK_EQUAL(pl.nbBinsY(), 2);
  BOOST_CHECK_CLOSE(pl.binX(0, 0), 1.0f, 0.0001f);
  BOOST_CHECK_CLOSE(pl.binX(1, 0), 3.5f, 0.0001f);
  BOOST_CHECK_CLOSE(pl.binX(2, 0), 6.0f, 0.0001f);
  BOOST_CHECK_CLOSE(pl.binX(3, 0), 8.5f, 0.0001f);
  BOOST_CHECK_CLOSE(pl.binY(0, 0), 1.0f, 0.0001f);
  BOOST_CHECK_CLOSE(pl.binY(0, 1), 3.5f, 0.0001f);
}

BOOST_AUTO_TEST_CASE(UnitArea) {
  Rectangle area(0, 10, 0, 5);
  DensityGrid grid(1, area);
  DensityLegalizer pl(grid, std::vector<int>(50, 1));
  pl.run();
  pl.check();
  BOOST_CHECK_EQUAL(pl.totalOverflow(), 0);
  pl.coarsenFully();
  pl.check();
}
