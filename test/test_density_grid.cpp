#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE CPP_API

#include <boost/test/unit_test.hpp>

#include "place_global/density_grid.hpp"

using namespace coloquinte;

BOOST_AUTO_TEST_CASE(BasicGrid) {
  Rectangle area(0, 49, 0, 70);
  DensityGrid grid(7, area);
  grid.check();
  BOOST_CHECK_EQUAL(grid.nbBinsX(), 7);
  BOOST_CHECK_EQUAL(grid.nbBinsY(), 10);
}

BOOST_AUTO_TEST_CASE(BasicHierarchicalPlacement1) {
  Rectangle area(0, 49, 0, 70);
  DensityGrid grid(7, area);
  std::vector<int> demands = {1, 2, 3};
  HierarchicalDensityPlacement hpl(grid, demands);
  hpl.check();
  BOOST_CHECK_EQUAL(hpl.nbCells(), 3);
  BOOST_CHECK_EQUAL(hpl.nbBinsX(), 1);
  BOOST_CHECK_EQUAL(hpl.nbBinsY(), 1);
  BOOST_CHECK_EQUAL(hpl.levelX(), hpl.nbLevelX() - 1);
  BOOST_CHECK_EQUAL(hpl.levelY(), hpl.nbLevelY() - 1);
  for (int i = 0; i + 1 < hpl.nbLevelX(); ++i) {
    hpl.refineX();
    hpl.check();
  }
  for (int i = 0; i + 1 < hpl.nbLevelY(); ++i) {
    hpl.refineY();
    hpl.check();
  }
  for (int i = 0; i + 1 < hpl.nbLevelX(); ++i) {
    hpl.coarsenX();
    hpl.check();
  }
  for (int i = 0; i + 1 < hpl.nbLevelY(); ++i) {
    hpl.coarsenY();
    hpl.check();
  }
}

BOOST_AUTO_TEST_CASE(HierarchicalPlacementFindCoord) {
  Rectangle area(0, 49, -20, 70);
  DensityGrid grid(7, area);
  std::vector<int> demands = {1, 2, 3};
  HierarchicalDensityPlacement hpl(grid, demands);
  hpl.check();
  hpl.refineX();
  hpl.refineX();
  hpl.refineY();
  for (int i = -10; i < 60; ++i) {
    hpl.findBinByX(i);
  }
  for (int i = -30; i < 80; ++i) {
    hpl.findBinByY(i);
  }
}