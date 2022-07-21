#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE CPP_API

#include <boost/test/unit_test.hpp>

#include "place_global/density_grid.hpp"

BOOST_AUTO_TEST_CASE(BasicGrid) {
  Rectangle area(0, 49, 0, 70);
  DensityGrid grid(7, area);
  grid.check();
  BOOST_CHECK_EQUAL(grid.nbBinsX(), 7);
  BOOST_CHECK_EQUAL(grid.nbBinsY(), 10);
}

BOOST_AUTO_TEST_CASE(BasicPlacement) {
  Rectangle area(0, 49, 0, 70);
  DensityGrid grid(7, area);
  grid.check();
  std::vector<int> demands = {1, 2, 3};
  DensityPlacement pl(grid, demands);
  pl.check();
  BOOST_CHECK_EQUAL(pl.nbCells(), 3);
  BOOST_CHECK_EQUAL(pl.nbBinsX(), 7);
  BOOST_CHECK_EQUAL(pl.nbBinsY(), 10);
  pl.binCells(0, 0) = std::vector<int>({0, 1, 2});
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

BOOST_AUTO_TEST_CASE(BasicHierarchicalPlacement2) {
  Rectangle area(0, 49, 0, 70);
  DensityGrid grid(7, area);
  std::vector<int> demands = {1, 2, 3};
  DensityPlacement pl(grid, demands);
  pl.binCells(0, 0) = std::vector<int>({0, 1});
  pl.binCells(0, 1) = std::vector<int>({2});
  HierarchicalDensityPlacement hpl(pl);
  hpl.check();
  BOOST_CHECK_EQUAL(hpl.nbCells(), 3);
  BOOST_CHECK_EQUAL(hpl.nbBinsX(), 7);
  BOOST_CHECK_EQUAL(hpl.nbBinsY(), 10);
  BOOST_CHECK_EQUAL(hpl.levelX(), 0);
  BOOST_CHECK_EQUAL(hpl.levelY(), 0);
  for (int i = 0; i + 1 < hpl.levelX(); ++i) {
    hpl.coarsenX();
    hpl.check();
  }
  for (int i = 0; i + 1 < hpl.levelY(); ++i) {
    hpl.coarsenY();
    hpl.check();
  }
  for (int i = 0; i + 1 < hpl.levelX(); ++i) {
    hpl.refineX();
    hpl.check();
  }
  for (int i = 0; i + 1 < hpl.levelY(); ++i) {
    hpl.refineY();
    hpl.check();
  }
}
