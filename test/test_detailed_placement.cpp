#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE CPP_API

#include <boost/test/unit_test.hpp>

#include "place_detailed/detailed_placement.hpp"

using namespace coloquinte;

BOOST_AUTO_TEST_CASE(TestConstruction) {
  // 3 cells in two rows:
  // Row 0, x from 2 to 12, y 10, with cell 1 from 2 to 8 and cell 0 from 8 to
  // 12 Row 1, x from 2 to 12, y 20, with cell 2 from 2 to 7
  std::vector<int> widths = {4, 6, 5};
  std::vector<int> cellX = {8, 2, 2};
  std::vector<int> cellY = {10, 10, 20};
  std::vector<Rectangle> rows = {{2, 12, 10, 20}, {2, 12, 20, 30}};
  DetailedPlacement pl(rows, widths, cellX, cellY);
  pl.check();
}

BOOST_AUTO_TEST_CASE(TestBeforeRow) {
  std::vector<int> widths = {2};
  std::vector<int> cellX = {5};
  std::vector<int> cellY = {10};
  std::vector<Rectangle> rows = {{6, 12, 10, 20}};
  BOOST_CHECK_THROW(DetailedPlacement(rows, widths, cellX, cellY),
                    std::runtime_error);
}

BOOST_AUTO_TEST_CASE(TestAfterRow) {
  std::vector<int> widths = {2};
  std::vector<int> cellX = {5};
  std::vector<int> cellY = {10};
  std::vector<Rectangle> rows = {{-10, 6, 10, 20}};
  BOOST_CHECK_THROW(DetailedPlacement(rows, widths, cellX, cellY),
                    std::runtime_error);
}

BOOST_AUTO_TEST_CASE(TestBadY) {
  std::vector<int> widths = {2};
  std::vector<int> cellX = {5};
  std::vector<int> cellY = {11};
  std::vector<Rectangle> rows = {{-20, 20, 10, 20}};
  BOOST_CHECK_THROW(DetailedPlacement(rows, widths, cellX, cellY),
                    std::runtime_error);
}

BOOST_AUTO_TEST_CASE(TestOverlap) {
  std::vector<int> widths = {2, 2};
  std::vector<int> cellX = {5, 6};
  std::vector<int> cellY = {10, 10};
  std::vector<Rectangle> rows = {{-20, 20, 10, 20}};
  BOOST_CHECK_THROW(DetailedPlacement(rows, widths, cellX, cellY),
                    std::runtime_error);
}

BOOST_AUTO_TEST_CASE(TestChanges) {
  std::vector<int> widths = {8, 8, 8};
  std::vector<int> cellX = {0, 8, 16};
  std::vector<int> cellY = {10, 10, 10};
  std::vector<Rectangle> rows = {{0, 32, 10, 20}, {0, 32, 20, 30}};
  DetailedPlacement pl(rows, widths, cellX, cellY);
  pl.check();
  pl.unplace(0);
  pl.check();
  pl.place(0, 0, 2, 24);
  pl.check();
  pl.unplace(0);
  pl.check();
  pl.unplace(2);
  pl.check();
  pl.place(2, 1, -1, 12);
  pl.check();
  pl.place(0, 1, 2, 22);
  pl.check();
}
