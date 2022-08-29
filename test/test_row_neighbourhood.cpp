#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE CPP_API

#include <boost/test/unit_test.hpp>

#include "place_detailed/row_neighbourhood.hpp"

using namespace coloquinte;

BOOST_AUTO_TEST_CASE(TestBasic1) {
  std::vector<Rectangle> rows = {{0, 100, 0, 12}, {10, 90, 12, 24}};
  RowNeighbourhood neigh(rows);

  BOOST_CHECK_EQUAL(neigh.rowsAbove(0).size(), 1);
  BOOST_CHECK_EQUAL(neigh.rowsBelow(0).size(), 0);
  BOOST_CHECK_EQUAL(neigh.rowsLeft(0).size(), 0);
  BOOST_CHECK_EQUAL(neigh.rowsRight(0).size(), 0);
  BOOST_CHECK_EQUAL(neigh.rowsAbove(1).size(), 0);
  BOOST_CHECK_EQUAL(neigh.rowsBelow(1).size(), 1);
  BOOST_CHECK_EQUAL(neigh.rowsLeft(1).size(), 0);
  BOOST_CHECK_EQUAL(neigh.rowsRight(1).size(), 0);

  BOOST_CHECK_EQUAL(neigh.rowsAbove(0)[0], 1);
  BOOST_CHECK_EQUAL(neigh.rowsBelow(1)[0], 0);
}

BOOST_AUTO_TEST_CASE(TestBasic2) {
  std::vector<Rectangle> rows = {{0, 100, 0, 12},
                                 {10, 90, 12, 24},
                                 {-200, -10, 12, 24},
                                 {-100, 5, 24, 36}};
  RowNeighbourhood neigh(rows, 3);

  BOOST_CHECK_EQUAL(neigh.rowsAbove(0).size(), 2);
  BOOST_CHECK_EQUAL(neigh.rowsBelow(0).size(), 0);
  BOOST_CHECK_EQUAL(neigh.rowsLeft(0).size(), 1);
  BOOST_CHECK_EQUAL(neigh.rowsRight(0).size(), 0);

  BOOST_CHECK_EQUAL(neigh.rowsAbove(1).size(), 0);
  BOOST_CHECK_EQUAL(neigh.rowsBelow(1).size(), 1);
  BOOST_CHECK_EQUAL(neigh.rowsLeft(1).size(), 2);
  BOOST_CHECK_EQUAL(neigh.rowsRight(1).size(), 0);

  BOOST_CHECK_EQUAL(neigh.rowsAbove(2).size(), 1);
  BOOST_CHECK_EQUAL(neigh.rowsBelow(2).size(), 0);
  BOOST_CHECK_EQUAL(neigh.rowsLeft(2).size(), 0);
  BOOST_CHECK_EQUAL(neigh.rowsRight(2).size(), 2);

  BOOST_CHECK_EQUAL(neigh.rowsAbove(3).size(), 0);
  BOOST_CHECK_EQUAL(neigh.rowsBelow(3).size(), 2);
  BOOST_CHECK_EQUAL(neigh.rowsLeft(3).size(), 0);
  BOOST_CHECK_EQUAL(neigh.rowsRight(3).size(), 1);
}

/**
 * @brief Only valid if the rows don't overlap
 */
void checkNeighbourhood(const RowNeighbourhood &neigh, int nbRows) {
  BOOST_CHECK_EQUAL(neigh.nbRows(), nbRows);
  for (int i = 0; i < nbRows; ++i) {
    int nb = 0;
    nb += neigh.rowsBelow(i).size();
    nb += neigh.rowsAbove(i).size();
    nb += neigh.rowsLeft(i).size();
    nb += neigh.rowsRight(i).size();
    BOOST_CHECK_EQUAL(nb, nbRows - 1);
  }
}

BOOST_AUTO_TEST_CASE(TestAbove) {
  std::vector<Rectangle> rows = {{0, 100, 0, 12},
                                 {10, 90, 12, 24},
                                 {-200, 1, 12, 24},
                                 {-100, 5, 24, 36},
                                 {99, 105, 10, 22}};
  int nbRows = rows.size();
  RowNeighbourhood neigh(rows, nbRows);
  checkNeighbourhood(neigh, nbRows);
  BOOST_CHECK_EQUAL(neigh.rowsAbove(0).size(), nbRows - 1);
}

BOOST_AUTO_TEST_CASE(TestBelow) {
  std::vector<Rectangle> rows = {{0, 100, 0, 12},
                                 {10, 90, -12, 0},
                                 {-200, 1, -12, 0},
                                 {-100, 5, -36, 0},
                                 {99, 105, -24, 12}};
  int nbRows = rows.size();
  RowNeighbourhood neigh(rows, nbRows);
  checkNeighbourhood(neigh, nbRows);
  BOOST_CHECK_EQUAL(neigh.rowsBelow(0).size(), nbRows - 1);
}

BOOST_AUTO_TEST_CASE(TestLeft) {
  std::vector<Rectangle> rows = {{0, 100, 0, 12},
                                 {-10, 0, -12, 0},
                                 {-100, -10, -12, 0},
                                 {-90, -1, 0, 12},
                                 {-20, -5, 24, 36}};
  int nbRows = rows.size();
  RowNeighbourhood neigh(rows, nbRows);
  checkNeighbourhood(neigh, nbRows);
  BOOST_CHECK_EQUAL(neigh.rowsLeft(0).size(), nbRows - 1);
}

BOOST_AUTO_TEST_CASE(TestRight) {
  std::vector<Rectangle> rows = {{0, 100, 0, 12},
                                 {110, 130, -12, 0},
                                 {130, 300, -12, 0},
                                 {100, 201, -36, -24},
                                 {104, 105, 0, 12}};
  int nbRows = rows.size();
  RowNeighbourhood neigh(rows, nbRows);
  checkNeighbourhood(neigh, nbRows);
  BOOST_CHECK_EQUAL(neigh.rowsRight(0).size(), nbRows - 1);
}