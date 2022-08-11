#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE CPP_API

#include <boost/test/unit_test.hpp>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <vector>

#include "place_global/net_model.hpp"

using namespace coloquinte;

BOOST_AUTO_TEST_CASE(b2b_singlepin) {
  int nbCells = 4;
  NetModel model(nbCells);
  model.addNet({0}, {0.0f}, 4.0f, 4.0f);
  model.addNet({1}, {5.0f}, 1.0f, 1.0f);
  model.addNet({2}, {-6.0f}, 2.0f, 2.0f);
  model.addNet({3}, {2.0f}, 3.0f, 3.0f);
  BOOST_CHECK_EQUAL(model.nbCells(), 4);
  BOOST_CHECK_EQUAL(model.nbNets(), 4);
  BOOST_CHECK_EQUAL(model.nbPins(0), 2);
  BOOST_CHECK_EQUAL(model.nbPins(1), 2);
  BOOST_CHECK_EQUAL(model.nbPins(2), 2);
  BOOST_CHECK_EQUAL(model.nbPins(3), 2);
  std::vector<float> place = {10.0f, 20.0f, 30.0f, 40.0f};
  auto res = model.solveB2B(place, 1.0e-8);
  BOOST_CHECK_CLOSE(res[0], 4.0, 0.001);
  BOOST_CHECK_CLOSE(res[1], -4.0, 0.001);
  BOOST_CHECK_CLOSE(res[2], 8.0, 0.001);
  BOOST_CHECK_CLOSE(res[3], 1.0, 0.001);
}
BOOST_AUTO_TEST_CASE(b2b_mid) {
  int nbCells = 2;
  NetModel model(nbCells);
  model.addNet({0, 1}, {0.0f, 0.0f}, 0.0f, 4.0f);
  std::vector<float> place = {2.0f, 3.0f};
  auto res = model.solveB2B(place, 1.0e-8);
  BOOST_CHECK_CLOSE(res[0], 2.0, 0.001);
  BOOST_CHECK_CLOSE(res[1], 3.0, 0.001);
}

BOOST_AUTO_TEST_CASE(b2b_two_pins) {
  int nbCells = 1;
  NetModel model(nbCells);
  model.addNet({0}, {0.0f}, 0.0f, 4.0f);
  std::vector<float> place = {8.0f};
  auto res = model.solveB2B(place, 1.0e-8);
  BOOST_CHECK_CLOSE(res[0], 8.0 / 3, 0.001);
}

BOOST_AUTO_TEST_CASE(star_singlepin) {
  int nbCells = 4;
  NetModel model(nbCells);
  model.addNet({0}, {0.0f}, 4.0f, 4.0f);
  model.addNet({1}, {5.0f}, 1.0f, 1.0f);
  model.addNet({2}, {-6.0f}, 2.0f, 2.0f);
  model.addNet({3}, {2.0f}, 3.0f, 3.0f);
  std::vector<float> place = {10.0f, 20.0f, 30.0f, 40.0f};
  auto res = model.solveStar(place, 1.0e-8);
  BOOST_CHECK_CLOSE(res[0], 4.0, 0.001);
  BOOST_CHECK_CLOSE(res[1], -4.0, 0.001);
  BOOST_CHECK_CLOSE(res[2], 8.0, 0.001);
  BOOST_CHECK_CLOSE(res[3], 1.0, 0.001);
}

BOOST_AUTO_TEST_CASE(star_mid) {
  int nbCells = 2;
  NetModel model(nbCells);
  model.addNet({0, 1}, {0.0f, 0.0f}, 0.0f, 4.0f);
  std::vector<float> place = {2.0f, 3.0f};
  auto res = model.solveStar(place, 1.0e-8);
  BOOST_CHECK_CLOSE(res[0], 2.0, 0.001);
  BOOST_CHECK_CLOSE(res[1], 3.0, 0.001);
}