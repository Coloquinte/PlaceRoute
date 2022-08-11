#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE CPP_API

#include <boost/test/unit_test.hpp>

#include "place_detailed/row_legalizer.hpp"

using namespace coloquinte;

BOOST_AUTO_TEST_CASE(TestMiddle) {
  RowLegalizer leg(-10, 20);
  leg.push(8, -5);
  leg.push(5, 4);
  BOOST_CHECK_EQUAL(leg.getPlacement().size(), 2);
  BOOST_CHECK_EQUAL(leg.getPlacement()[0], -5);
  BOOST_CHECK_EQUAL(leg.getPlacement()[1], 4);
}

BOOST_AUTO_TEST_CASE(TestConstrainLeft) {
  RowLegalizer leg(0, 20);
  leg.push(8, -5);
  BOOST_CHECK_EQUAL(leg.getPlacement().size(), 1);
  BOOST_CHECK_EQUAL(leg.getPlacement()[0], 0);
}

BOOST_AUTO_TEST_CASE(TestConstrainRight) {
  RowLegalizer leg(0, 20);
  leg.push(8, 15);
  BOOST_CHECK_EQUAL(leg.getPlacement().size(), 1);
  BOOST_CHECK_EQUAL(leg.getPlacement()[0], 12);
}

BOOST_AUTO_TEST_CASE(TestEquilibrium1) {
  RowLegalizer leg(-10, 20);
  leg.push(8, 4);
  leg.push(7, 4);
  BOOST_CHECK_EQUAL(leg.getPlacement().size(), 2);
  BOOST_CHECK_EQUAL(leg.getPlacement()[0], 4);
  BOOST_CHECK_EQUAL(leg.getPlacement()[1], 12);
}

BOOST_AUTO_TEST_CASE(TestEquilibrium2) {
  RowLegalizer leg(-10, 20);
  leg.push(7, 4);
  leg.push(8, 4);
  BOOST_CHECK_EQUAL(leg.getPlacement().size(), 2);
  BOOST_CHECK_EQUAL(leg.getPlacement()[0], -3);
  BOOST_CHECK_EQUAL(leg.getPlacement()[1], 4);
}

BOOST_AUTO_TEST_CASE(TestBothSides) {
  RowLegalizer leg(-10, 20);
  leg.push(7, -20);
  leg.push(8, 30);
  BOOST_CHECK_EQUAL(leg.getPlacement().size(), 2);
  BOOST_CHECK_EQUAL(leg.getPlacement()[0], -10);
  BOOST_CHECK_EQUAL(leg.getPlacement()[1], 12);
}

BOOST_AUTO_TEST_CASE(TestWinLeft) {
  RowLegalizer leg(-10, 20);
  leg.push(3, 50);
  leg.push(8, -5);
  BOOST_CHECK_EQUAL(leg.getPlacement().size(), 2);
  BOOST_CHECK_EQUAL(leg.getPlacement()[0], -8);
  BOOST_CHECK_EQUAL(leg.getPlacement()[1], -5);
}

BOOST_AUTO_TEST_CASE(TestWinRight) {
  RowLegalizer leg(-10, 20);
  leg.push(8, 5);
  leg.push(3, -50);
  BOOST_CHECK_EQUAL(leg.getPlacement().size(), 2);
  BOOST_CHECK_EQUAL(leg.getPlacement()[0], 5);
  BOOST_CHECK_EQUAL(leg.getPlacement()[1], 13);
}

BOOST_AUTO_TEST_CASE(TestKeepAll) {
  RowLegalizer leg(0, 100);
  for (int i = 0; i < 20; ++i) {
    leg.push(1, 2 * i);
  }
  auto pl = leg.getPlacement();
  BOOST_CHECK_EQUAL(pl.size(), 20);
  for (int i = 0; i < 20; ++i) {
    BOOST_CHECK_EQUAL(pl[i], 2 * i);
  }
}

BOOST_AUTO_TEST_CASE(TestEmpty) {
  RowLegalizer leg(-10, 20);
  BOOST_CHECK(leg.getPlacement().empty());
}