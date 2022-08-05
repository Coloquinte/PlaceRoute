#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE CPP_API

#include <boost/test/unit_test.hpp>

#include "place_detailed/row_legalizer.hpp"

BOOST_AUTO_TEST_CASE(TestMiddle) {
  RowLegalizer leg(-10, 20);
  leg.push(8, -5);
  leg.push(5, 4);
  BOOST_CHECK_EQUAL(leg.getPlacement()[0], -5);
  BOOST_CHECK_EQUAL(leg.getPlacement()[1], 4);
}

BOOST_AUTO_TEST_CASE(TestConstrainLeft) {
  RowLegalizer leg(0, 20);
  leg.push(8, -5);
  BOOST_CHECK_EQUAL(leg.getPlacement()[0], 0);
}

BOOST_AUTO_TEST_CASE(TestConstrainRight) {
  RowLegalizer leg(0, 20);
  leg.push(8, 15);
  BOOST_CHECK_EQUAL(leg.getPlacement()[0], 12);
}

BOOST_AUTO_TEST_CASE(TestEquilibrium1) {
  RowLegalizer leg(-10, 20);
  leg.push(8, 4);
  leg.push(7, 4);
  BOOST_CHECK_EQUAL(leg.getPlacement()[0], 4);
  BOOST_CHECK_EQUAL(leg.getPlacement()[1], 12);
}

BOOST_AUTO_TEST_CASE(TestEquilibrium2) {
  RowLegalizer leg(-10, 20);
  leg.push(7, 4);
  leg.push(8, 4);
  BOOST_CHECK_EQUAL(leg.getPlacement()[0], -3);
  BOOST_CHECK_EQUAL(leg.getPlacement()[1], 4);
}
