#define BOOST_TEST_MODULE ROW_LEGALIZER

#include <boost/test/unit_test.hpp>

#include "place_detailed/row_legalizer.hpp"

using namespace coloquinte;

BOOST_AUTO_TEST_CASE(TestMiddle) {
  RowLegalizer leg(-10, 20);
  int c1 = leg.push(8, -5);
  int c2 = leg.push(5, 4);
  BOOST_CHECK_EQUAL(leg.getPlacement().size(), 2);
  BOOST_CHECK_EQUAL(leg.getPlacement()[0], -5);
  BOOST_CHECK_EQUAL(leg.getPlacement()[1], 4);
  BOOST_CHECK_EQUAL(c1, 0);
  BOOST_CHECK_EQUAL(c2, 0);
}

BOOST_AUTO_TEST_CASE(TestConstrainLeft) {
  RowLegalizer leg(0, 20);
  int c1 = leg.push(8, -5);
  BOOST_CHECK_EQUAL(leg.getPlacement().size(), 1);
  BOOST_CHECK_EQUAL(leg.getPlacement()[0], 0);
  BOOST_CHECK_EQUAL(c1, 40);
}

BOOST_AUTO_TEST_CASE(TestConstrainRight) {
  RowLegalizer leg(0, 20);
  int c1 = leg.push(8, 15);
  BOOST_CHECK_EQUAL(leg.getPlacement().size(), 1);
  BOOST_CHECK_EQUAL(leg.getPlacement()[0], 12);
  BOOST_CHECK_EQUAL(c1, 24);
}

BOOST_AUTO_TEST_CASE(TestEquilibrium1) {
  RowLegalizer leg(-10, 20);
  int c1 = leg.push(8, 4);
  int c2 = leg.push(7, 4);
  BOOST_CHECK_EQUAL(leg.getPlacement().size(), 2);
  BOOST_CHECK_EQUAL(leg.getPlacement()[0], 4);
  BOOST_CHECK_EQUAL(leg.getPlacement()[1], 12);
  BOOST_CHECK_EQUAL(c1, 0);
  BOOST_CHECK_EQUAL(c2, 56);
}

BOOST_AUTO_TEST_CASE(TestEquilibrium2) {
  RowLegalizer leg(-10, 20);
  int c1 = leg.push(7, 4);
  int c2 = leg.push(8, 4);
  BOOST_CHECK_EQUAL(leg.getPlacement().size(), 2);
  BOOST_CHECK_EQUAL(leg.getPlacement()[0], -3);
  BOOST_CHECK_EQUAL(leg.getPlacement()[1], 4);
  BOOST_CHECK_EQUAL(c1, 0);
  BOOST_CHECK_EQUAL(c2, 49);
}

BOOST_AUTO_TEST_CASE(TestBothSides) {
  RowLegalizer leg(-10, 20);
  int c1 = leg.push(7, -20);
  int c2 = leg.push(8, 30);
  BOOST_CHECK_EQUAL(leg.getPlacement().size(), 2);
  BOOST_CHECK_EQUAL(leg.getPlacement()[0], -10);
  BOOST_CHECK_EQUAL(leg.getPlacement()[1], 12);
  BOOST_CHECK_EQUAL(c1, 70);
  BOOST_CHECK_EQUAL(c2, 144);
}

BOOST_AUTO_TEST_CASE(TestWinLeft) {
  RowLegalizer leg(-10, 20);
  int c1 = leg.push(3, 50);
  int c2 = leg.push(8, -5);
  BOOST_CHECK_EQUAL(leg.getPlacement().size(), 2);
  BOOST_CHECK_EQUAL(leg.getPlacement()[0], -8);
  BOOST_CHECK_EQUAL(leg.getPlacement()[1], -5);
  BOOST_CHECK_EQUAL(c1, 99);
  BOOST_CHECK_EQUAL(c2, 75);
}

BOOST_AUTO_TEST_CASE(TestWinRight) {
  RowLegalizer leg(-10, 20);
  int c1 = leg.push(8, 5);
  int c2 = leg.push(3, -50);
  BOOST_CHECK_EQUAL(leg.getPlacement().size(), 2);
  BOOST_CHECK_EQUAL(leg.getPlacement()[0], 5);
  BOOST_CHECK_EQUAL(leg.getPlacement()[1], 13);
  BOOST_CHECK_EQUAL(c1, 0);
  BOOST_CHECK_EQUAL(c2, 189);
}

BOOST_AUTO_TEST_CASE(TestKeepAll) {
  RowLegalizer leg(0, 100);
  for (int i = 0; i < 20; ++i) {
    int c = leg.push(1, 2 * i);
    BOOST_CHECK_EQUAL(c, 0);
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

int computeCost(const std::vector<int>& widths, const std::vector<int>& targets,
                const std::vector<int>& positions) {
  assert(widths.size() == targets.size());
  assert(widths.size() == positions.size());
  int ret = 0;
  for (int i = 0; i < widths.size(); ++i) {
    ret += widths[i] * std::abs(targets[i] - positions[i]);
  }
  return ret;
}

void fuzzRowLegalizer(int begin, int end, const std::vector<int>& widths,
                      const std::vector<int>& targets) {
  assert(widths.size() == targets.size());

  RowLegalizer leg(begin, end);
  int cumCost = 0;
  for (int i = 0; i < widths.size(); ++i) {
    int predictedCost = leg.getCost(widths[i], targets[i]);
    int cost = leg.push(widths[i], targets[i]);
    BOOST_CHECK_EQUAL(predictedCost, cost);
    cumCost += cost;
    std::vector<int> curWidths(widths.begin(), widths.begin() + i + 1);
    std::vector<int> curTargets(targets.begin(), targets.begin() + i + 1);
    BOOST_CHECK_EQUAL(cumCost,
                      computeCost(curWidths, curTargets, leg.getPlacement()));
    leg.check();
  }
}

BOOST_AUTO_TEST_CASE(TestComplex1) {
  std::vector<int> widths = {5, 5, 8, 7, 6, 3};
  std::vector<int> targets = {20, 10, 45, 67, 32, 49};
  fuzzRowLegalizer(0, 100, widths, targets);
}

BOOST_AUTO_TEST_CASE(TestComplex2) {
  std::vector<int> widths = {5, 5, 8, 7, 6, 3, 23, 4, 3, 6, 7, 3};
  std::vector<int> targets = {20, -10, -45, 67,  32,  49,
                              20, 80,  36,  120, 110, 90};
  fuzzRowLegalizer(0, 100, widths, targets);
}