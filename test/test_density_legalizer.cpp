#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE CPP_API

#include <boost/test/unit_test.hpp>
#include <cmath>
#include <iostream>
#include <limits>
#include <random>
#include <stdexcept>
#include <vector>

#include "place_global/density_legalizer.hpp"

BOOST_AUTO_TEST_CASE(OneSide1) {
  Rectangle area(0, 1, 0, 1);
  int nbCells = 3;
  std::vector<int> cellDemand = {4, 4, 4};
  std::vector<float> cellTargetX = {0.0, 0.0, 0.0};
  std::vector<float> cellTargetY = {0.0, 0.0, 0.0};
  DensityLegalizer leg(area, nbCells);
  leg.updateCellDemand(cellDemand);
  leg.updateCellTargetX(cellTargetX);
  leg.updateCellTargetY(cellTargetY);
  for (LegalizationModel mod :
       {LegalizationModel::L1, LegalizationModel::L2,
        LegalizationModel::L2Squared, LegalizationModel::LInf}) {
    auto bisect = leg.bisect(0.0, 0.0, 1.0, 0.0, 12LL, 12LL, {0, 1, 2}, mod);
    BOOST_CHECK_EQUAL(bisect.first.size(), 3);
    BOOST_CHECK_EQUAL(bisect.second.size(), 0);
  }
}

BOOST_AUTO_TEST_CASE(OneSide2) {
  Rectangle area(0, 1, 0, 1);
  int nbCells = 3;
  std::vector<int> cellDemand = {4, 4, 4};
  std::vector<float> cellTargetX = {1.0, 1.0, 1.0};
  std::vector<float> cellTargetY = {0.0, 0.0, 0.0};
  DensityLegalizer leg(area, nbCells);
  leg.updateCellDemand(cellDemand);
  leg.updateCellTargetX(cellTargetX);
  leg.updateCellTargetY(cellTargetY);
  for (LegalizationModel mod :
       {LegalizationModel::L1, LegalizationModel::L2,
        LegalizationModel::L2Squared, LegalizationModel::LInf}) {
    auto bisect = leg.bisect(0.0, 0.0, 1.0, 0.0, 12LL, 12LL, {0, 1, 2}, mod);
    BOOST_CHECK_EQUAL(bisect.first.size(), 0);
    BOOST_CHECK_EQUAL(bisect.second.size(), 3);
  }
}

BOOST_AUTO_TEST_CASE(OneSide3) {
  Rectangle area(0, 1, 0, 1);
  int nbCells = 3;
  std::vector<int> cellDemand = {4, 4, 4};
  std::vector<float> cellTargetX = {0.0, 0.0, 0.0};
  std::vector<float> cellTargetY = {0.0, 0.0, 0.0};
  DensityLegalizer leg(area, nbCells);
  leg.updateCellDemand(cellDemand);
  leg.updateCellTargetX(cellTargetX);
  leg.updateCellTargetY(cellTargetY);
  for (LegalizationModel mod :
       {LegalizationModel::L1, LegalizationModel::L2,
        LegalizationModel::L2Squared, LegalizationModel::LInf}) {
    auto bisect = leg.bisect(0.0, 0.0, 0.0, 1.0, 12LL, 12LL, {0, 1, 2}, mod);
    BOOST_CHECK_EQUAL(bisect.first.size(), 3);
    BOOST_CHECK_EQUAL(bisect.second.size(), 0);
  }
}

BOOST_AUTO_TEST_CASE(OneSide4) {
  Rectangle area(0, 1, 0, 1);
  int nbCells = 3;
  std::vector<int> cellDemand = {4, 4, 4};
  std::vector<float> cellTargetX = {0.0, 0.0, 0.0};
  std::vector<float> cellTargetY = {1.0, 1.0, 1.0};
  DensityLegalizer leg(area, nbCells);
  leg.updateCellDemand(cellDemand);
  leg.updateCellTargetX(cellTargetX);
  leg.updateCellTargetY(cellTargetY);
  for (LegalizationModel mod :
       {LegalizationModel::L1, LegalizationModel::L2,
        LegalizationModel::L2Squared, LegalizationModel::LInf}) {
    auto bisect = leg.bisect(0.0, 0.0, 0.0, 1.0, 12LL, 12LL, {0, 1, 2}, mod);
    BOOST_CHECK_EQUAL(bisect.first.size(), 0);
    BOOST_CHECK_EQUAL(bisect.second.size(), 3);
  }
}

BOOST_AUTO_TEST_CASE(NoRoom1) {
  Rectangle area(0, 1, 0, 1);
  int nbCells = 3;
  std::vector<int> cellDemand = {4, 4, 4};
  std::vector<float> cellTargetX = {0.0, 0.0, 0.0};
  std::vector<float> cellTargetY = {0.0, 0.0, 0.0};
  DensityLegalizer leg(area, nbCells);
  leg.updateCellDemand(cellDemand);
  leg.updateCellTargetX(cellTargetX);
  leg.updateCellTargetY(cellTargetY);
  for (LegalizationModel mod :
       {LegalizationModel::L1, LegalizationModel::L2,
        LegalizationModel::L2Squared, LegalizationModel::LInf}) {
    auto bisect = leg.bisect(0.0, 0.0, 1.0, 0.0, 0LL, 12LL, {0, 1, 2}, mod);
    BOOST_CHECK_EQUAL(bisect.first.size(), 0);
    BOOST_CHECK_EQUAL(bisect.second.size(), 3);
  }
}

BOOST_AUTO_TEST_CASE(NoRoom2) {
  Rectangle area(0, 1, 0, 1);
  int nbCells = 3;
  std::vector<int> cellDemand = {4, 4, 4};
  std::vector<float> cellTargetX = {1.0, 1.0, 1.0};
  std::vector<float> cellTargetY = {0.0, 0.0, 0.0};
  DensityLegalizer leg(area, nbCells);
  leg.updateCellDemand(cellDemand);
  leg.updateCellTargetX(cellTargetX);
  leg.updateCellTargetY(cellTargetY);
  for (LegalizationModel mod :
       {LegalizationModel::L1, LegalizationModel::L2,
        LegalizationModel::L2Squared, LegalizationModel::LInf}) {
    auto bisect = leg.bisect(0.0, 0.0, 1.0, 0.0, 12LL, 0L, {0, 1, 2}, mod);
    BOOST_CHECK_EQUAL(bisect.first.size(), 3);
    BOOST_CHECK_EQUAL(bisect.second.size(), 0);
  }
}

BOOST_AUTO_TEST_CASE(Sorted1) {
  Rectangle area(0, 1, 0, 1);
  int nbCells = 4;
  std::vector<int> cellDemand = {4, 4, 4, 4};
  std::vector<float> cellTargetX = {0.0, 0.0, 1.0, 1.0};
  std::vector<float> cellTargetY = {0.0, 0.0, 0.0, 0.0};
  DensityLegalizer leg(area, nbCells);
  leg.updateCellDemand(cellDemand);
  leg.updateCellTargetX(cellTargetX);
  leg.updateCellTargetY(cellTargetY);
  for (LegalizationModel mod :
       {LegalizationModel::L1, LegalizationModel::L2,
        LegalizationModel::L2Squared, LegalizationModel::LInf}) {
    auto bisect = leg.bisect(0.0, 0.0, 1.0, 0.0, 16LL, 16LL, {0, 1, 2, 3}, mod);
    BOOST_CHECK_EQUAL(bisect.first.size(), 2);
    BOOST_CHECK_EQUAL(bisect.second.size(), 2);
    BOOST_CHECK(bisect.first[0] == 0 || bisect.first[1] == 0);
    BOOST_CHECK(bisect.first[0] == 1 || bisect.first[1] == 1);
    BOOST_CHECK(bisect.second[0] == 2 || bisect.second[1] == 2);
    BOOST_CHECK(bisect.second[0] == 3 || bisect.second[1] == 3);
  }
}

BOOST_AUTO_TEST_CASE(Sorted2) {
  Rectangle area(0, 1, 0, 1);
  int nbCells = 4;
  std::vector<int> cellDemand = {4, 4, 4, 4};
  std::vector<float> cellTargetX = {1.0, 1.0, 0.0, 0.0};
  std::vector<float> cellTargetY = {0.0, 0.0, 0.0, 0.0};
  DensityLegalizer leg(area, nbCells);
  leg.updateCellDemand(cellDemand);
  leg.updateCellTargetX(cellTargetX);
  leg.updateCellTargetY(cellTargetY);
  for (LegalizationModel mod :
       {LegalizationModel::L1, LegalizationModel::L2,
        LegalizationModel::L2Squared, LegalizationModel::LInf}) {
    auto bisect = leg.bisect(0.0, 0.0, 1.0, 0.0, 16LL, 16LL, {0, 1, 2, 3}, mod);
    BOOST_CHECK_EQUAL(bisect.first.size(), 2);
    BOOST_CHECK_EQUAL(bisect.second.size(), 2);
    BOOST_CHECK(bisect.second[0] == 0 || bisect.second[1] == 0);
    BOOST_CHECK(bisect.second[0] == 1 || bisect.second[1] == 1);
    BOOST_CHECK(bisect.first[0] == 2 || bisect.first[1] == 2);
    BOOST_CHECK(bisect.first[0] == 3 || bisect.first[1] == 3);
  }
}

BOOST_AUTO_TEST_CASE(LeastOverflow1) {
  Rectangle area(0, 1, 0, 1);
  int nbCells = 3;
  std::vector<int> cellDemand = {4, 4, 4};
  std::vector<float> cellTargetX = {0.0, 0.0, 0.0};
  std::vector<float> cellTargetY = {0.0, 0.0, 0.0};
  DensityLegalizer leg(area, nbCells);
  leg.updateCellDemand(cellDemand);
  leg.updateCellTargetX(cellTargetX);
  leg.updateCellTargetY(cellTargetY);
  auto bisect = leg.bisect(0.0, 0.0, 1.0, 0.0, 6LL, 7LL, {0, 1, 2},
                           LegalizationModel::L1);
  BOOST_CHECK_EQUAL(bisect.first.size(), 1);
  BOOST_CHECK_EQUAL(bisect.second.size(), 2);
}

BOOST_AUTO_TEST_CASE(LeastOverflow2) {
  Rectangle area(0, 1, 0, 1);
  int nbCells = 4;
  std::vector<int> cellDemand = {4, 4, 4, 4};
  std::vector<float> cellTargetX = {0.0, 0.0, 0.0, 0.0};
  std::vector<float> cellTargetY = {0.0, 0.0, 0.0, 0.0};
  DensityLegalizer leg(area, nbCells);
  leg.updateCellDemand(cellDemand);
  leg.updateCellTargetX(cellTargetX);
  leg.updateCellTargetY(cellTargetY);
  auto bisect = leg.bisect(0.0, 0.0, 1.0, 0.0, 6LL, 7LL, {0, 1, 2, 3},
                           LegalizationModel::L1);
  BOOST_CHECK_EQUAL(bisect.first.size(), 2);
  BOOST_CHECK_EQUAL(bisect.second.size(), 2);
}

BOOST_AUTO_TEST_CASE(LeastOverflow3) {
  Rectangle area(0, 1, 0, 1);
  int nbCells = 1;
  std::vector<int> cellDemand = {4};
  std::vector<float> cellTargetX = {0.0};
  std::vector<float> cellTargetY = {0.0};
  DensityLegalizer leg(area, nbCells);
  leg.updateCellDemand(cellDemand);
  leg.updateCellTargetX(cellTargetX);
  leg.updateCellTargetY(cellTargetY);
  auto bisect =
      leg.bisect(0.0, 0.0, 1.0, 0.0, 1LL, 2LL, {0}, LegalizationModel::L1);
  BOOST_CHECK_EQUAL(bisect.first.size(), 0);
  BOOST_CHECK_EQUAL(bisect.second.size(), 1);
}

BOOST_AUTO_TEST_CASE(LargeAssign) {
  Rectangle area(-100, 50, 100, 200);
  int nbCells = 10000;
  int nbBinsX = 20;
  int nbBinsY = 10;
  std::vector<int> cellDemand;
  std::vector<float> cellTargetX;
  std::vector<float> cellTargetY;
  std::mt19937 rgen;
  for (int i = 0; i < nbCells; ++i) {
    cellDemand.push_back(i);
    cellTargetX.push_back(
        std::uniform_real_distribution<float>(area.minX, area.maxX)(rgen));
    cellTargetY.push_back(
        std::uniform_real_distribution<float>(area.minY, area.maxY)(rgen));
  }
  DensityLegalizer leg(area, nbCells);
  leg.updateBins(nbBinsX, nbBinsY);
  leg.updateCellDemand(cellDemand);
  leg.updateCellTargetX(cellTargetX);
  leg.updateCellTargetY(cellTargetY);
  leg.assign();
}

BOOST_AUTO_TEST_CASE(LargeBisect) {
  Rectangle area(-100, 50, 100, 200);
  int nbCells = 10000;
  int nbBinsX = 19;
  int nbBinsY = 11;
  std::vector<int> cellDemand;
  std::vector<float> cellTargetX;
  std::vector<float> cellTargetY;
  std::mt19937 rgen;
  for (int i = 0; i < nbCells; ++i) {
    cellDemand.push_back(1);
    cellTargetX.push_back(
        std::uniform_real_distribution<float>(area.minX, area.maxX)(rgen));
    cellTargetY.push_back(
        std::uniform_real_distribution<float>(area.minY, area.maxY)(rgen));
  }
  DensityLegalizer leg(area, nbCells);
  leg.updateBins(nbBinsX, nbBinsY);
  leg.updateCellDemand(cellDemand);
  leg.updateCellTargetX(cellTargetX);
  leg.updateCellTargetY(cellTargetY);
  leg.bisect(LegalizationModel::L2);
}
