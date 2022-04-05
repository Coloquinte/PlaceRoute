#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE CPP_API

#include <boost/test/unit_test.hpp>

#include "place_global/density_legalizer.hpp"

#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>
#include <iostream>

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
    auto bisect = leg.bisect(0.0, 0.0, 1.0, 0.0, 12LL, 12L, {0, 1, 2}, LegalizationModel::L1);
    BOOST_CHECK_EQUAL(bisect.first.size(), 3);
    BOOST_CHECK_EQUAL(bisect.second.size(), 0);
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
    auto bisect = leg.bisect(0.0, 0.0, 1.0, 0.0, 12LL, 12L, {0, 1, 2}, LegalizationModel::L1);
    BOOST_CHECK_EQUAL(bisect.first.size(), 0);
    BOOST_CHECK_EQUAL(bisect.second.size(), 3);
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
    auto bisect = leg.bisect(0.0, 0.0, 0.0, 1.0, 12LL, 12L, {0, 1, 2}, LegalizationModel::L1);
    BOOST_CHECK_EQUAL(bisect.first.size(), 3);
    BOOST_CHECK_EQUAL(bisect.second.size(), 0);
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
    auto bisect = leg.bisect(0.0, 0.0, 0.0, 1.0, 12LL, 12L, {0, 1, 2}, LegalizationModel::L1);
    BOOST_CHECK_EQUAL(bisect.first.size(), 0);
    BOOST_CHECK_EQUAL(bisect.second.size(), 3);
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
    auto bisect = leg.bisect(0.0, 0.0, 1.0, 0.0, 0LL, 12L, {0, 1, 2}, LegalizationModel::L1);
    BOOST_CHECK_EQUAL(bisect.first.size(), 0);
    BOOST_CHECK_EQUAL(bisect.second.size(), 3);
}

BOOST_AUTO_TEST_CASE(NoRoom2) {
    Rectangle area(0, 1, 0, 1);
    int nbCells = 3;
    std::vector<int> cellDemand = {4, 4, 4};
    std::vector<float> cellTargetX = {0.0, 0.0, 0.0};
    std::vector<float> cellTargetY = {0.0, 0.0, 0.0};
    DensityLegalizer leg(area, nbCells);
    leg.updateCellDemand(cellDemand);
    leg.updateCellTargetX(cellTargetX);
    leg.updateCellTargetY(cellTargetY);
    auto bisect = leg.bisect(0.0, 0.0, 1.0, 0.0, 12LL, 0L, {0, 1, 2}, LegalizationModel::L1);
    BOOST_CHECK_EQUAL(bisect.first.size(), 3);
    BOOST_CHECK_EQUAL(bisect.second.size(), 0);
}


