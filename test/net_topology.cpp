#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE CPP_API

#include <boost/test/unit_test.hpp>

#include "place_global/topology.hpp"

#include <xtensor/xio.hpp>

#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>
#include <iostream>

BOOST_AUTO_TEST_CASE(valueHPWLTwoPins) {
    int nbCells = 2;
    xt::xtensor<int, 2> pinCells = {{0, 1}};
    xt::xtensor<float, 2> pinOffsets = {{0.0, 0.0}};
    xt::xtensor<float, 1> place = {-1000.0, 1000.0};
    auto topo = NetTopologyFixedSize(nbCells, pinCells, pinOffsets);
    BOOST_CHECK_CLOSE(topo.valueHPWL(place), 2000.0, 0.001);
}

BOOST_AUTO_TEST_CASE(valueLSETwoPins) {
    int nbCells = 2;
    xt::xtensor<int, 2> pinCells = {{0, 1}};
    xt::xtensor<float, 2> pinOffsets = {{0.0, 0.0}};
    xt::xtensor<float, 1> place = {-1000.0, 1000.0};
    auto topo = NetTopologyFixedSize(nbCells, pinCells, pinOffsets);
    BOOST_CHECK_CLOSE(topo.valueLSE(place, 0.01), 2000.0, 0.001);
}

BOOST_AUTO_TEST_CASE(valueWATwoPins) {
    int nbCells = 2;
    xt::xtensor<int, 2> pinCells = {{0, 1}};
    xt::xtensor<float, 2> pinOffsets = {{0.0, 0.0}};
    xt::xtensor<float, 1> place = {-1000.0, 1000.0};
    auto topo = NetTopologyFixedSize(nbCells, pinCells, pinOffsets);
    BOOST_CHECK_CLOSE(topo.valueWA(place, 0.01), 2000.0, 0.001);
}

BOOST_AUTO_TEST_CASE(valueHPWLThreePins) {
    int nbCells = 3;
    xt::xtensor<int, 2> pinCells = {{2, 0, 1}};
    xt::xtensor<float, 2> pinOffsets = {{0.0, 0.0, 0.0}};
    xt::xtensor<float, 1> place = {0.0, -1000.0, 1000.0};
    auto topo = NetTopologyFixedSize(nbCells, pinCells, pinOffsets);
    BOOST_CHECK_CLOSE(topo.valueHPWL(place), 2000.0, 0.001);
}

BOOST_AUTO_TEST_CASE(valueLSEThreePins) {
    int nbCells = 3;
    xt::xtensor<int, 2> pinCells = {{2, 0, 1}};
    xt::xtensor<float, 2> pinOffsets = {{0.0, 0.0, 0.0}};
    xt::xtensor<float, 1> place = {0.0, -1000.0, 1000.0};
    auto topo = NetTopologyFixedSize(nbCells, pinCells, pinOffsets);
    BOOST_CHECK_CLOSE(topo.valueLSE(place, 0.01), 2000.0, 0.001);
}

BOOST_AUTO_TEST_CASE(valueWAThreePins) {
    int nbCells = 3;
    xt::xtensor<int, 2> pinCells = {{2, 0, 1}};
    xt::xtensor<float, 2> pinOffsets = {{0.0, 0.0, 0.0}};
    xt::xtensor<float, 1> place = {0.0, -1000.0, 1000.0};
    auto topo = NetTopologyFixedSize(nbCells, pinCells, pinOffsets);
    BOOST_CHECK_CLOSE(topo.valueWA(place, 0.01), 2000.0, 0.001);
}

BOOST_AUTO_TEST_CASE(valueHPWLRepeat) {
    int nbCells = 1;
    xt::xtensor<int, 2> pinCells = {{0, 0}};
    xt::xtensor<float, 2> pinOffsets = {{1000.0, -1000.0}};
    xt::xtensor<float, 1> place = {5000.0};
    auto topo = NetTopologyFixedSize(nbCells, pinCells, pinOffsets);
    BOOST_CHECK_CLOSE(topo.valueHPWL(place), 2000.0, 0.001);
}

BOOST_AUTO_TEST_CASE(valueHPWLMultiNet) {
    int nbCells = 3;
    xt::xtensor<int, 2> pinCells = {{0, 1}, {2, 1}, {0, 2}};
    xt::xtensor<float, 2> pinOffsets = {{0.0, 0.0}, {0.0, 0.0}, {-5.0, 5.0}};
    xt::xtensor<float, 1> place = {0.0, 2.0, 3.0};
    auto topo = NetTopologyFixedSize(nbCells, pinCells, pinOffsets);
    BOOST_CHECK_CLOSE(topo.valueHPWL(place), 16.0, 0.001);
}

BOOST_AUTO_TEST_CASE(gradHPWLTwoPins) {
    int nbCells = 2;
    xt::xtensor<int, 2> pinCells = {{0, 1}};
    xt::xtensor<float, 2> pinOffsets = {{0.0, 0.0}};
    xt::xtensor<float, 1> place = {-1000.0, 1000.0};
    auto topo = NetTopologyFixedSize(nbCells, pinCells, pinOffsets);
    xt::xtensor<float, 1> grad = topo.gradHPWL(place);
    BOOST_CHECK_CLOSE(grad[0], -1.0, 0.001);
    BOOST_CHECK_CLOSE(grad[1], 1.0, 0.001);
}

BOOST_AUTO_TEST_CASE(gradLSETwoPins) {
    int nbCells = 2;
    xt::xtensor<int, 2> pinCells = {{0, 1}};
    xt::xtensor<float, 2> pinOffsets = {{0.0, 0.0}};
    xt::xtensor<float, 1> place = {-1000.0, 1000.0};
    auto topo = NetTopologyFixedSize(nbCells, pinCells, pinOffsets);
    xt::xtensor<float, 1> grad = topo.gradLSE(place, 0.01);
    BOOST_CHECK_CLOSE(grad[0], -1.0, 0.001);
    BOOST_CHECK_CLOSE(grad[1], 1.0, 0.001);
}

BOOST_AUTO_TEST_CASE(gradWATwoPins) {
    int nbCells = 2;
    xt::xtensor<int, 2> pinCells = {{0, 1}};
    xt::xtensor<float, 2> pinOffsets = {{0.0, 0.0}};
    xt::xtensor<float, 1> place = {-1000.0, 1000.0};
    auto topo = NetTopologyFixedSize(nbCells, pinCells, pinOffsets);
    xt::xtensor<float, 1> grad = topo.gradWA(place, 0.01);
    BOOST_CHECK_CLOSE(grad[0], -1.0, 0.001);
    BOOST_CHECK_CLOSE(grad[1], 1.0, 0.001);
}

BOOST_AUTO_TEST_CASE(gradHPWLThreePins) {
    int nbCells = 3;
    xt::xtensor<int, 2> pinCells = {{2, 0, 1}};
    xt::xtensor<float, 2> pinOffsets = {{0.0, 0.0, 0.0}};
    xt::xtensor<float, 1> place = {0.0, -1000.0, 1000.0};
    auto topo = NetTopologyFixedSize(nbCells, pinCells, pinOffsets);
    xt::xtensor<float, 1> grad = topo.gradHPWL(place);
    BOOST_CHECK_SMALL(grad[0], 0.001f);
    BOOST_CHECK_CLOSE(grad[1], -1.0, 0.001);
    BOOST_CHECK_CLOSE(grad[2], 1.0, 0.001);
}

BOOST_AUTO_TEST_CASE(gradLSEThreePins) {
    int nbCells = 3;
    xt::xtensor<int, 2> pinCells = {{2, 0, 1}};
    xt::xtensor<float, 2> pinOffsets = {{0.0, 0.0, 0.0}};
    xt::xtensor<float, 1> place = {0.0, -1000.0, 1000.0};
    auto topo = NetTopologyFixedSize(nbCells, pinCells, pinOffsets);
    xt::xtensor<float, 1> grad = topo.gradLSE(place, 0.01);
    BOOST_CHECK_SMALL(grad[0], 0.001f);
    BOOST_CHECK_CLOSE(grad[1], -1.0, 0.001);
    BOOST_CHECK_CLOSE(grad[2], 1.0, 0.001);
}

BOOST_AUTO_TEST_CASE(gradWAThreePins) {
    int nbCells = 3;
    xt::xtensor<int, 2> pinCells = {{2, 0, 1}};
    xt::xtensor<float, 2> pinOffsets = {{0.0, 0.0, 0.0}};
    xt::xtensor<float, 1> place = {0.0, -1000.0, 1000.0};
    auto topo = NetTopologyFixedSize(nbCells, pinCells, pinOffsets);
    xt::xtensor<float, 1> grad = topo.gradWA(place, 0.01);
    BOOST_CHECK_SMALL(grad[0], 0.001f);
    BOOST_CHECK_CLOSE(grad[1], -1.0, 0.001);
    BOOST_CHECK_CLOSE(grad[2], 1.0, 0.001);
}

BOOST_AUTO_TEST_CASE(gradHPWLMultiNet) {
    int nbCells = 3;
    xt::xtensor<int, 2> pinCells = {{0, 1}, {2, 1}, {0, 2}};
    xt::xtensor<float, 2> pinOffsets = {{0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}};
    xt::xtensor<float, 1> place = {0.0, 1000.0, 2000.0};
    auto topo = NetTopologyFixedSize(nbCells, pinCells, pinOffsets);
    xt::xtensor<float, 1> grad = topo.gradHPWL(place);
    BOOST_CHECK_CLOSE(grad[0], -2.0, 0.001);
    BOOST_CHECK_SMALL(grad[1], 0.001f);
    BOOST_CHECK_CLOSE(grad[2], 2.0, 0.001);
}

BOOST_AUTO_TEST_CASE(gradLSEMultiNet) {
    int nbCells = 3;
    xt::xtensor<int, 2> pinCells = {{0, 1}, {2, 1}, {0, 2}};
    xt::xtensor<float, 2> pinOffsets = {{0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}};
    xt::xtensor<float, 1> place = {0.0, 1000.0, 2000.0};
    auto topo = NetTopologyFixedSize(nbCells, pinCells, pinOffsets);
    xt::xtensor<float, 1> grad = topo.gradLSE(place, 0.01);
    BOOST_CHECK_CLOSE(grad[0], -2.0, 0.001);
    BOOST_CHECK_SMALL(grad[1], 0.001f);
    BOOST_CHECK_CLOSE(grad[2], 2.0, 0.001);
}

BOOST_AUTO_TEST_CASE(gradWAMultiNet) {
    int nbCells = 3;
    xt::xtensor<int, 2> pinCells = {{0, 1}, {2, 1}, {0, 2}};
    xt::xtensor<float, 2> pinOffsets = {{0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}};
    xt::xtensor<float, 1> place = {0.0, 1000.0, 2000.0};
    auto topo = NetTopologyFixedSize(nbCells, pinCells, pinOffsets);
    xt::xtensor<float, 1> grad = topo.gradWA(place, 0.01);
    BOOST_CHECK_CLOSE(grad[0], -2.0, 0.001);
    BOOST_CHECK_SMALL(grad[1], 0.001f);
    BOOST_CHECK_CLOSE(grad[2], 2.0, 0.001);
}

BOOST_AUTO_TEST_CASE(proximalTwoPins) {
    int nbCells = 2;
    xt::xtensor<int, 2> pinCells = {{0, 1}};
    xt::xtensor<float, 2> pinOffsets = {{0.0, 0.0}};
    xt::xtensor<float, 1> place = {-1000.0, 1000.0};
    auto topo = NetTopologyFixedSize(nbCells, pinCells, pinOffsets);
    xt::xtensor<float, 1> grad = topo.proximalStep(place, 1.0);
    BOOST_CHECK_CLOSE(grad[0], 1.0, 0.001);
    BOOST_CHECK_CLOSE(grad[1], -1.0, 0.001);
}

BOOST_AUTO_TEST_CASE(proximalThreePins1) {
    int nbCells = 3;
    xt::xtensor<int, 2> pinCells = {{0, 1, 2}};
    xt::xtensor<float, 2> pinOffsets = {{0.0, 0.0, 0.0}};
    xt::xtensor<float, 1> place = {-1000.0, 1000.0, 0.0};
    auto topo = NetTopologyFixedSize(nbCells, pinCells, pinOffsets);
    xt::xtensor<float, 1> grad = topo.proximalStep(place, 1.0);
    BOOST_CHECK_CLOSE(grad[0], 1.0, 0.001);
    BOOST_CHECK_CLOSE(grad[1], -1.0, 0.001);
    BOOST_CHECK_SMALL(grad[2], 0.001f);
}

BOOST_AUTO_TEST_CASE(proximalThreePins2) {
    int nbCells = 3;
    xt::xtensor<int, 2> pinCells = {{0, 1, 2}};
    xt::xtensor<float, 2> pinOffsets = {{0.0, 0.0, 0.0}};
    xt::xtensor<float, 1> place = {-1000.0, 0.0, 0.0};
    auto topo = NetTopologyFixedSize(nbCells, pinCells, pinOffsets);
    xt::xtensor<float, 1> grad = topo.proximalStep(place, 1.0);
    BOOST_CHECK_CLOSE(grad[0], 1.0, 0.001);
    BOOST_CHECK_CLOSE(grad[1], -0.5, 0.001);
    BOOST_CHECK_CLOSE(grad[2], -0.5, 0.001);
}

BOOST_AUTO_TEST_CASE(valueHPWLOneTerminal) {
    int nbCells = 1;
    xt::xtensor<int, 2> pinCells = {{0}};
    xt::xtensor<float, 2> pinOffsets = {{0.0}};
    xt::xtensor<float, 1> minFixed = {2000.0};
    xt::xtensor<float, 1> maxFixed = {2000.0};
    xt::xtensor<float, 1> place = {1000.0};
    auto topo = NetTopologyFixedSizeTerminals(nbCells, pinCells, pinOffsets, minFixed, maxFixed);
    BOOST_CHECK_CLOSE(topo.valueHPWL(place), 1000.0, 0.001);
}

BOOST_AUTO_TEST_CASE(valueLSEOneTerminal) {
    int nbCells = 1;
    xt::xtensor<int, 2> pinCells = {{0}};
    xt::xtensor<float, 2> pinOffsets = {{0.0}};
    xt::xtensor<float, 1> minFixed = {2000.0};
    xt::xtensor<float, 1> maxFixed = {2000.0};
    xt::xtensor<float, 1> place = {1000.0};
    auto topo = NetTopologyFixedSizeTerminals(nbCells, pinCells, pinOffsets, minFixed, maxFixed);
    BOOST_CHECK_CLOSE(topo.valueLSE(place, 0.01), 1000.0, 0.001);
}

BOOST_AUTO_TEST_CASE(valueWAOneTerminal) {
    int nbCells = 1;
    xt::xtensor<int, 2> pinCells = {{0}};
    xt::xtensor<float, 2> pinOffsets = {{0.0}};
    xt::xtensor<float, 1> minFixed = {2000.0};
    xt::xtensor<float, 1> maxFixed = {2000.0};
    xt::xtensor<float, 1> place = {1000.0};
    auto topo = NetTopologyFixedSizeTerminals(nbCells, pinCells, pinOffsets, minFixed, maxFixed);
    BOOST_CHECK_CLOSE(topo.valueWA(place, 0.01), 1000.0, 0.001);
}

BOOST_AUTO_TEST_CASE(valueHPWLTwoTerminals){
    int nbCells = 1;
    xt::xtensor<int, 2> pinCells = {{0}};
    xt::xtensor<float, 2> pinOffsets = {{0.0}};
    xt::xtensor<float, 1> minFixed = {0.0};
    xt::xtensor<float, 1> maxFixed = {2000.0};
    xt::xtensor<float, 1> place = {1000.0};
    auto topo = NetTopologyFixedSizeTerminals(nbCells, pinCells, pinOffsets, minFixed, maxFixed);
    BOOST_CHECK_CLOSE(topo.valueHPWL(place), 2000.0, 0.001);
}

BOOST_AUTO_TEST_CASE(valueLSETwoTerminals) {
    int nbCells = 1;
    xt::xtensor<int, 2> pinCells = {{0}};
    xt::xtensor<float, 2> pinOffsets = {{0.0}};
    xt::xtensor<float, 1> minFixed = {0.0};
    xt::xtensor<float, 1> maxFixed = {2000.0};
    xt::xtensor<float, 1> place = {1000.0};
    auto topo = NetTopologyFixedSizeTerminals(nbCells, pinCells, pinOffsets, minFixed, maxFixed);
    BOOST_CHECK_CLOSE(topo.valueLSE(place, 0.01), 2000.0, 0.001);
}

BOOST_AUTO_TEST_CASE(valueWATwoTerminals) {
    int nbCells = 1;
    xt::xtensor<int, 2> pinCells = {{0}};
    xt::xtensor<float, 2> pinOffsets = {{0.0}};
    xt::xtensor<float, 1> minFixed = {0.0};
    xt::xtensor<float, 1> maxFixed = {2000.0};
    xt::xtensor<float, 1> place = {1000.0};
    auto topo = NetTopologyFixedSizeTerminals(nbCells, pinCells, pinOffsets, minFixed, maxFixed);
    BOOST_CHECK_CLOSE(topo.valueWA(place, 0.01), 2000.0, 0.001);
}

BOOST_AUTO_TEST_CASE(gradHPWLTwoTerminals){
    int nbCells = 1;
    xt::xtensor<int, 2> pinCells = {{0}};
    xt::xtensor<float, 2> pinOffsets = {{0.0}};
    xt::xtensor<float, 1> minFixed = {0.0};
    xt::xtensor<float, 1> maxFixed = {2000.0};
    xt::xtensor<float, 1> place = {1000.0};
    auto topo = NetTopologyFixedSizeTerminals(nbCells, pinCells, pinOffsets, minFixed, maxFixed);
    xt::xtensor<float, 1> grad = topo.gradHPWL(place);
    BOOST_CHECK_SMALL(grad[0], 0.001f);
}

BOOST_AUTO_TEST_CASE(gradLSETwoTerminals) {
    int nbCells = 1;
    xt::xtensor<int, 2> pinCells = {{0}};
    xt::xtensor<float, 2> pinOffsets = {{0.0}};
    xt::xtensor<float, 1> minFixed = {0.0};
    xt::xtensor<float, 1> maxFixed = {2000.0};
    xt::xtensor<float, 1> place = {1000.0};
    auto topo = NetTopologyFixedSizeTerminals(nbCells, pinCells, pinOffsets, minFixed, maxFixed);
    xt::xtensor<float, 1> grad = topo.gradLSE(place, 0.01);
    BOOST_CHECK_SMALL(grad[0], 0.001f);
}

BOOST_AUTO_TEST_CASE(gradWATwoTerminals) {
    int nbCells = 1;
    xt::xtensor<int, 2> pinCells = {{0}};
    xt::xtensor<float, 2> pinOffsets = {{0.0}};
    xt::xtensor<float, 1> minFixed = {0.0};
    xt::xtensor<float, 1> maxFixed = {2000.0};
    xt::xtensor<float, 1> place = {1000.0};
    auto topo = NetTopologyFixedSizeTerminals(nbCells, pinCells, pinOffsets, minFixed, maxFixed);
    xt::xtensor<float, 1> grad = topo.gradWA(place, 0.01);
    BOOST_CHECK_SMALL(grad[0], 0.001f);
}

BOOST_AUTO_TEST_CASE(proximalTerminal) {
    int nbCells = 2;
    xt::xtensor<int, 2> pinCells = {{0, 1}};
    xt::xtensor<float, 2> pinOffsets = {{0.0, 0.0}};
    xt::xtensor<float, 1> place = {-500.0, 2000.0};
    xt::xtensor<float, 1> minFixed = {-1000.0};
    xt::xtensor<float, 1> maxFixed = {1000.0};
    auto topo = NetTopologyFixedSizeTerminals(nbCells, pinCells, pinOffsets, minFixed, maxFixed);
    xt::xtensor<float, 1> grad = topo.proximalStep(place, 10000.0);
    BOOST_CHECK_SMALL(grad[0], 0.001f);
    BOOST_CHECK_CLOSE(grad[1], -1000.0, 0.001);
}

BOOST_AUTO_TEST_CASE(b2b_singlepin) {
    int nbCells = 4;
    NetTopologyBuilder bd(nbCells);
    bd.push({0}, {0.0f}, 4.0f, 4.0f);
    bd.push({1}, {5.0f}, 1.0f, 1.0f);
    bd.push({2}, {-6.0f}, 2.0f, 2.0f);
    bd.push({3}, {2.0f}, 3.0f, 3.0f);
    auto topo = bd.build();
    xt::xtensor<float, 1> place = {10.0f, 20.0f, 30.0f, 40.0f};
    auto res = topo.b2bSolve(place);
    BOOST_CHECK_CLOSE(res[0], 4.0, 0.001);
    BOOST_CHECK_CLOSE(res[1], -4.0, 0.001);
    BOOST_CHECK_CLOSE(res[2], 8.0, 0.001);
    BOOST_CHECK_CLOSE(res[3], 1.0, 0.001);
}

BOOST_AUTO_TEST_CASE(b2b_mid) {
    int nbCells = 2;
    NetTopologyBuilder bd(nbCells);
    bd.push({0, 1}, {0.0f, 0.0f}, 0.0f, 4.0f);
    auto topo = bd.build();
    xt::xtensor<float, 1> place = {2.0f, 3.0f};
    auto res = topo.b2bSolve(place);
    BOOST_CHECK_CLOSE(res[0], 2.0, 0.001);
    BOOST_CHECK_CLOSE(res[1], 3.0, 0.001);
}

