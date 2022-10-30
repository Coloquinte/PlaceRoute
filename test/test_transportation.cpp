#define BOOST_TEST_MODULE TRANSPORTATION

#include <boost/test/unit_test.hpp>

#include "place_global/transportation.hpp"

using namespace coloquinte;

BOOST_AUTO_TEST_CASE(TestAddSink) {
  std::vector<long long> capacities = {5, 6};
  std::vector<long long> demands = {1, 2, 3, 4, 5};
  std::vector<std::vector<float> > costs = {{0.0f, 1.0f, 2.0f, 3.0f, 4.0f},
                                            {5.0f, 6.0f, 7.0f, 8.0f, 9.0f}};
  TransportationProblem solver(capacities, demands, costs);
  solver.check();
  BOOST_CHECK(!solver.isBalanced());
  BOOST_CHECK_EQUAL(solver.nbSinks(), 2);
  BOOST_CHECK_EQUAL(solver.nbSources(), 5);

  solver.addDummyCapacity();
  solver.addDummyDemand();
  BOOST_CHECK_EQUAL(solver.nbSinks(), 3);
  BOOST_CHECK_EQUAL(solver.nbSources(), 5);
  BOOST_CHECK(solver.isBalanced());
  BOOST_CHECK(!solver.isFeasible());

  solver.makeFeasible();
  BOOST_CHECK(solver.isFeasible());
  solver.check();
}

BOOST_AUTO_TEST_CASE(TestAddSource) {
  std::vector<long long> capacities = {5, 6};
  std::vector<long long> demands = {1, 2, 3, 4};
  std::vector<std::vector<float> > costs = {{0.0f, 1.0f, 2.0f, 3.0f},
                                            {4.0f, 5.0f, 6.0f, 7.0f}};
  TransportationProblem solver(capacities, demands, costs);
  solver.check();
  BOOST_CHECK(!solver.isBalanced());
  BOOST_CHECK_EQUAL(solver.nbSinks(), 2);
  BOOST_CHECK_EQUAL(solver.nbSources(), 4);

  solver.addDummyCapacity();
  solver.addDummyDemand();
  BOOST_CHECK_EQUAL(solver.nbSinks(), 2);
  BOOST_CHECK_EQUAL(solver.nbSources(), 5);
  BOOST_CHECK(solver.isBalanced());
  BOOST_CHECK(!solver.isFeasible());

  solver.makeFeasible();
  BOOST_CHECK(solver.isFeasible());
  solver.check();
}

BOOST_AUTO_TEST_CASE(TestAlreadyFeasible) {
  std::vector<long long> capacities = {5, 5};
  std::vector<long long> demands = {1, 2, 3, 4};
  std::vector<std::vector<float> > costs = {{0.0f, 1.0f, 2.0f, 3.0f},
                                            {4.0f, 5.0f, 6.0f, 7.0f}};
  std::vector<std::vector<long long> > allocations = {{1, 2, 2, 0},
                                                      {0, 0, 1, 4}};
  TransportationProblem solver(capacities, demands, costs);
  solver.setAllocations(allocations);
  solver.check();
  BOOST_CHECK(solver.isBalanced());
  BOOST_CHECK(solver.isFeasible());
  BOOST_CHECK_EQUAL(solver.nbSinks(), 2);
  BOOST_CHECK_EQUAL(solver.nbSources(), 4);

  solver.addDummyCapacity();
  solver.addDummyDemand();
  BOOST_CHECK_EQUAL(solver.nbSinks(), 2);
  BOOST_CHECK_EQUAL(solver.nbSources(), 4);

  solver.makeFeasible();
  BOOST_CHECK(solver.isFeasible());
}

BOOST_AUTO_TEST_CASE(TestRun) {
  std::vector<long long> capacities = {12, 12, 12};
  std::vector<long long> demands = {1, 2, 3, 4, 5, 6, 7, 8};
  std::vector<std::vector<float> > costs = {
      {0.0f, 1.0f, 2.0f, 3.0f, 4.0f, 5.0f, 6.0f, 7.0f},
      {7.0f, 6.0f, 5.0f, 4.0f, 3.0f, 2.0f, 1.0f, 0.0f},
      {4.0f, 4.0f, 4.0f, 4.0f, 0.0f, 0.0f, 0.0f, 0.0f}};
  std::vector<std::vector<long long> > allocations = {{1, 2, 3, 4, 2, 0, 0, 0},
                                                      {0, 0, 0, 0, 3, 6, 3, 0},
                                                      {0, 0, 0, 0, 0, 0, 4, 8}};
  TransportationProblem solver(capacities, demands, costs);
  solver.setAllocations(allocations);
  solver.check();
  BOOST_CHECK(solver.isBalanced());
  BOOST_CHECK(solver.isFeasible());
  BOOST_CHECK_EQUAL(solver.nbSinks(), 3);
  BOOST_CHECK_EQUAL(solver.nbSources(), 8);

  solver.solve();
}
