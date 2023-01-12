#define BOOST_TEST_MODULE ROW_NEIGHBOURHOOD

#include <boost/test/unit_test.hpp>

#include "place_global/interior_point_transportation.hpp"

BOOST_AUTO_TEST_CASE(TestConstruction) {
  InteriorPointTransportation<double> pb =
      InteriorPointTransportation<double>::makeRandom(100, 10);
  pb.check();
  auto c = pb.costVector();
  auto b = pb.constraintVector();
  BOOST_CHECK_EQUAL(c.rows(), pb.nbVariables());
  BOOST_CHECK_EQUAL(b.rows(), pb.nbConstraints());
  pb.checkPrimalSolution(pb.initialPrimalSolution());
  pb.checkDualSolution(pb.initialDualSolution(1.0e-6));
  pb.checkDualSolution(pb.initialDualSolution(1.0e-2));
  pb.checkDualSolution(pb.initialDualSolution(1.0e2));
}

BOOST_AUTO_TEST_CASE(TestManhattan) {
  InteriorPointTransportation<double> pb =
      InteriorPointTransportation<double>::makeRandomFromDistance<2, 1>(100,
                                                                        10);
  pb.check();
}

BOOST_AUTO_TEST_CASE(TestEuclidean) {
  InteriorPointTransportation<double> pb =
      InteriorPointTransportation<double>::makeRandomFromDistance<2, 1>(100,
                                                                        10);
  pb.check();
}

BOOST_AUTO_TEST_CASE(TestApplyA) {
  InteriorPointTransportation<double> pb =
      InteriorPointTransportation<double>::makeRandom(20, 7);
  auto a = pb.makeDenseA();
  auto x = Eigen::VectorXd::LinSpaced(140, 0.1, 1.0);
  Eigen::VectorXd expected = a * x;
  Eigen::VectorXd ax = pb.applyA(x);
  BOOST_CHECK(ax.isApprox(expected));
}

BOOST_AUTO_TEST_CASE(TestApplyAt) {
  InteriorPointTransportation<double> pb =
      InteriorPointTransportation<double>::makeRandom(20, 7);
  auto a = pb.makeDenseA();
  auto x = Eigen::VectorXd::LinSpaced(26, 0.1, 1.0);
  Eigen::VectorXd expected = a.transpose() * x;
  Eigen::VectorXd atx = pb.applyAt(x);
  BOOST_CHECK(atx.isApprox(expected));
}

BOOST_AUTO_TEST_CASE(TestApplyAGAt) {
  InteriorPointTransportation<double> pb =
      InteriorPointTransportation<double>::makeRandom(20, 7);
  auto a = pb.makeDenseA();
  auto x = Eigen::VectorXd::LinSpaced(26, 0.1, 1.0);
  auto g = Eigen::VectorXd::LinSpaced(140, 0.1, 1.0);
  Eigen::VectorXd expected = a * g.asDiagonal() * a.transpose() * x;
  Eigen::VectorXd agatx = pb.applyAGAt(g, x);
  BOOST_CHECK(agatx.isApprox(expected));
}

BOOST_AUTO_TEST_CASE(TestApplyAGAtInv) {
  InteriorPointTransportation<double> pb =
      InteriorPointTransportation<double>::makeRandom(20, 7);
  auto a = pb.makeDenseA();
  auto x = Eigen::VectorXd::LinSpaced(26, 0.1, 1.0);
  auto g = Eigen::VectorXd::LinSpaced(140, 0.1, 1.0);
  Eigen::MatrixXd mat = a * g.asDiagonal() * a.transpose();
  Eigen::VectorXd expected = mat.lu().solve(x);
  Eigen::VectorXd agatInvx = pb.applyAGAtInv(g, x);
  BOOST_CHECK(agatInvx.isApprox(expected));
}

BOOST_AUTO_TEST_CASE(TestSolve) {
  Eigen::VectorXd dem = Eigen::VectorXd::LinSpaced(20, 0.1, 1.0);
  Eigen::VectorXd cap = Eigen::VectorXd::LinSpaced(7, 0.1, 1.0);
  cap = cap / cap.sum();
  cap = cap * dem.sum();
  Eigen::MatrixXd costs = Eigen::VectorXd::LinSpaced(20, 0.1, 1.0) *
                          Eigen::VectorXd::LinSpaced(7, 0.1, 1.0).transpose();

  InteriorPointTransportation<double> pb(dem, cap, costs);
  pb.check();
  Eigen::VectorXd x = pb.solve();
  pb.checkPrimalSolution(x);
}
