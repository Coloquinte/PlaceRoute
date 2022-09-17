#define BOOST_TEST_MODULE INCR_NET_MODEL

#include <boost/test/unit_test.hpp>

#include "place_detailed/incr_net_model.hpp"

using namespace coloquinte;

BOOST_AUTO_TEST_CASE(TestConstruction1) {
  IncrNetModelBuilder bd(0);
  IncrNetModel model = bd.build();
  model.check();
  BOOST_CHECK_EQUAL(model.value(), 0LL);
}

BOOST_AUTO_TEST_CASE(TestConstruction2) {
  IncrNetModelBuilder bd(0);
  bd.addNet({}, {});
  IncrNetModel model = bd.build();
  model.check();
  BOOST_CHECK_EQUAL(model.value(), 0LL);
  BOOST_CHECK_EQUAL(model.nbNets(), 0);
}

BOOST_AUTO_TEST_CASE(TestConstruction3) {
  IncrNetModelBuilder bd(3);
  bd.addNet({1, 2}, {0, 0});
  bd.addNet({1, 0}, {0, 0});
  IncrNetModel model = bd.build({0, 1, 2});
  model.check();
  BOOST_CHECK_EQUAL(model.value(), 2LL);
  BOOST_CHECK_EQUAL(model.nbNets(), 2);
}

BOOST_AUTO_TEST_CASE(TestConstruction4) {
  IncrNetModelBuilder bd(3);
  bd.addNet({1, 2, 0}, {8, 0, 2});
  bd.addNet({1, 0}, {0, 4});
  IncrNetModel model = bd.build({0, 1, 2});
  model.check();
  BOOST_CHECK_EQUAL(model.value(), 10LL);
  BOOST_CHECK_EQUAL(model.nbNets(), 2);
}

BOOST_AUTO_TEST_CASE(TestMove1) {
  IncrNetModelBuilder bd(3);
  bd.addNet({0, 1, 2}, {2, 4, 6});
  IncrNetModel model = bd.build({0, 1, 2});
  model.check();
  BOOST_CHECK_EQUAL(model.value(), 6LL);
  BOOST_CHECK_EQUAL(model.nbNets(), 1);
  model.updateCellPos(0, 2);
  model.check();
  model.updateCellPos(0, 2);
  model.check();
  model.updateCellPos(1, -2);
  model.check();
  model.updateCellPos(1, 8);
  model.check();
  model.updateCellPos(2, 1);
  model.check();
}

BOOST_AUTO_TEST_CASE(TestMove2) {
  IncrNetModelBuilder bd(3);
  bd.addNet({0, 1, 2}, {0, 0, 0});
  IncrNetModel model = bd.build({0, 1, 2});
  model.check();
  BOOST_CHECK_EQUAL(model.nbNets(), 1);
  BOOST_CHECK_EQUAL(model.value(), 2LL);
  model.updateCellPos(0, 2);
  model.check();
  BOOST_CHECK_EQUAL(model.value(), 1LL);
  model.updateCellPos(0, 2);
  model.check();
  BOOST_CHECK_EQUAL(model.value(), 1LL);
  model.updateCellPos(1, -2);
  model.check();
  BOOST_CHECK_EQUAL(model.value(), 4LL);
  model.updateCellPos(1, 8);
  model.check();
  BOOST_CHECK_EQUAL(model.value(), 6LL);
  model.updateCellPos(2, 1);
  model.check();
  BOOST_CHECK_EQUAL(model.value(), 7LL);
}
