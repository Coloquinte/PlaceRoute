#define BOOST_TEST_MODULE ROW_LEGALIZER

#include <boost/test/unit_test.hpp>

#include "place_detailed/legalizer.hpp"

using namespace coloquinte;

BOOST_AUTO_TEST_CASE(TestStats) {
  std::vector<int> widths = {};
  std::vector<int> heights = {};
  std::vector<int> cellX = {};
  std::vector<int> cellY = {};
  std::vector<CellRowPolarity> polarities = {};
  std::vector<CellOrientation> orientations = {};
  std::vector<Row> rows = {Row(0, 15, 10, 20, CellOrientation::N),
                           Row(0, 15, 20, 30, CellOrientation::N)};

  ColoquinteParameters params;
  Legalizer leg(rows, widths, heights, polarities, cellX, cellY, orientations);
  leg.check();
  BOOST_CHECK_EQUAL(leg.rowHeight(), 10);
}

BOOST_AUTO_TEST_CASE(TestBasic1) {
  std::vector<int> widths = {5};
  std::vector<int> heights = {10};
  std::vector<int> cellX = {0};
  std::vector<int> cellY = {0};
  std::vector<CellRowPolarity> polarities = {CellRowPolarity::SAME};
  std::vector<CellOrientation> orientations = {CellOrientation::N};
  std::vector<Row> rows = {Row(0, 15, 10, 20, CellOrientation::N),
                           Row(0, 15, 20, 30, CellOrientation::N)};

  ColoquinteParameters params;
  Legalizer leg(rows, widths, heights, polarities, cellX, cellY, orientations);
  leg.run(params);
  BOOST_CHECK_EQUAL(leg.cellLegalX()[0], 0);
  BOOST_CHECK_EQUAL(leg.cellLegalY()[0], 10);
  BOOST_CHECK_EQUAL(leg.cellLegalOrientation()[0], CellOrientation::N);
}

BOOST_AUTO_TEST_CASE(TestBasic2) {
  std::vector<int> widths = {5};
  std::vector<int> heights = {10};
  std::vector<int> cellX = {20};
  std::vector<int> cellY = {16};
  std::vector<CellRowPolarity> polarities = {CellRowPolarity::SAME};
  std::vector<CellOrientation> orientations = {CellOrientation::N};
  std::vector<Row> rows = {Row(0, 15, 10, 20, CellOrientation::N),
                           Row(0, 15, 20, 30, CellOrientation::N)};

  ColoquinteParameters params;
  Legalizer leg(rows, widths, heights, polarities, cellX, cellY, orientations);
  leg.run(params);
  BOOST_CHECK_EQUAL(leg.cellLegalX()[0], 10);
  BOOST_CHECK_EQUAL(leg.cellLegalY()[0], 20);
  BOOST_CHECK_EQUAL(leg.cellLegalOrientation()[0], CellOrientation::N);
}

BOOST_AUTO_TEST_CASE(TestSameOrientation1) {
  std::vector<int> widths = {5, 5};
  std::vector<int> heights = {10, 10};
  std::vector<int> cellX = {0, 0};
  std::vector<int> cellY = {0, 20};
  std::vector<CellRowPolarity> polarities = {CellRowPolarity::SAME, CellRowPolarity::SAME};
  std::vector<CellOrientation> orientations = {CellOrientation::N, CellOrientation::N};
  std::vector<Row> rows = {Row(0, 15, 10, 20, CellOrientation::S),
                           Row(0, 15, 20, 30, CellOrientation::N)};

  ColoquinteParameters params;
  Legalizer leg(rows, widths, heights, polarities, cellX, cellY, orientations);
  leg.run(params);
  BOOST_CHECK_EQUAL(leg.cellLegalX()[0], 0);
  BOOST_CHECK_EQUAL(leg.cellLegalY()[0], 10);
  BOOST_CHECK_EQUAL(leg.cellLegalOrientation()[0], CellOrientation::S);
  BOOST_CHECK_EQUAL(leg.cellLegalX()[1], 0);
  BOOST_CHECK_EQUAL(leg.cellLegalY()[1], 20);
  BOOST_CHECK_EQUAL(leg.cellLegalOrientation()[1], CellOrientation::N);
}

BOOST_AUTO_TEST_CASE(TestOppositeOrientation1) {
  std::vector<int> widths = {5, 5};
  std::vector<int> heights = {10, 10};
  std::vector<int> cellX = {0, 0};
  std::vector<int> cellY = {0, 20};
  std::vector<CellRowPolarity> polarities = {CellRowPolarity::OPPOSITE, CellRowPolarity::OPPOSITE};
  std::vector<CellOrientation> orientations = {CellOrientation::N, CellOrientation::N};
  std::vector<Row> rows = {Row(0, 15, 10, 20, CellOrientation::S),
                           Row(0, 15, 20, 30, CellOrientation::N)};

  ColoquinteParameters params;
  Legalizer leg(rows, widths, heights, polarities, cellX, cellY, orientations);
  leg.run(params);
  BOOST_CHECK_EQUAL(leg.cellLegalX()[0], 0);
  BOOST_CHECK_EQUAL(leg.cellLegalY()[0], 10);
  BOOST_CHECK_EQUAL(leg.cellLegalOrientation()[0], CellOrientation::FN);
  BOOST_CHECK_EQUAL(leg.cellLegalX()[1], 0);
  BOOST_CHECK_EQUAL(leg.cellLegalY()[1], 20);
  BOOST_CHECK_EQUAL(leg.cellLegalOrientation()[1], CellOrientation::FS);
}

BOOST_AUTO_TEST_CASE(TestTwoRowCell1) {
  std::vector<int> widths = {5};
  std::vector<int> heights = {20};
  std::vector<int> cellX = {0};
  std::vector<int> cellY = {0};
  std::vector<int> cellIndex = {0, 1, 2, 3, 4};
  std::vector<CellRowPolarity> polarities = {CellRowPolarity::SAME};
  std::vector<CellOrientation> orientations = {CellOrientation::N};
  std::vector<Row> rows = {Row(0, 15, 10, 20, CellOrientation::N),
                           Row(0, 15, 20, 30, CellOrientation::S)};

  ColoquinteParameters params;
  Legalizer leg(rows, widths, heights, polarities, cellX, cellY, orientations);
  leg.run(params);
  BOOST_CHECK_EQUAL(leg.cellLegalX()[0], 0);
  BOOST_CHECK_EQUAL(leg.cellLegalY()[0], 10);
  BOOST_CHECK_EQUAL(leg.cellLegalOrientation()[0], CellOrientation::N);
}

BOOST_AUTO_TEST_CASE(TestTwoRowCell2) {
  std::vector<int> widths = {5};
  std::vector<int> heights = {20};
  std::vector<int> cellX = {0};
  std::vector<int> cellY = {20};
  std::vector<int> cellIndex = {0, 1, 2, 3, 4};
  std::vector<CellRowPolarity> polarities = {CellRowPolarity::SAME};
  std::vector<CellOrientation> orientations = {CellOrientation::N};
  std::vector<Row> rows = {Row(0, 15, 10, 20, CellOrientation::N),
                           Row(0, 15, 20, 30, CellOrientation::S)};

  ColoquinteParameters params;
  Legalizer leg(rows, widths, heights, polarities, cellX, cellY, orientations);
  leg.run(params);
  BOOST_CHECK_EQUAL(leg.cellLegalX()[0], 0);
  BOOST_CHECK_EQUAL(leg.cellLegalY()[0], 10);
  BOOST_CHECK_EQUAL(leg.cellLegalOrientation()[0], CellOrientation::N);
}

BOOST_AUTO_TEST_CASE(TestTwoRowCell3) {
  std::vector<int> widths = {5};
  std::vector<int> heights = {20};
  std::vector<int> cellX = {0};
  std::vector<int> cellY = {20};
  std::vector<CellRowPolarity> polarities = {CellRowPolarity::SAME};
  std::vector<CellOrientation> orientations = {CellOrientation::N};
  std::vector<Row> rows = {Row(0, 15, 10, 20, CellOrientation::N),
                           Row(0, 15, 20, 30, CellOrientation::S),
                           Row(0, 15, 30, 40, CellOrientation::N)};

  ColoquinteParameters params;
  Legalizer leg(rows, widths, heights, polarities, cellX, cellY, orientations);
  leg.run(params);
  BOOST_CHECK_EQUAL(leg.cellLegalX()[0], 0);
  BOOST_CHECK_EQUAL(leg.cellLegalY()[0], 20);
  BOOST_CHECK_EQUAL(leg.cellLegalOrientation()[0], CellOrientation::S);
}

BOOST_AUTO_TEST_CASE(TestTwoRowCell4) {
  std::vector<int> widths = {5};
  std::vector<int> heights = {20};
  std::vector<int> cellX = {0};
  std::vector<int> cellY = {20};
  std::vector<int> cellIndex = {0, 1, 2, 3, 4};
  std::vector<CellRowPolarity> polarities = {CellRowPolarity::OPPOSITE};
  std::vector<CellOrientation> orientations = {CellOrientation::N};
  std::vector<Row> rows = {Row(0, 15, 10, 20, CellOrientation::N),
                           Row(0, 15, 20, 30, CellOrientation::S),
                           Row(0, 15, 30, 40, CellOrientation::N)};

  ColoquinteParameters params;
  Legalizer leg(rows, widths, heights, polarities, cellX, cellY, orientations);
  leg.run(params);
  BOOST_CHECK_EQUAL(leg.cellLegalX()[0], 0);
  BOOST_CHECK_EQUAL(leg.cellLegalY()[0], 20);
  BOOST_CHECK_EQUAL(leg.cellLegalOrientation()[0], CellOrientation::FN);
}

BOOST_AUTO_TEST_CASE(TestTwoRowCell5) {
  std::vector<int> widths = {5};
  std::vector<int> heights = {20};
  std::vector<int> cellX = {0};
  std::vector<int> cellY = {20};
  std::vector<CellRowPolarity> polarities = {CellRowPolarity::NW};
  std::vector<CellOrientation> orientations = {CellOrientation::N};
  std::vector<Row> rows = {Row(0, 15, 10, 20, CellOrientation::N),
                           Row(0, 15, 20, 30, CellOrientation::S),
                           Row(0, 15, 30, 40, CellOrientation::N)};

  ColoquinteParameters params;
  Legalizer leg(rows, widths, heights, polarities, cellX, cellY, orientations);
  leg.run(params);
  BOOST_CHECK_EQUAL(leg.cellLegalX()[0], 0);
  BOOST_CHECK_EQUAL(leg.cellLegalY()[0], 10);
  BOOST_CHECK_EQUAL(leg.cellLegalOrientation()[0], CellOrientation::N);
}

BOOST_AUTO_TEST_CASE(TestTwoRowCell6) {
  std::vector<int> widths = {5};
  std::vector<int> heights = {20};
  std::vector<int> cellX = {0};
  std::vector<int> cellY = {20};
  std::vector<CellRowPolarity> polarities = {CellRowPolarity::SE};
  std::vector<CellOrientation> orientations = {CellOrientation::N};
  std::vector<Row> rows = {Row(0, 15, 10, 20, CellOrientation::N),
                           Row(0, 15, 20, 30, CellOrientation::S),
                           Row(0, 15, 30, 40, CellOrientation::N)};

  ColoquinteParameters params;
  Legalizer leg(rows, widths, heights, polarities, cellX, cellY, orientations);
  leg.run(params);
  BOOST_CHECK_EQUAL(leg.cellLegalX()[0], 0);
  BOOST_CHECK_EQUAL(leg.cellLegalY()[0], 20);
  BOOST_CHECK_EQUAL(leg.cellLegalOrientation()[0], CellOrientation::S);
}

BOOST_AUTO_TEST_CASE(TestThreeRowCell1) {
  std::vector<int> widths = {5};
  std::vector<int> heights = {30};
  std::vector<int> cellX = {8};
  std::vector<int> cellY = {30};
  std::vector<int> cellIndex = {0, 1, 2, 3, 4};
  std::vector<CellRowPolarity> polarities = {CellRowPolarity::SAME};
  std::vector<CellOrientation> orientations = {CellOrientation::N};
  std::vector<Row> rows = {Row(0, 15, 10, 20, CellOrientation::N),
                           Row(0, 15, 20, 30, CellOrientation::S),
                           Row(0, 15, 30, 40, CellOrientation::N),
                           Row(0, 15, 40, 50, CellOrientation::S)};

  ColoquinteParameters params;
  Legalizer leg(rows, widths, heights, polarities, cellX, cellY, orientations);
  leg.run(params);
  BOOST_CHECK_EQUAL(leg.cellLegalX()[0], 8);
  BOOST_CHECK_EQUAL(leg.cellLegalY()[0], 20);
  BOOST_CHECK_EQUAL(leg.cellLegalOrientation()[0], CellOrientation::S);
}

BOOST_AUTO_TEST_CASE(TestThreeRowCell2) {
  std::vector<int> widths = {5};
  std::vector<int> heights = {30};
  std::vector<int> cellX = {8};
  std::vector<int> cellY = {0};
  std::vector<int> cellIndex = {0, 1, 2, 3, 4};
  std::vector<CellRowPolarity> polarities = {CellRowPolarity::OPPOSITE};
  std::vector<CellOrientation> orientations = {CellOrientation::N};
  std::vector<Row> rows = {Row(0, 15, 10, 20, CellOrientation::N),
                           Row(0, 15, 20, 30, CellOrientation::S),
                           Row(0, 15, 30, 40, CellOrientation::N),
                           Row(0, 15, 40, 50, CellOrientation::S)};

  ColoquinteParameters params;
  Legalizer leg(rows, widths, heights, polarities, cellX, cellY, orientations);
  leg.run(params);
  BOOST_CHECK_EQUAL(leg.cellLegalX()[0], 8);
  BOOST_CHECK_EQUAL(leg.cellLegalY()[0], 10);
  BOOST_CHECK_EQUAL(leg.cellLegalOrientation()[0], CellOrientation::FS);
}