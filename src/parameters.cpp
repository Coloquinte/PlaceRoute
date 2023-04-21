

#include <cmath>
#include <sstream>

#include "coloquinte.hpp"

namespace coloquinte {
std::string toString(LegalizationModel model) {
  switch (model) {
    case LegalizationModel::L1:
      return "L1";
    case LegalizationModel::L2:
      return "L2";
    case LegalizationModel::LInf:
      return "LInf";
    case LegalizationModel::L1Squared:
      return "L1Squared";
    case LegalizationModel::L2Squared:
      return "L2Squared";
    case LegalizationModel::LInfSquared:
      return "LInfSquared";
    default:
      return "UnknownLegalizationModel";
  }
}

std::string toString(CellOrientation o) {
  switch (o) {
    case CellOrientation::N:
      return "N";
    case CellOrientation::S:
      return "S";
    case CellOrientation::E:
      return "E";
    case CellOrientation::W:
      return "W";
    case CellOrientation::FN:
      return "FN";
    case CellOrientation::FS:
      return "FS";
    case CellOrientation::FE:
      return "FE";
    case CellOrientation::FW:
      return "FW";
    default:
      return "UnknownCellOrientation";
  }
}

CellOrientation oppositeRowOrientation(CellOrientation o) {
  switch (o) {
    case CellOrientation::N:
      return CellOrientation::FS;
    case CellOrientation::S:
      return CellOrientation::FN;
    case CellOrientation::E:
      return CellOrientation::FW;
    case CellOrientation::W:
      return CellOrientation::FE;
    case CellOrientation::FN:
      return CellOrientation::S;
    case CellOrientation::FS:
      return CellOrientation::N;
    case CellOrientation::FE:
      return CellOrientation::W;
    case CellOrientation::FW:
      return CellOrientation::E;
    default:
      return CellOrientation::INVALID;
  }
}

bool isTurn(CellOrientation orient) {
  return orient == CellOrientation::E || orient == CellOrientation::W ||
         orient == CellOrientation::FW || orient == CellOrientation::FE;
}

std::string toString(NetModelOption model) {
  switch (model) {
    case NetModelOption::BoundToBound:
      return "BoundToBound";
    case NetModelOption::Star:
      return "Star";
    default:
      return "UnknownNetModel";
  }
}

std::string toString(CellRowPolarity pol) {
  switch (pol) {
    case CellRowPolarity::SAME:
      return "SAME";
    case CellRowPolarity::OPPOSITE:
      return "OPPOSITE";
    case CellRowPolarity::ANY:
      return "ANY";
    default:
      return "UnknownCellRowPolarity";
  }
}

std::string Rectangle::toString() const {
  std::stringstream ss;
  ss << "Rectangle " << minX << ".." << maxX << " x " << minY << ".." << maxY;
  return ss.str();
}

namespace {

double interpolateEffort(double minVal, double maxVal, int effort,
                         int minEffort = 1, int maxEffort = 9) {
  assert(minEffort < maxEffort);
  assert(effort >= minEffort && effort <= maxEffort);
  double fact = (effort - minEffort) / (float)(maxEffort - minEffort);
  return maxVal * fact + minVal * (1.0 - fact);
}

double interpolateLogEffort(double minVal, double maxVal, int effort,
                            int minEffort = 1, int maxEffort = 9) {
  return std::exp(interpolateEffort(std::log(minVal), std::log(maxVal), effort,
                                    minEffort, maxEffort));
}
}  // namespace

ColoquinteParameters::ColoquinteParameters(int effort, int seed)
    : seed(seed), global(effort), legalization(effort), detailed(effort) {
  if (effort < 1 || effort > 9) {
    throw std::runtime_error("Placement effort must be between 1 and 9");
  }
}

RoughLegalizationParameters::RoughLegalizationParameters(int effort) {
  costModel = LegalizationModel::L1;
  nbSteps = 1;
  // TODO: find best parameter
  binSize = 5.0;
  lineReoptSize = 2;
  lineReoptOverlap = 1;
  diagReoptSize = 2;
  diagReoptOverlap = 1;
  unidimensionalTransport = false;
  // TODO: find best parameter
  sideMargin = 0.9;
  coarseningLimit = 100.0;
  quadraticPenalty = 0.001;
  // TODO: find best parameter
  targetBlending = 0.0;
  int squareSizeArray[9] = {1, 2, 3, 3, 3, 4, 4, 4, 5};
  squareReoptSize = squareSizeArray[effort - 1];
  squareReoptOverlap = 1;
}

PenaltyParameters::PenaltyParameters(int effort) {
  // TODO: make cutoff distance smaller at small effort
  cutoffDistance = 40.0;
  cutoffDistanceUpdateFactor = 1.0;
  areaExponent = 0.5;
  // TODO: make initial penalty bigger at small effort
  initialValue = 0.03;
  // TODO: find best parameter
  targetBlending = 1.0;
  double updateFactorArray[9] = {1.23, 1.23, 1.23, 1.22, 1.22,
                                 1.22, 1.22, 1.17, 1.07};
  updateFactor = updateFactorArray[effort - 1];
}

ContinuousModelParameters::ContinuousModelParameters(int effort) {
  netModel = NetModelOption::BoundToBound;
  approximationDistance = 2.0;
  approximationDistanceUpdateFactor = 1.0;
  maxNbConjugateGradientSteps = 1000;
  conjugateGradientErrorTolerance = 1.0e-6;
}

GlobalPlacerParameters::GlobalPlacerParameters(int effort)
    : continuousModel(effort), roughLegalization(effort), penalty(effort) {
  maxNbSteps = 400;
  nbInitialSteps = 0;
  nbStepsBeforeRoughLegalization = 1;
  distanceTolerance = 2.0;
  // TODO: find best parameter
  exportBlending = 0.99;
  noise = 1.0e-4;
  // Parameters that vary with effort here
  double gapToleranceArray[9] = {0.13,  0.13,  0.058, 0.038, 0.026,
                                 0.026, 0.026, 0.026, 0.026};
  gapTolerance = gapToleranceArray[effort - 1];
  check();
}

std::string ColoquinteParameters::toString() const {
  return global.toString() + legalization.toString() + detailed.toString();
}

std::string GlobalPlacerParameters::toString() const {
  std::stringstream ss;
  ss << "Global placer params:"
     << "\n\tGap tolerance: " << gapTolerance
     << "\n\tDistance tolerance: " << distanceTolerance
     << "\n\tMax nb steps: " << maxNbSteps
     << "\n\tInitial placement steps: " << nbInitialSteps
     << "\n\tPlacement steps per legalization: "
     << nbStepsBeforeRoughLegalization
     << "\n\tExport blending: " << exportBlending;
  ss << std::endl;
  return ss.str();
}

DetailedPlacerParameters::DetailedPlacerParameters(int effort) {
  nbPasses = std::round(interpolateLogEffort(2.0, 8.0, effort));
  localSearchNbNeighbours = std::round(interpolateLogEffort(2.0, 16.0, effort));
  localSearchNbRows = std::round(interpolateEffort(1.0, 4.0, effort));
  shiftNbRows = 3;
  shiftMaxNbCells = std::round(interpolateLogEffort(50, 120.0, effort));
  reorderingNbRows = 1;
  reorderingMaxNbCells = 1;
  check();
}

LegalizationParameters::LegalizationParameters(int effort) {
  costModel = LegalizationModel::L1;
  orderingHeight = -1.0;
  // TODO: find best parameter
  orderingWidth = 0.2;
  orderingY = 0.0;
  check();
}

void PenaltyParameters::check() const {
  if (cutoffDistance < 1.0e-6) {
    throw std::runtime_error("Too small cutoff distance may lead to issues");
  }
  if (cutoffDistanceUpdateFactor < 0.8 || cutoffDistanceUpdateFactor > 1.2) {
    throw std::runtime_error(
        "Penalty cutoff update factor should be close to 1");
  }
  if (areaExponent < 0.49 || areaExponent > 1.01) {
    throw std::runtime_error(
        "Penalty area exponent should be between 0.5 and 1");
  }
  if (initialValue <= 0.0f) {
    throw std::runtime_error("Initial penalty should be positive");
  }
  if (updateFactor <= 1.0f || updateFactor >= 2.0f) {
    throw std::runtime_error(
        "Penalty update factor should be between one and two");
  }
  if (targetBlending < 0.1f || targetBlending > 1.1f) {
    throw std::runtime_error(
        "Penalty target blending should generally be between 0.5 and 1");
  }
}

void ContinuousModelParameters::check() const {
  if (approximationDistance < 1.0e-6) {
    throw std::runtime_error(
        "Too small approximation distance may lead to issues");
  }
  if (approximationDistanceUpdateFactor < 0.8 ||
      approximationDistanceUpdateFactor > 1.2) {
    throw std::runtime_error(
        "Approximation distance update factor should be close to 1");
  }
  if (approximationDistance > 1.0e3) {
    throw std::runtime_error(
        "Too large approximation distance is highly imprecise");
  }
  if (maxNbConjugateGradientSteps <= 0) {
    throw std::runtime_error(
        "Must have positive number of steps during conjugate gradients");
  }
  if (conjugateGradientErrorTolerance < 1.0e-8) {
    throw std::runtime_error("Too small error tolerance may lead to issues");
  }
  if (conjugateGradientErrorTolerance > 1.0) {
    throw std::runtime_error("Too large error tolerance is highly imprecise");
  }
}

void RoughLegalizationParameters::check() const {
  if (nbSteps < 0) {
    throw std::runtime_error(
        "Must have non-negative number of steps for rough legalization");
  }
  if (binSize < 1.0f) {
    throw std::runtime_error(
        "Bin size should generally be larger than 1 (one standard cell)");
  }
  if (binSize > 25.0f) {
    throw std::runtime_error(
        "Bin size should not be too large (10 should be enough)");
  }
  if (lineReoptSize < 1 || diagReoptSize < 1 || squareReoptSize < 1) {
    throw std::runtime_error(
        "Rough legalization reopt size should be at least 1");
  }
  if (lineReoptOverlap < 1 || diagReoptOverlap < 1 || squareReoptOverlap < 1) {
    throw std::runtime_error(
        "Rough legalization reopt overlap should be at least 1");
  }
  if (lineReoptSize > 64 || diagReoptSize > 64 || squareReoptSize > 8) {
    throw std::runtime_error("Rough legalization reopt should be small");
  }
  if (lineReoptSize < 2 && diagReoptSize < 2 && squareReoptSize < 2 &&
      !unidimensionalTransport) {
    throw std::runtime_error(
        "At least one rough legalization reopt value should be 2 or more");
  }
  if (lineReoptSize > 1 && lineReoptOverlap >= lineReoptSize) {
    throw std::runtime_error(
        "Rough legalization reopt overlap should be smaller than reopt size");
  }
  if (diagReoptSize > 1 && diagReoptOverlap >= diagReoptSize) {
    throw std::runtime_error(
        "Rough legalization reopt overlap should be smaller than reopt size");
  }
  if (squareReoptSize > 1 && squareReoptOverlap >= squareReoptSize) {
    throw std::runtime_error(
        "Rough legalization reopt overlap should be smaller than reopt size");
  }
  if (quadraticPenalty < 0.0 || quadraticPenalty > 1.0) {
    throw std::runtime_error(
        "Rough legalization quadratic penalty should be non-negative and small "
        "(< 1.0)");
  }
  if (targetBlending < -0.1 || targetBlending > 0.9f) {
    throw std::runtime_error(
        "Rough legalization target blending should generally be between 0 and "
        "0.5");
  }

  if (costModel != LegalizationModel::L1 && unidimensionalTransport) {
    throw std::runtime_error(
        "Unidimensional transport can only be used with L1 cost model");
  }
}

void GlobalPlacerParameters::check() const {
  roughLegalization.check();
  continuousModel.check();
  penalty.check();
  if (maxNbSteps < 0) {
    throw std::runtime_error("Invalid number of steps");
  }
  if (nbInitialSteps < 0) {
    throw std::runtime_error("Invalid number of initial steps");
  }
  if (nbInitialSteps >= maxNbSteps) {
    throw std::runtime_error(
        "Number of initial steps should be lower than max number");
  }
  if (nbStepsBeforeRoughLegalization < 1) {
    throw std::runtime_error(
        "Number of steps per legalization should be positive");
  }
  if (gapTolerance < 0.0f || gapTolerance > 1.0f) {
    throw std::runtime_error(
        "Invalid gap tolerance (should be between 0 and 1)");
  }
  if (distanceTolerance < 0.0f) {
    throw std::runtime_error(
        "Invalid distance tolerance (should be non-negative");
  }
  if (exportBlending < -0.5f || exportBlending > 1.5f) {
    throw std::runtime_error(
        "Export blending should generally be between 0 and 1");
  }
  if (noise < 0.0 || noise > 2.0) {
    throw std::runtime_error(
        "Noise should be a very small non-negative number");
  }
}

void LegalizationParameters::check() const {
  if (costModel != LegalizationModel::L1) {
    throw std::runtime_error("Only L1 legalization model is supported");
  }
  if (orderingWidth > 2.0 || orderingWidth < -1.0) {
    throw std::runtime_error(
        "Legalization ordering width should be small (0 < ... < 1)");
  }
  if (orderingY > 0.2 || orderingY < -0.2) {
    throw std::runtime_error(
        "Legalization ordering y should be small (-0.1 < ... < 0.1)");
  }
}

void DetailedPlacerParameters::check() const {
  if (nbPasses < 0) {
    throw std::runtime_error(
        "Number of detailed placement passes must be non-negative");
  }
  if (localSearchNbNeighbours < 0) {
    throw std::runtime_error(
        "Number of detailed placement neighbour cells must be non-negative");
  }
  if (localSearchNbRows < 0) {
    throw std::runtime_error(
        "Number of detailed placement rows must be non-negative");
  }
  if (shiftNbRows <= 0) {
    throw std::runtime_error(
        "Number of detailed placement shift rows must be positive");
  }
  if (shiftMaxNbCells < 0) {
    throw std::runtime_error(
        "Number of detailed placement shift cells must be non-negative");
  }
  if (reorderingNbRows <= 0) {
    throw std::runtime_error(
        "Number of detailed placement reordering rows must be positive");
  }
  if (reorderingMaxNbCells < 0) {
    throw std::runtime_error(
        "Number of detailed placement reordering cells must be non-negative");
  }
}

void ColoquinteParameters::check() const {
  global.check();
  legalization.check();
  detailed.check();
}

std::string LegalizationParameters::toString() const {
  std::stringstream ss;
  ss << "legalization:"
     << "\n\tLegalization cost model: " << coloquinte::toString(costModel)
     << "\n\tordering width: " << orderingWidth
     << "\n\t ordering y: " << orderingY;
  ss << std::endl;
  return ss.str();
}

std::string DetailedPlacerParameters::toString() const {
  std::stringstream ss;
  ss << "detailed:"
     << "\n\tNb passes: " << nbPasses
     << "\n\tLocal search nb neighbours: " << localSearchNbNeighbours
     << "\n\tLocal search nb rows: " << localSearchNbRows
     << "\n\tReordering max nb rows: " << reorderingNbRows
     << "\n\tReordering max nb cells: " << reorderingMaxNbCells
     << "\n\tShift nb rows: " << shiftNbRows
     << "\n\tShift max nb cells: " << shiftMaxNbCells;
  ss << std::endl;
  return ss.str();
}

std::string ContinuousModelParameters::toString() const {
  std::stringstream ss;
  ss << "continuous model:"
     << "\n\tNet model: " << coloquinte::toString(netModel)
     << "\n\tApproximation distance: " << approximationDistance
     << "\n\tMax nb CG steps: " << maxNbConjugateGradientSteps
     << "\n\tCG error tolerance: " << conjugateGradientErrorTolerance
     << std::endl;
  return ss.str();
}

std::string RoughLegalizationParameters::toString() const {
  std::stringstream ss;
  ss << "rough legalization:"
     << "\n\tCost model: " << coloquinte::toString(costModel)
     << "\n\tNb steps: " << nbSteps << "\n\tBin size: " << binSize
     << "\n\tLine reopt length: " << lineReoptSize
     << "\n\tLine reopt overlap: " << lineReoptOverlap
     << "\n\tDiag reopt size: " << diagReoptSize
     << "\n\tDiag reopt overlap: " << diagReoptOverlap
     << "\n\tSquare reopt size: " << squareReoptSize
     << "\n\tSquare reopt overlap: " << squareReoptOverlap
     << "\n\t1D transport: " << unidimensionalTransport
     << "\n\tCoarsening limit: " << coarseningLimit
     << "\n\tQuadratic penalty: " << quadraticPenalty
     << "\n\tSide margin: " << sideMargin
     << "\n\tTarget blending: " << targetBlending << std::endl;
  return ss.str();
}

std::string PenaltyParameters::toString() const {
  std::stringstream ss;
  ss << "penalty:"
     << "\n\tCutoff distance: " << cutoffDistance
     << "\n\tArea exponent: " << areaExponent
     << "\n\tInitial value: " << initialValue
     << "\n\tUpdate factor: " << updateFactor
     << "\n\tTarget blending: " << targetBlending << std::endl;
  return ss.str();
}
}  // namespace coloquinte
