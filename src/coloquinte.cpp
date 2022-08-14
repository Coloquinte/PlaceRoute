
#include "coloquinte.hpp"

#include <boost/polygon/polygon.hpp>
#include <iostream>

#include "coloquinte.hpp"
#include "place_detailed/place_detailed.hpp"
#include "place_global/density_legalizer.hpp"
#include "place_global/net_model.hpp"
#include "place_global/place_global.hpp"

namespace bpl = boost::polygon;
namespace coloquinte {

Circuit::Circuit(int nbCells) {
  cellWidths.resize(nbCells);
  cellHeights.resize(nbCells);
  cellFixed.resize(nbCells);
  cellX.resize(nbCells);
  cellY.resize(nbCells);
  cellOrientation.resize(nbCells);
  netLimits.push_back(0);
  check();
}

void Circuit::setNets(const std::vector<int> &limits,
                      const std::vector<int> &cells,
                      const std::vector<int> &xOffsets,
                      const std::vector<int> &yOffsets) {
  assert(!limits.empty());
  assert(limits.front() == 0);
  assert(limits.back() == cells.size());
  assert(limits.back() == xOffsets.size());
  assert(limits.back() == yOffsets.size());
  netLimits = limits;
  pinCells = cells;
  pinXOffsets = xOffsets;
  pinYOffsets = yOffsets;
}

int Circuit::width(int cell) const {
  CellOrientation orient = orientation(cell);
  bool turned = orient == CellOrientation::E || orient == CellOrientation::W ||
                orient == CellOrientation::FW || orient == CellOrientation::FE;
  return turned ? cellHeights[cell] : cellWidths[cell];
}

int Circuit::height(int cell) const {
  CellOrientation orient = orientation(cell);
  bool turned = orient == CellOrientation::E || orient == CellOrientation::W ||
                orient == CellOrientation::FW || orient == CellOrientation::FE;
  return turned ? cellWidths[cell] : cellHeights[cell];
}

int Circuit::pinXOffset(int net, int i) const {
  int cell = pinCell(net, i);
  CellOrientation orient = orientation(cell);
  bool turned = orient == CellOrientation::E || orient == CellOrientation::W ||
                orient == CellOrientation::FW || orient == CellOrientation::FE;
  int offs = turned ? pinYOffsets[netLimits[net] + i]
                    : pinXOffsets[netLimits[net] + i];
  bool flipped = orient == CellOrientation::S || orient == CellOrientation::W ||
                 orient == CellOrientation::FN || orient == CellOrientation::FE;
  return flipped ? width(cell) - offs : offs;
}

int Circuit::pinYOffset(int net, int i) const {
  int cell = pinCell(net, i);
  CellOrientation orient = orientation(cell);
  bool turned = orient == CellOrientation::E || orient == CellOrientation::W ||
                orient == CellOrientation::FW || orient == CellOrientation::FE;
  int offs = turned ? pinXOffsets[netLimits[net] + i]
                    : pinYOffsets[netLimits[net] + i];
  bool flipped = orient == CellOrientation::S || orient == CellOrientation::E ||
                 orient == CellOrientation::FS || orient == CellOrientation::FE;
  return flipped ? height(cell) - offs : offs;
}

long long Circuit::hpwl() const {
  long long ret = 0;
  for (int net = 0; net < nbNets(); ++net) {
    if (nbPins(net) == 0) {
      continue;
    }
    int minX = std::numeric_limits<int>::max();
    int maxX = std::numeric_limits<int>::min();
    int minY = std::numeric_limits<int>::max();
    int maxY = std::numeric_limits<int>::min();
    for (int pin = 0; pin < nbPins(net); ++pin) {
      int cell = pinCell(net, pin);
      int px = x(cell) + pinXOffset(net, pin);
      int py = y(cell) + pinYOffset(net, pin);
      minX = std::min(px, minX);
      maxX = std::max(px, maxX);
      minY = std::min(py, minY);
      maxY = std::max(py, maxY);
    }
    ret += (maxX - minX);
    ret += (maxY - minY);
  }
  return ret;
}

Rectangle Circuit::computePlacementArea() const {
  int minX = std::numeric_limits<int>::max();
  int maxX = std::numeric_limits<int>::min();
  int minY = std::numeric_limits<int>::max();
  int maxY = std::numeric_limits<int>::min();
  if (rows.empty()) {
    return Rectangle(0, 0, 0, 0);
  }
  for (Rectangle row : rows) {
    minX = std::min(row.minX, minX);
    maxX = std::max(row.maxX, maxX);
    minY = std::min(row.minY, minY);
    maxY = std::max(row.maxY, maxY);
  }
  return Rectangle(minX, maxX, minY, maxY);
}

std::vector<Rectangle> Circuit::computeRows() const {
  std::vector<Rectangle> obstacles;
  for (int i = 0; i < nbCells(); ++i) {
    if (!fixed(i)) continue;
    obstacles.emplace_back(x(i), x(i) + width(i), y(i), y(i) + height(i));
  }
  std::vector<Rectangle> ret;
  for (Rectangle row : rows) {
    bpl::polygon_90_set_data<int> row_set;
    row_set.insert(
        bpl::rectangle_data<int>(row.minX, row.minY, row.maxX, row.maxY));
    for (Rectangle r : obstacles) {
      row_set.insert(bpl::rectangle_data<int>(r.minX, r.minY, r.maxX, r.maxY),
                     true);
    }

    std::vector<bpl::rectangle_data<int> > diff;
    bpl::get_rectangles(diff, row_set);
    for (auto r : diff) {
      Rectangle newRow(bpl::xl(r), bpl::xh(r), bpl::yl(r), bpl::yh(r));
      // Filter out partially covered rows
      if (newRow.height() == row.height()) {
        ret.push_back(newRow);
      }
    }
  }

  return ret;
}

void Circuit::check() const {
  assert(cellWidths.size() == nbCells());
  assert(cellHeights.size() == nbCells());
  assert(cellFixed.size() == nbCells());
  assert(cellX.size() == nbCells());
  assert(cellY.size() == nbCells());
  assert(cellOrientation.size() == nbCells());
  assert(!netLimits.empty());
  assert(netLimits.front() == 0);
  assert(pinCells.size() == nbPins());
  assert(pinXOffsets.size() == nbPins());
  assert(pinYOffsets.size() == nbPins());
}

void place(Circuit &circuit, int effort) {
  std::cout << "Placing circuit with " << circuit.nbCells() << " cells, "
            << circuit.nbNets() << " nets and " << circuit.nbPins() << " pins."
            << std::endl;
  GlobalPlacer::place(circuit, effort);
  std::cout << "Wirelength after global placement: " << circuit.hpwl()
            << std::endl;
  DetailedPlacer::place(circuit, effort);
  std::cout << "Wirelength after detailed placement: " << circuit.hpwl()
            << std::endl;
}

extern "C" {
int place_ispd(int nb_cells, int nb_nets, int *cell_widths, int *cell_heights,
                char *cell_fixed, int *net_limits, int *pin_cells,
                int *pin_x_offsets, int *pin_y_offsets, int *cell_x,
                int *cell_y, int *cell_orientation, int nb_rows, int *row_min_x,
                int *row_max_x, int *row_min_y, int *row_max_y, int effort) {
  
  Circuit circuit(nb_cells);
    
  try {
    
    circuit.setCellWidths(std::vector<int>(cell_widths, cell_widths + nb_cells));
    circuit.setCellHeights(
        std::vector<int>(cell_heights, cell_heights + nb_cells));
    circuit.setCellFixed(std::vector<char>(cell_fixed, cell_fixed + nb_cells));
    int nb_pins = net_limits[nb_nets];
    circuit.setNets(std::vector<int>(net_limits, net_limits + nb_nets + 1),
                    std::vector<int>(pin_cells, pin_cells + nb_pins),
                    std::vector<int>(pin_x_offsets, pin_x_offsets + nb_pins),
                    std::vector<int>(pin_y_offsets, pin_y_offsets + nb_pins));
    circuit.setCellX(std::vector<int>(cell_x, cell_x + nb_cells));
    circuit.setCellY(std::vector<int>(cell_y, cell_y + nb_cells));
    std::vector<CellOrientation> orient;
    for (int i = 0; i < nb_cells; ++i) {
      orient.push_back(static_cast<CellOrientation>(cell_orientation[i]));
    }
    circuit.setOrientation(orient);
    std::vector<Rectangle> rows;
    for (int i = 0; i < nb_rows; ++i) {
      rows.emplace_back(row_min_x[i], row_max_x[i], row_min_y[i], row_max_y[i]);
    }
    circuit.setRows(rows);
    
  } catch (const std::exception& e) {
    std::cout << "initialization terminated with exception: " << e.what() << std::endl;
    return -1;
  }
  
  try {
    
    place(circuit, effort);
  
  } catch (const std::exception& e) {
    std::cout << "placement terminated with exception: " << e.what() << std::endl;
    return -1;
  }
  
  for (int i = 0; i < nb_cells; ++i) {
    cell_x[i] = circuit.cellX[i];
    cell_y[i] = circuit.cellY[i];
  }
  
  // Return with success.
  return 0;
}
}
}
