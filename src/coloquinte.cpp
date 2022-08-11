
#include "coloquinte.hpp"

#include <iostream>

#include "place_detailed/place_detailed.hpp"
#include "place_global/density_legalizer.hpp"
#include "place_global/net_model.hpp"
#include "place_global/place_global.hpp"

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
void place_ispd(int nb_cells, int nb_nets, int *cell_widths, int *cell_heights,
                char *cell_fixed, int *net_limits, int *pin_cells,
                int *pin_x_offsets, int *pin_y_offsets, int *cell_x,
                int *cell_y, int *cell_orientation, int nb_rows, int *row_min_x,
                int *row_max_x, int *row_min_y, int *row_max_y, int effort) {
  Circuit circuit(nb_cells);
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
  place(circuit, effort);
  for (int i = 0; i < nb_cells; ++i) {
    cell_x[i] = circuit.cellX[i];
    cell_y[i] = circuit.cellY[i];
  }
}
}
