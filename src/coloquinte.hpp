#pragma once

#include "circuit.hpp"

extern "C" {
void place_ispd(int nb_cells, int nb_nets, int *cell_widths, int *cell_heights,
                char *cell_fixed, int *net_limits, int *pin_cells,
                int *pin_x_offsets, int *pin_y_offsets, int *cell_x,
                int *cell_y, int *cell_orientation, int nb_rows, int *row_min_x,
                int *row_max_x, int *row_min_y, int *row_max_y);
void benchmark_quadratic_models(int nb_cells, int nb_nets, int *cell_widths,
                                int *cell_heights, char *cell_fixed,
                                int *net_limits, int *pin_cells,
                                int *pin_x_offsets, int *pin_y_offsets,
                                int *cell_x, int *cell_y, int *cell_orientation,
                                int model_type, int nb_steps, float epsilon,
                                float relaxation);
}
