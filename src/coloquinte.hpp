#pragma once

#include <vector>

extern "C" {
    void place_ispd(
        int nb_cells,
        int nb_nets,
        int *cell_widths,
        int *cell_heights,
        char *cell_fixed,
        int *net_limits,
        int *pin_cells,
        int *pin_x_offsets,
        int *pin_y_offsets,
        int *cell_x,
        int *cell_y,
        char *cell_flip_x,
        char *cell_flip_y
    );
}


struct Circuit {
    /**
     * Representation of a circuit
     */

    std::vector<int> cellWidths;
    std::vector<int> cellHeights;
    std::vector<char> cellFixed;
    std::vector<int> netLimits;
    std::vector<int> pinCells;
    std::vector<int> pinXOffsets;
    std::vector<int> pinYOffsets;
    std::vector<int> cellX;
    std::vector<int> cellY;
    std::vector<char> cellFlipX;
    std::vector<char> cellFlipY;

    int nbCells() const {
        return cellWidths.size();
    }

    int nbNets() const {
        return netLimits.size() - 1;
    }

    int nbPins() const {
        return netLimits.back();
    }

    static Circuit create_ispd(
        int nb_cells,
        int nb_nets,
        int *cell_widths,
        int *cell_heights,
        char *cell_fixed,
        int *net_limits,
        int *pin_cells,
        int *pin_x_offsets,
        int *pin_y_offsets,
        int *cell_x,
        int *cell_y,
        char *cell_flip_x,
        char *cell_flip_y
    );

    void check() const;
};

