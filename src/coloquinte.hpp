#pragma once

#include <vector>


struct Rectangle {
    Rectangle() : minX(0), maxX(0), minY(0), maxY(0) {}
    Rectangle(int minX, int maxX, int minY, int maxY) : minX(minX), maxX(maxX), minY(minY), maxY(maxY) {}
    int minX;
    int maxX;
    int minY;
    int maxY;

    int width() const { return maxX - minX; }
    int height() const { return maxY - minY; }
    long long area() const { return (long long) width() * (long long) height(); }
};

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
        char *cell_flip_y,
        int min_x, int max_x, int min_y, int max_y
    );
    void benchmark_quadratic_models(
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
        char *cell_flip_y,
        int model_type,
        int nb_steps,
        float epsilon,
        float relaxation
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
    Rectangle placementArea;

    int nbCells() const {
        return cellWidths.size();
    }

    int nbNets() const {
        return netLimits.size() - 1;
    }

    int nbPins() const {
        return netLimits.back();
    }

    bool isFixed(int cell) {
        return cellFixed[cell];
    }

    long long getArea(int cell) {
        return static_cast<long long>(cellWidths[cell]) * static_cast<long long>(cellHeights[cell]);
    }

    static Circuit createIspd(
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
        char *cell_flip_y,
        int min_x, int max_x, int min_y, int max_y
    );

    void check() const;
};

