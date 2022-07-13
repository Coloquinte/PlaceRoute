#pragma once

#include <vector>


/**
 * 2D point
 */
struct Point {
    int x;
    int y;
};

/**
 * 2D rectangle
 */
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
    bool intersects(Rectangle o) const {
        return minX < o.maxX
            && o.minX < maxX
            && minY < o.maxY
            && o.minY < maxY;
    }
};

/**
 * Cost model to use when doing legalization
 */
enum LegalizationModel { L1, L2, LInf, L2Squared };
float norm(float x, float y, LegalizationModel leg);
long long norm(int x, int y, LegalizationModel leg);

/**
 * Representation of a flat circuit from an ISPD benchmark
 */

struct Circuit {
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

    long long getArea(int cell) const {
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


