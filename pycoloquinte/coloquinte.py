"""
Coloquinte VLSI placer
"""


import gzip
import lzma
import os

import coloquinte_pybind
from coloquinte_pybind import (CellOrientation, DetailedPlacerParameters,
                               GlobalPlacerParameters, LegalizationModel,
                               Rectangle)


def _open_file(name, write=False):
    """
    Open the file with the appropriate decompression method. In read mode, search for compressed versions if the exact name does not exist.
    """
    mode = "wt" if write else "rt"
    if name.endswith(".gz"):
        return gzip.open(name, mode=mode)
    elif name.endswith(".xz") or name.endswith(".lzma"):
        return lzma.open(name, mode=mode)
    elif write:
        return open(name, mode=mode)
    elif os.path.exists(name):
        return open(name, mode=mode)
    elif os.path.exists(name + ".gz"):
        return gzip.open(name + ".gz", mode=mode)
    elif os.path.exists(name + ".xz"):
        return lzma.open(name + ".xz", mode=mode)
    elif os.path.exists(name + ".lzma"):
        return lzma.open(name + ".lzma", mode=mode)
    else:
        raise RuntimeError(f"Could not find file {name}")


def _read_aux(filename):
    if os.path.isdir(filename):
        dir_list = os.listdir(filename)
        all_files = [f for f in dir_list if f.endswith(".aux")]
        default_name = os.path.basename(filename) + ".aux"
        if len(all_files) == 1:
            filename = os.path.join(filename, all_files[0])
        elif len(all_files) > 1 and default_name in all_files:
            filename = os.path.join(filename, default_name)
        else:
            raise RuntimeError(
                f"There should be one file ending with .aux, got {', '.join(all_files)}")
    elif not os.path.exists(filename):
        filename = filename + ".aux"
    dirname = os.path.dirname(filename)
    files = []
    with open(filename) as f:
        for line in f:
            for name in line.split():
                files.append(name)
    node_files = [n for n in files if n.endswith(".nodes")]
    net_files = [n for n in files if n.endswith(".nets")]
    pl_files = [n for n in files if n.endswith(".pl")]
    scl_files = [n for n in files if n.endswith(".scl")]
    if len(node_files) != 1:
        raise RuntimeError("There should be a .nodes file in .aux")
    if len(net_files) != 1:
        raise RuntimeError("There should be a .nets file in .aux")
    if len(pl_files) != 1:
        raise RuntimeError("There should be a .pl file in .aux")
    if len(scl_files) != 1:
        raise RuntimeError("There should be a .scl file in .aux")
    return (os.path.join(dirname, node_files[0]),
            os.path.join(dirname, net_files[0]),
            os.path.join(dirname, pl_files[0]),
            os.path.join(dirname, scl_files[0]),
            )


def _parse_num_line(line):
    tokens = line.split(':')
    if len(tokens) != 2:
        raise RuntimeError(f"Couldn't interpret <{line}> as <Key : Value>")
    return int(tokens[1].strip())


def _read_nodes(filename):
    nb_nodes = None
    nb_terminals = None
    nodes = []
    with _open_file(filename) as f:
        first_line_found = False
        for line in f:
            line = line.strip()
            if len(line) == 0:
                continue
            if line.startswith("#"):
                continue
            if line.startswith("UCLA") and not first_line_found:
                first_line_found = True
                continue
            if line.startswith("NumNodes"):
                assert nb_nodes is None
                nb_nodes = _parse_num_line(line)
                continue
            if line.startswith("NumTerminals"):
                assert nb_terminals is None
                nb_terminals = _parse_num_line(line)
                continue
            vals = line.split()
            assert 3 <= len(vals) <= 4
            name, width, height = vals[:3]
            width = int(width)
            height = int(height)
            fixed = False
            if len(vals) == 4:
                assert vals[3] == "terminal"
                fixed = True
            nodes.append((name, width, height, fixed))
    if nb_nodes is not None:
        assert len(nodes) == nb_nodes
    if nb_terminals is not None:
        assert len([n for n in nodes if n[3]]) == nb_terminals
    names = [n[0] for n in nodes]
    widths = [n[1] for n in nodes]
    heights = [n[2] for n in nodes]
    fixed = [n[3] for n in nodes]
    return names, widths, heights, fixed


def _read_nets(filename, cell_names, cell_widths, cell_heights, sort_entries=False):
    name_dir = dict((name, i) for i, name in enumerate(cell_names))
    nb_nets = None
    nb_pins = None
    net_degree = None
    nets = []
    with _open_file(filename) as f:
        first_line_found = False
        for line in f:
            line = line.strip()
            if len(line) == 0:
                continue
            if line.startswith("#"):
                continue
            if line.startswith("UCLA") and not first_line_found:
                first_line_found = True
                continue
            if line.startswith("NumNets"):
                assert nb_nets is None
                nb_nets = _parse_num_line(line)
                continue
            if line.startswith("NumPins"):
                assert nb_pins is None
                nb_pins = _parse_num_line(line)
                continue
            line = line.replace(":", " ")
            if line.startswith("NetDegree"):
                vals = line.split()
                assert 2 <= len(vals) <= 3
                net_degree = int(vals[1])
                if len(vals) == 3:
                    name = vals[2]
                else:
                    name = f"n{len(nets)}"
                nets.append((name, net_degree, []))
                continue
            vals = line.split()
            assert len(vals) == 4 or len(vals) == 2
            if len(vals) == 4:
                cell, direction, x, y = vals
                x = float(x)
                y = float(y)
            else:
                cell, direction = vals
                x = 0.0
                y = 0.0
            assert cell in name_dir
            nets[-1][2].append((name_dir[cell], direction, x, y))
    total_pins = 0
    for name, net_degree, pins in nets:
        if net_degree != len(pins):
            raise RuntimeError(
                f"Net degree for {name} is {len(pins)}; expected {net_degree}")
        total_pins += net_degree
    if nb_nets is not None:
        assert len(nets) == nb_nets
    if nb_pins is not None:
        assert total_pins == nb_pins
    if sort_entries:
        # Sort so that same-size nets are contiguous
        nets.sort(key=lambda net: len(net[-1]))
    cell_x_offset = [0.5 * c for c in cell_widths]
    cell_y_offset = [0.5 * c for c in cell_heights]
    ret = []
    for name, _, pins in nets:
        cells = []
        pin_x = []
        pin_y = []
        for cell, _, x, y in pins:
            cells.append(cell)
            pin_x.append(int(round(cell_x_offset[cell] + x)))
            pin_y.append(int(round(cell_y_offset[cell] + y)))
        ret.append((name, cells, pin_x, pin_y))
    return ret


def _read_place(filename, cell_names):
    name_dir = dict((name, i) for i, name in enumerate(cell_names))
    cell_x = [0 for i in cell_names]
    cell_y = [0 for i in cell_names]
    cell_orient = [None for i in cell_names]

    with _open_file(filename) as f:
        first_line_found = False
        for line in f:
            line = line.strip()
            if len(line) == 0:
                continue
            if line.startswith("#"):
                continue
            if line.startswith("UCLA") and not first_line_found:
                first_line_found = True
                continue
            line = line.replace(":", " ")
            vals = line.split()
            assert len(vals) >= 3
            cell, x, y, orient = vals[:4]
            assert cell in name_dir
            cell_ind = name_dir[cell]
            cell_x[cell_ind] = int(x)
            cell_y[cell_ind] = int(y)
            if orient in CellOrientation.__members__:
                cell_orient[cell_ind] = CellOrientation.__members__[orient]
            else:
                raise RuntimeError(f"Unknown orientation encountered {orient}")
    return cell_x, cell_y, cell_orient


def _read_rows(filename):
    nb_rows = None
    rows = []
    with _open_file(filename) as f:
        lines = [l.strip() for l in f]
        for line in lines:
            if line.startswith("NumRows"):
                assert nb_rows is None
                nb_rows = _parse_num_line(line)
        row_descs = []
        in_row = False
        for line in lines:
            if line.startswith("CoreRow"):
                row_descs.append([])
                in_row = True
            elif line.startswith("End"):
                in_row = False
            elif in_row:
                row_descs[-1].extend(line.replace(":", " ").split())
        for desc in row_descs:
            min_x = None
            min_y = None
            width = None
            height = None
            for i in range(1, len(desc)):
                if desc[i-1].lower() == "coordinate":
                    min_y = int(desc[i])
                if desc[i-1].lower() == "subroworigin":
                    min_x = int(desc[i])
                if desc[i-1].lower() == "numsites":
                    width = int(desc[i])
                if desc[i-1].lower() == "height":
                    height = int(desc[i])

            assert min_x is not None
            assert min_y is not None
            assert width is not None
            assert height is not None
            rows.append((min_x, min_x + width, min_y, min_y + height))
        return rows


def _rows_to_area(rows):
    min_x = [row[0] for row in rows]
    min_y = [row[1] for row in rows]
    max_x = [row[2] for row in rows]
    max_y = [row[3] for row in rows]

    # Various checks, since our placement is simplified

    # All rows have same min x, otherwise we take the min
    min_x = list(sorted(set(min_x)))
    if len(min_x) != 1:
        print(
            f"Only one row origin is supported, got {', '.join([str(i) for i in min_x])}")

    # All rows have same max x, otherwise we take the min
    max_x = list(sorted(set(max_x)))
    if len(max_x) != 1:
        print(
            f"Only one row end is supported, got {', '.join([str(i) for i in max_x])}")

    # Rows are contiguous
    row_y = [p for p in zip(min_y, max_y)]
    row_y.sort()
    for i in range(len(row_y)-1):
        y_b = row_y[i+1][0]
        y_e = row_y[i][1]
        if y_b != y_e:
            print(f"Hole between rows at coordinates {y_b} and {y_e}")

    return (min(min_x), min(max_x), min(min_y), max(max_y))


def _rows_to_height(rows):
    heights = [row[3] - row[1] for row in rows]
    # All rows have same height
    heights = list(sorted(set(heights)))
    if len(set(heights)) != 1:
        raise RuntimeError(
            f"Only one row height is supported, got {', '.join([str(i) for i in heights])}")
    return heights[0]


class Circuit(coloquinte_pybind.Circuit):
    def __init__(self, nb_cells):
        super(Circuit, self).__init__(nb_cells)
        self._filename = None
        self._cell_name = None
        self._net_name = None

    @staticmethod
    def read_ispd(filename):
        """
        Read an ISPD benchmark from its .aux file
        """
        node_filename, net_filename, pl_filename, scl_filename = _read_aux(
            filename)
        cell_names, cell_widths, cell_heights, cell_fixed = _read_nodes(
            node_filename)
        nets = _read_nets(
            net_filename, cell_names, cell_widths, cell_heights)
        cell_x, cell_y, cell_orient = _read_place(pl_filename, cell_names)
        rows = _read_rows(scl_filename)

        ret = Circuit(len(cell_names))
        ret._filename = filename

        # Setup cell properties
        ret._cell_name = cell_names
        ret.cell_width = cell_widths
        ret.cell_height = cell_heights
        ret.cell_fixed = cell_fixed

        # Setup nets and pins
        ret._net_name = []
        for name, cells, pins_x, pins_y in nets:
            ret._net_name.append(name)
            ret.add_net(cells, pins_x, pins_y)

        # Setup initial cell placement
        ret.cell_x = cell_x
        ret.cell_y = cell_y
        ret._cell_orientation = cell_orient

        # Setup rows
        ret.rows = [coloquinte_pybind.Rectangle(*row) for row in rows]
        ret.check()
        return ret

    def write_placement(self, filename):
        """
        Write the placement result in ISPD file format
        """
        if filename is None:
            if self._filename is None:
                raise RuntimeError("No filename to export placement to")
            filename = os.path.splitext(self._filename)[0] + ".sol.pl"
        with _open_file(filename, True) as f:
            print("UCLA pl 1.0", file=f)
            print("# Created by Coloquinte", file=f)
            print("# https://github.com/Coloquinte/PlaceRoute", file=f)
            print("", file=f)
            cell_x = self.cell_x
            cell_y = self.cell_y
            cell_orientation = self.cell_orientation
            cell_fixed = self.cell_fixed
            for i in range(self.nb_cells):
                name = self._cell_name[i]
                x = cell_x[i]
                y = cell_y[i]
                orient = cell_orientation[i].name
                if cell_fixed[i]:
                    orient += " /FIXED"
                print(f"{name}\t{x}\t{y}\t: {orient}", file=f)


def main():
    """
    Run the whole placement algorithm from the command line
    """
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("instance", help="Benchmark instance")
    parser.add_argument("--solution", help="Placement result")
    parser.add_argument("--effort", help="Placement effort",
                        type=int, default=3)
    args = parser.parse_args()

    circuit = Circuit.read_ispd(args.instance)
    circuit.place(args.effort)
    circuit.write_placement(args.solution)


__all__ = [
    "Circuit",
    "GlobalPlacerParameters",
    "DetailedPlacerParameters",
    "Rectangle",
    "LegalizationModel",
    "CellOrientation",
    "main",
]
