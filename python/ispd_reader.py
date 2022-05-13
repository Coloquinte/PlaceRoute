
from . import circuit

import gzip
import lzma
import os
import numpy as np


def read_aux(filename):
    if os.path.isdir(filename):
        dir_list = os.listdir(filename)
        all_files = [f for f in dir_list if f.endswith(".aux")]
        if len(all_files) != 1:
            raise RuntimeError(f"There should be one file ending with .aux, got {len(all_files)}")
        filename = os.path.join(filename, all_files[0])
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
    if  len(node_files) != 1:
        raise RuntimeError(f"There should be a .nodes file in .aux")
    if  len(net_files) != 1:
        raise RuntimeError(f"There should be a .nets file in .aux")
    if  len(pl_files) != 1:
        raise RuntimeError(f"There should be a .pl file in .aux")
    if  len(scl_files) != 1:
        raise RuntimeError(f"There should be a .scl file in .aux")
    return (os.path.join(dirname, node_files[0]),
            os.path.join(dirname, net_files[0]),
            os.path.join(dirname, pl_files[0]),
            os.path.join(dirname, scl_files[0]),
            )


def open_file(name):
    if os.path.exists(name):
        return open(name)
    elif os.path.exists(name + ".gz"):
        return gzip.open(name + ".gz", mode="rt")
    elif os.path.exists(name + ".xz"):
        return lzma.open(name + ".xz", mode="rt")
    elif os.path.exists(name + ".lzma"):
        return lzma.open(name + ".lzma", mode="rt")
    else:
        raise RuntimeError(f"Could not find file {name}")


def parse_num_line(line):
    tokens = line.split(':')
    if len(tokens) != 2:
        raise RuntimeError(f"Couldn't interpret <{line}> as <Key : Value>")
    return int(tokens[1].strip())


def read_nodes(filename):
    nb_nodes = None
    nb_terminals = None
    nodes = []
    with open_file(filename) as f:
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
                nb_nodes = parse_num_line(line)
                continue
            if line.startswith("NumTerminals"):
                assert nb_terminals is None
                nb_terminals = parse_num_line(line)
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
    widths = np.array([n[1] for n in nodes], dtype=np.int32)
    heights = np.array([n[2] for n in nodes], dtype=np.int32)
    fixed = np.array([n[3] for n in nodes], dtype=np.bool)
    return names, widths, heights, fixed


def read_nets(filename, cell_names, cell_widths, cell_heights):
    name_dir = dict((name, i) for i, name in enumerate(cell_names))
    nb_nets = None
    nb_pins = None
    net_degree = None
    nets = []
    with open_file(filename) as f:
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
                nb_nets = parse_num_line(line)
                continue
            if line.startswith("NumPins"):
                assert nb_pins is None
                nb_pins = parse_num_line(line)
                continue
            line = line.replace(":", " ")
            if line.startswith("NetDegree"):
                vals = line.split()
                assert len(vals) == 3
                net_degree = int(vals[-2])
                name = vals[-1]
                nets.append((name, net_degree, []))
                continue
            vals = line.split()
            assert len(vals) == 4
            cell, direction, x, y = vals
            x = float(x)
            y = float(y)
            assert cell in name_dir
            nets[-1][2].append((name_dir[cell], direction, x, y))
    total_pins = 0
    for name, net_degree, pins in nets:
        if net_degree != len(pins):
            raise RuntimeError(f"Net degree for {name} is {len(pins)}; expected {net_degree}")
        total_pins += net_degree
    if nb_nets is not None:
        assert len(nets) == nb_nets
    if nb_pins is not None:
        assert total_pins == nb_pins
    # Sort so that same-size nets are contiguous
    cell_x_offset = 0.5 * cell_widths
    cell_y_offset = 0.5 * cell_heights
    nets.sort(key=lambda net: len(net[-1]))
    names = [n[0] for n in nets]
    degrees = [n[1] for n in nets]
    net_limits = [0]
    for d in degrees:
        net_limits.append(net_limits[-1] + d)
    pin_cell = []
    pin_x = []
    pin_y = []
    for n in nets:
        for cell, _, x, y in n[2]:
            pin_cell.append(cell)
            pin_x.append(int(round(cell_x_offset[cell] + x)))
            pin_y.append(int(round(cell_y_offset[cell] + y)))
    return (names,
            np.array(net_limits, dtype=np.int32),
            np.array(pin_cell, dtype=np.int32),
            np.array(pin_x, dtype=np.int32),
            np.array(pin_y, dtype=np.int32)
           )


def read_place(filename, cell_names):
    name_dir = dict((name, i) for i, name in enumerate(cell_names))
    cell_x = np.zeros(len(cell_names), dtype=np.int32)
    cell_y = np.zeros(len(cell_names), dtype=np.int32)
    flip_x = np.zeros(len(cell_names), dtype=np.bool)
    flip_y = np.zeros(len(cell_names), dtype=np.bool)
    with open_file(filename) as f:
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
            assert orient == "N"
    return cell_x, cell_y, flip_x, flip_y


def read_rows(filename):
    nb_rows = None
    rows = []
    with open_file(filename) as f:
        lines = [l.strip() for l in f]
        for line in lines:
            if line.startswith("NumRows"):
                assert nb_rows is None
                nb_rows = parse_num_line(line)
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
                if desc[i-1] == "Coordinate":
                    min_y = int(desc[i])
                if desc[i-1] == "SubrowOrigin":
                    min_x = int(desc[i])
                if desc[i-1] == "NumSites":
                    width = int(desc[i])
                if desc[i-1] == "Height":
                    height = int(desc[i])

            assert min_x is not None
            assert min_y is not None
            assert width is not None
            assert height is not None
            rows.append((min_x, min_y, min_x + width, min_y + height))

        min_x = [row[0] for row in rows]
        min_y = [row[1] for row in rows]
        max_x = [row[2] for row in rows]
        max_y = [row[3] for row in rows]
        heights = [row[3] - row[1] for row in rows]

        # Various checks, since our placement is simplified

        # All rows have same height
        heights = list(sorted(set(heights)))
        if len(set(heights)) != 1:
            raise RuntimeError(f"Only one row height is supported, got {', '.join([str(i) for i in heights])}")
        row_height = heights[0]

        # All rows have same min x, otherwise we take the min
        min_x = list(sorted(set(min_x)))
        if len(min_x) != 1:
            print(f"Only one row origin is supported, got {', '.join([str(i) for i in min_x])}")

        # All rows have same max x, otherwise we take the min
        max_x = list(sorted(set(max_x)))
        if len(max_x) != 1:
            print(f"Only one row end is supported, got {', '.join([str(i) for i in max_x])}")

        # Rows are contiguous
        min_y.sort()
        for y, nxt_y in zip(min_y[:-1], min_y[1:]):
            if nxt_y - y  != row_height:
                raise RuntimeError(f"Hole between rows at coordinates {y} and {nxt_y}")

        return (min(min_x), min(min_y), max(max_x), max(max_y)), row_height


def read_ispd(filename):
    node_filename, net_filename, pl_filename, scl_filename = read_aux(filename)
    cell_names, cell_widths, cell_heights, cell_fixed = read_nodes(node_filename)
    net_names, net_limits, pin_cell, pin_x, pin_y = read_nets(net_filename, cell_names, cell_widths, cell_heights)
    cell_x, cell_y, flip_x, flip_y = read_place(pl_filename, cell_names)
    pl_area, row_height = read_rows(scl_filename)
    ret = circuit.Circuit()
    ret._nb_cells = len(cell_names)
    ret._nb_nets = len(net_names)
    ret._cell_name = cell_names
    ret._cell_width = cell_widths
    ret._cell_height = cell_heights
    ret._cell_fixed = cell_fixed
    ret._net_name = net_names
    ret._net_limits = net_limits
    ret._pin_cells = pin_cell
    ret._pin_x = pin_x
    ret._pin_y = pin_y
    ret._cell_x = cell_x
    ret._cell_y = cell_y
    ret._cell_flip_x = flip_x
    ret._cell_flip_y = flip_y
    ret._pl_area = pl_area
    ret._row_height = row_height
    ret.check()
    return ret

