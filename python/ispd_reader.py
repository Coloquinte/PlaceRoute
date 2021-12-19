
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
    if  len(node_files) != 1 and len(net_files) != 1 and len(pl_files) != 1:
        raise RuntimeError(f"There should be one file of each type in .aux")
    return (os.path.join(dirname, node_files[0]),
            os.path.join(dirname, net_files[0]),
            os.path.join(dirname, pl_files[0]))


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


def read_ispd(filename):
    node_filename, net_filename, pl_filename = read_aux(filename)
    cell_names, cell_widths, cell_heights, cell_fixed = read_nodes(node_filename)
    net_names, net_limits, pin_cell, pin_x, pin_y = read_nets(net_filename, cell_names, cell_widths, cell_heights)
    cell_x, cell_y, flip_x, flip_y = read_place(pl_filename, cell_names)
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
    ret.check()
    return ret

