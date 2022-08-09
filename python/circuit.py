
import gzip
import lzma
import os
import numpy as np


orient_map = {
    "N": 0, "S": 1, "W": 2, "E": 3,
    "FN": 4, "FS": 5, "FW": 6, "FE": 7
}


def open_file(name, write=False):
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


def read_aux(filename):
    if os.path.isdir(filename):
        dir_list = os.listdir(filename)
        all_files = [f for f in dir_list if f.endswith(".aux")]
        if len(all_files) != 1:
            raise RuntimeError(
                f"There should be one file ending with .aux, got {len(all_files)}")
        filename = os.path.join(filename, all_files[0])
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
        raise RuntimeError(f"There should be a .nodes file in .aux")
    if len(net_files) != 1:
        raise RuntimeError(f"There should be a .nets file in .aux")
    if len(pl_files) != 1:
        raise RuntimeError(f"There should be a .pl file in .aux")
    if len(scl_files) != 1:
        raise RuntimeError(f"There should be a .scl file in .aux")
    return (os.path.join(dirname, node_files[0]),
            os.path.join(dirname, net_files[0]),
            os.path.join(dirname, pl_files[0]),
            os.path.join(dirname, scl_files[0]),
            )


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


def read_nets(filename, cell_names, cell_widths, cell_heights, sort_entries=False):
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
    cell_x_offset = 0.5 * cell_widths
    cell_y_offset = 0.5 * cell_heights
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
    cell_orient = np.zeros(len(cell_names), dtype=np.int32)

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
            if orient in orient_map:
                cell_orient[cell_ind] = orient_map[orient]
            else:
                print(f"Unknown orientation encountered {orient}")
    return cell_x, cell_y, cell_orient


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
            rows.append((min_x, min_y, min_x + width, min_y + height))
        return rows


def rows_to_numpy(rows):
    min_x = np.array([row[0] for row in rows], dtype=np.int32)
    min_y = np.array([row[1] for row in rows], dtype=np.int32)
    max_x = np.array([row[2] for row in rows], dtype=np.int32)
    max_y = np.array([row[3] for row in rows], dtype=np.int32)
    return min_x, min_y, max_x, max_y


def rows_to_area(rows):
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


def rows_to_height(rows):
    heights = [row[3] - row[1] for row in rows]
    # All rows have same height
    heights = list(sorted(set(heights)))
    if len(set(heights)) != 1:
        raise RuntimeError(
            f"Only one row height is supported, got {', '.join([str(i) for i in heights])}")
    return heights[0]


class Circuit:
    def __init__(self):
        self._nb_cells = None
        self._nb_nets = None
        self._cell_name = None
        self._cell_height = None
        self._cell_width = None
        self._cell_fixed = None
        self._cell_x = None
        self._cell_y = None
        self._cell_orient = None
        self._net_name = None
        self._net_limits = None
        self._pin_cells = None
        self._pin_x = None
        self._pin_y = None
        self._nb_rows = None
        self._row_min_x = None
        self._row_max_x = None
        self._row_min_y = None
        self._row_max_y = None

    @property
    def nb_cells(self):
        return self._nb_cells

    @property
    def nb_nets(self):
        return self._nb_nets

    @property
    def nb_rows(self):
        return self._nb_rows

    @property
    def cell_name(self):
        return self._cell_name

    @property
    def cell_height(self):
        return self._cell_height

    @property
    def cell_width(self):
        return self._cell_width

    @property
    def cell_fixed(self):
        return self._cell_fixed

    def net_pins_cell(self, i):
        return self._pin_cells[self._net_limits[i]:self._net_limits[i+1]]

    def net_pins_x(self, i):
        return self._pin_x[self._net_limits[i]:self._net_limits[i+1]]

    def net_pins_y(self, i):
        return self._pin_y[self._net_limits[i]:self._net_limits[i+1]]

    @property
    def cell_x(self):
        return self._cell_x

    @property
    def cell_y(self):
        return self._cell_y

    @property
    def cell_orient(self):
        return self._cell_orient

    def check(self):
        # No member is None
        assert self._nb_cells is not None
        assert self._nb_nets is not None
        assert self._nb_rows is not None
        assert self._cell_name is not None
        assert self._cell_height is not None
        assert self._cell_width is not None
        assert self._cell_fixed is not None
        assert self._cell_x is not None
        assert self._cell_y is not None
        assert self._cell_orient is not None
        assert self._net_name is not None
        assert self._net_limits is not None
        assert self._pin_cells is not None
        assert self._pin_x is not None
        assert self._pin_y is not None
        assert self._row_min_x is not None
        assert self._row_max_x is not None
        assert self._row_min_y is not None
        assert self._row_max_y is not None

        assert isinstance(self._cell_width, np.ndarray)
        assert isinstance(self._cell_height, np.ndarray)
        assert isinstance(self._cell_fixed, np.ndarray)
        assert isinstance(self._cell_x, np.ndarray)
        assert isinstance(self._cell_y, np.ndarray)
        assert isinstance(self._cell_orient, np.ndarray)
        assert isinstance(self._pin_cells, np.ndarray)
        assert isinstance(self._pin_x, np.ndarray)
        assert isinstance(self._pin_y, np.ndarray)
        assert isinstance(self._row_min_x, np.ndarray)
        assert isinstance(self._row_max_x, np.ndarray)
        assert isinstance(self._row_min_y, np.ndarray)
        assert isinstance(self._row_max_y, np.ndarray)

        assert self._cell_width.dtype == np.int32
        assert self._cell_height.dtype == np.int32
        assert self._cell_fixed.dtype == np.bool
        assert self._cell_x.dtype == np.int32
        assert self._cell_y.dtype == np.int32
        assert self._cell_orient.dtype == np.int32
        assert self._pin_cells.dtype == np.int32
        assert self._pin_x.dtype == np.int32
        assert self._pin_y.dtype == np.int32
        assert self._row_min_x.dtype == np.int32
        assert self._row_max_x.dtype == np.int32
        assert self._row_min_y.dtype == np.int32
        assert self._row_max_y.dtype == np.int32

        # Check the dimensions
        assert len(self._cell_name) == self.nb_cells
        assert len(self._cell_height) == self.nb_cells
        assert len(self._cell_width) == self.nb_cells
        assert len(self._cell_fixed) == self.nb_cells
        assert len(self._cell_x) == self.nb_cells
        assert len(self._cell_y) == self.nb_cells
        assert len(self._cell_orient) == self.nb_cells
        assert len(self._net_limits) == self.nb_nets + 1
        assert len(self._pin_cells) == len(self._pin_x)
        assert len(self._pin_cells) == len(self._pin_y)
        assert len(self._row_min_x) == self.nb_rows
        assert len(self._row_max_x) == self.nb_rows
        assert len(self._row_min_y) == self.nb_rows
        assert len(self._row_max_y) == self.nb_rows

        # Pin cells must be valid
        assert np.all(self._pin_cells >= 0)
        assert np.all(self._pin_cells < self.nb_cells)

        # Rows must have correct dimensions
        assert (self._row_min_x <= self._row_max_x).all()
        assert (self._row_min_y <= self._row_max_y).all()

    def place(self):
        """
        Call the C++ library to place the circuit
        """
        import ctypes
        dll = ctypes.CDLL("./build/libcoloquinte.so")
        c_bool_p = ctypes.POINTER(ctypes.c_char)
        c_int_p = ctypes.POINTER(ctypes.c_int)
        c_float_p = ctypes.POINTER(ctypes.c_float)
        dll.place_ispd(
            self.nb_cells,
            self.nb_nets,
            self._cell_width.ctypes.data_as(c_int_p),
            self._cell_height.ctypes.data_as(c_int_p),
            self._cell_fixed.ctypes.data_as(c_bool_p),
            self._net_limits.ctypes.data_as(c_int_p),
            self._pin_cells.ctypes.data_as(c_int_p),
            self._pin_x.ctypes.data_as(c_int_p),
            self._pin_y.ctypes.data_as(c_int_p),
            self._cell_x.ctypes.data_as(c_int_p),
            self._cell_y.ctypes.data_as(c_int_p),
            self._cell_orient.ctypes.data_as(c_int_p),
            self.nb_rows,
            self._row_min_x.ctypes.data_as(c_int_p),
            self._row_max_x.ctypes.data_as(c_int_p),
            self._row_min_y.ctypes.data_as(c_int_p),
            self._row_max_y.ctypes.data_as(c_int_p),
        )

    def benchmark_all_quadratic_models(self):
        for model_type in ["STAR", "BSTAR", "B2B"]:
            for epsilon in [0.1, 1.0, 10.0]:
                self.benchmark_quadratic_models(model_type, 20, epsilon, 0.0)

    def benchmark_quadratic_models(self, model_type, nb_steps, epsilon, relaxation):
        """
        Call the C++ library to o
        """
        model_types = {"STAR": 0, "BSTAR": 1, "B2B": 2}
        assert epsilon > 0.0
        assert 0.0 <= relaxation <= 1.0
        assert nb_steps > 0
        assert model_type in model_types.keys()

        import ctypes
        dll = ctypes.CDLL("./build/libcoloquinte.so")
        c_bool_p = ctypes.POINTER(ctypes.c_char)
        c_int_p = ctypes.POINTER(ctypes.c_int)
        c_float_p = ctypes.POINTER(ctypes.c_float)
        dll.benchmark_quadratic_models(
            self.nb_cells,
            self.nb_nets,
            self._cell_width.ctypes.data_as(c_int_p),
            self._cell_height.ctypes.data_as(c_int_p),
            self._cell_fixed.ctypes.data_as(c_bool_p),
            self._net_limits.ctypes.data_as(c_int_p),
            self._pin_cells.ctypes.data_as(c_int_p),
            self._pin_x.ctypes.data_as(c_int_p),
            self._pin_y.ctypes.data_as(c_int_p),
            self._cell_x.ctypes.data_as(c_int_p),
            self._cell_y.ctypes.data_as(c_int_p),
            self._cell_orient.ctypes.data_as(c_int_p),
            model_types[model_type],
            nb_steps,
            ctypes.c_float(epsilon),
            ctypes.c_float(relaxation)
        )

    def partition(self):
        """
        Try and partition the circuit using KaHypar's interface
        """
        import kahypar

        hyperedge_indices = [int(i) for i in self._net_limits]
        hyperedges = [int(i) for i in self._pin_cells]

        node_weights = [1 for i in range(self.nb_cells)]
        edge_weights = [1 for i in range(self.nb_nets)]

        k = 2
        hypergraph = kahypar.Hypergraph(
            self.nb_cells, self.nb_nets, hyperedge_indices, hyperedges, k, edge_weights, node_weights)

        context = kahypar.Context()
        context.loadINIconfiguration(
            "/home/alnurn/Documents/kahypar/config/km1_kKaHyPar_eco_sea20.ini")

        context.setK(k)
        context.setEpsilon(0.03)

        kahypar.partition(hypergraph, context)

    def export_hgr(self, filename):
        """
        Export the circuit in hgr format for hypergraph partitioning tools.

        Node areas are not exported
        """
        edges = []
        for i in range(self.nb_nets):
            b = self._net_limits[i]
            e = self._net_limits[i+1]
            pins = [int(self._pin_cells[j]) for j in range(b, e)]
            # Uniquify the pins and ensure that there are at least two
            pins = list(sorted(set(pins)))
            if len(pins) < 2:
                continue
            edges.append(pins)
        self.export_hgr_data(filename, self.nb_cells, edges)

    def export_edge_hgr(self, filename):
        """
        Export the circuit in hgr format to partition the edges, instead of the nodes
        """
        edges = [[] for n in range(self.nb_cells)]
        net_ind = 0
        for i in range(self.nb_nets):
            b = self._net_limits[i]
            e = self._net_limits[i+1]
            pins = [int(self._pin_cells[j]) for j in range(b, e)]
            # Uniquify the pins and ensure that there are at least two
            pins = list(sorted(set(pins)))
            if len(pins) < 2:
                continue
            for p in pins:
                edges[p].append(net_ind)
            net_ind += 1
        self.export_hgr_data(filename, net_ind, edges)

    def export_hgr_data(self, filename, nb_cells, edges):
        with open(filename, "w") as f:
            print(f"{len(edges)} {nb_cells}", file=f)
            for edge in edges:
                pins = " ".join(str(p + 1) for p in edge)
                print(pins, file=f)

    @staticmethod
    def read_ispd(filename):
        """
        Read an ISPD benchmark from its .aux file
        """
        node_filename, net_filename, pl_filename, scl_filename = read_aux(
            filename)
        cell_names, cell_widths, cell_heights, cell_fixed = read_nodes(
            node_filename)
        net_names, net_limits, pin_cell, pin_x, pin_y = read_nets(
            net_filename, cell_names, cell_widths, cell_heights)
        cell_x, cell_y, cell_orient = read_place(pl_filename, cell_names)
        rows = read_rows(scl_filename)

        ret = Circuit()
        ret._nb_cells = len(cell_names)
        ret._nb_nets = len(net_names)

        # Setup cell properties
        ret._cell_name = cell_names
        ret._cell_width = cell_widths
        ret._cell_height = cell_heights
        ret._cell_fixed = cell_fixed

        # Setup nets and pins
        ret._net_name = net_names
        ret._net_limits = net_limits
        ret._pin_cells = pin_cell
        ret._pin_x = pin_x
        ret._pin_y = pin_y

        # Setup initial cell placement
        ret._cell_x = cell_x
        ret._cell_y = cell_y
        ret._cell_orient = cell_orient

        # Setup rows
        r_min_x, r_min_y, r_max_x, r_max_y = rows_to_numpy(rows)
        ret._nb_rows = len(rows)
        ret._row_min_x = r_min_x
        ret._row_max_x = r_max_x
        ret._row_min_y = r_min_y
        ret._row_max_y = r_max_y
        ret.check()
        return ret

    def write_pl(self, filename):
        with open_file(filename, True) as f:
            print("UCLA pl 1.0", file=f)
            print("# Created by Coloquinte", file=f)
            print("# https://github.com/Coloquinte/PlaceRoute", file=f)
            print("", file=f)
            for i in range(self.nb_cells):
                name = self.cell_name[i]
                x = self.cell_x[i]
                y = self.cell_y[i]
                orient = "N"
                for k, v in orient_map.items():
                    if self.cell_orient[i] == v:
                        orient = k
                if self.cell_fixed[i]:
                    orient += " /FIXED"
                print(f"{name}\t{x}\t{y}\t: {orient}", file=f)
